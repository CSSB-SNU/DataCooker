"""Ray-backed fan-out for the LMDB build/rebuild loops.

Why this exists
---------------
The fork ("multiprocessing") backend sizes concurrency with a *predictive* memory
model (``peak = n_jobs x E x max_bytes``). That model has blind spots -- a deep MSA
compresses to a few MB on disk but its dense depth x length materialisation is 1-2 GB
live -- so a build that fits on paper OOMs a node, and the fix was hand-tuning
``n_jobs`` / ``E`` / shard counts per DB. Measured: on the 22 k deepest a3m records,
fork had to drop to ``n_jobs=16`` (to not OOM) and ran at ~15 rec/s.

Ray replaces prediction with **measurement + admission control**: its node memory
monitor kills a task that would push the node past ``RAY_memory_usage_threshold`` and
*automatically re-runs it* (``max_retries``), so concurrency self-regulates to what the
node can actually hold -- no per-DB tuning. We keep DataCooker's identity (recipes,
schemas, CLI, shard-merge, single-writer drain) unchanged: only the within-node fan-out
changed. ``build_lmdb`` / ``rebuild_lmdb`` use this unconditionally -- there is no
backend choice.

Design
------
* **Stateless batched tasks**, not actors. A task that the memory monitor kills is
  retried by Ray for free; an actor kill would lose its in-flight batch and needs manual
  restart handling. Batching (``batch`` keys/task) amortises the per-task state install
  and source-handle open over many records, recovering the throughput actors would give.
* The driver ``ray.put``s the build/rebuild state **once**; Ray stores it in the plasma
  object store (shared memory in ``/dev/shm``) and every task reads it zero-copy -- the
  same "shared read-only state, no per-task pickling" property the fork path got from
  copy-on-write.
* Results are regrouped into ``chunk_size`` lists and yielded, so the caller's
  single-writer drain loop (``_parallel_write``) is byte-for-byte unchanged.
* **Native Ray runtime policy.** We do not override Ray's temp dir or hand-manage
  ``/dev/shm``: Ray already isolates each ``ray.init()`` in its own ``session_*``
  dir (so concurrent SLURM shard tasks on one node don't collide) and cleans it on a
  graceful ``ray.shutdown``. A SIGKILL leak is bounded (the plasma file is sized to
  ``object_store_memory`` and named per-session, never colliding with a fresh run)
  and left to ops (``ray stop`` / SLURM epilog) rather than library-side rmtree.
"""
from __future__ import annotations

import logging
import os
from collections.abc import Iterator, Sequence
from importlib import import_module
from typing import Any

logger = logging.getLogger(__name__)


def _resolve(qualname: str) -> Any:
    module_name, attr = qualname.rsplit(".", 1)
    return getattr(import_module(module_name), attr)


def require_ray() -> Any:
    """Return the imported ``ray`` module or raise a clear error."""
    try:
        return import_module("ray")
    except ModuleNotFoundError as exc:  # pragma: no cover - env-dependent
        msg = (
            "DataCooker's build/rebuild fan-out requires the 'ray' dependency. "
            "Install DataCooker with ray support (e.g. `pixi add --pypi ray`)."
        )
        raise ModuleNotFoundError(msg) from exc


def _run_batch(
    worker_qualname: str,
    installer_qualname: str,
    state: Any,
    keys: Sequence[Any],
) -> list[Any]:
    """Ray task body: install the shipped state once per worker, then process a batch.

    ``state`` arrives already de-referenced from the plasma store (Ray resolves an
    ObjectRef argument before the call). We install it into the worker module's global
    only when this worker has not seen it yet -- Ray reuses worker processes across many
    tasks, so the install (and the source-handle open inside the worker) is paid once per
    worker, not once per batch.
    """
    worker = _resolve(worker_qualname)
    module_name = worker_qualname.rsplit(".", 1)[0]
    module = import_module(module_name)
    if not getattr(module, "_RAY_STATE_INSTALLED", False):
        _resolve(installer_qualname)(state)
        module._RAY_STATE_INSTALLED = True  # noqa: SLF001
    return [worker(k) for k in keys]


def _chunked(seq: Sequence[Any], size: int) -> Iterator[Sequence[Any]]:
    for i in range(0, len(seq), size):
        yield seq[i : i + size]


def iter_ray_chunks(
    *,
    items: Sequence[Any],
    chunk_size: int,
    worker_qualname: str,
    installer_qualname: str,
    state: Any,
    log_label: str,
    n_actors: int = 0,
    batch: int = 50,
    inflight_per_worker: int = 3,
    memory_threshold: float = 0.9,
    max_task_retries: int = 3,
) -> Iterator[list[Any]]:
    """Fan the build/rebuild loop out over Ray, yielding ``chunk_size``-sized result lists.

    Yields the SAME shape the fork path yields (a list whose entries are each worker's
    return value -- a record tuple, a list of them for 1:N, or ``None``), so the caller's
    drain loop is unchanged. ``memory_threshold`` sets ``RAY_memory_usage_threshold`` (the
    node fraction above which the monitor kills+retries the newest task).
    """
    ray = require_ray()
    # The monitor must be armed before ray.init(); it reads the env at startup.
    os.environ.setdefault("RAY_memory_usage_threshold", str(memory_threshold))
    started_here = not ray.is_initialized()
    if started_here:
        # Native Ray policy: each ray.init() gets its own ``session_<ts>_<pid>``
        # dir under the default temp root (plasma in /dev/shm), so concurrent SLURM
        # tasks on one node never collide, and a graceful ``ray.shutdown`` removes
        # our session. We do NOT override the temp dir or hand-clean /dev/shm --
        # a SIGKILL leak is bounded (plasma is sized to object_store_memory and
        # named per-session, so it never collides with a fresh run) and is an ops
        # concern (``ray stop`` / epilog), not the library's.
        ray.init(
            ignore_reinit_error=True,
            logging_level=logging.WARNING,
            include_dashboard=False,
        )
    resources = ray.cluster_resources()
    n_cpus = int(resources.get("CPU", os.cpu_count() or 8))
    n_workers = n_actors or n_cpus
    total = len(items)
    logger.info(
        "%s via ray: cpus=%d workers=%d batch=%d keys=%d thresh=%.2f",
        log_label, n_cpus, n_workers, batch, total, memory_threshold,
    )

    state_ref = ray.put(state)
    remote = ray.remote(max_retries=max_task_retries, num_cpus=1)(_run_batch)

    batches = list(_chunked(items, max(1, batch)))
    it = iter(batches)
    try:
        n_prime = min(n_workers * max(1, inflight_per_worker), len(batches))
        pending: list[Any] = [
            remote.remote(worker_qualname, installer_qualname, state_ref, next(it))
            for _ in range(n_prime)
        ]

        buffer: list[Any] = []
        done = 0
        while pending:
            ready, pending = ray.wait(pending, num_returns=1, timeout=None)
            for ref in ready:
                results = ray.get(ref)
                buffer.extend(results)
                done += len(results)
                nxt = next(it, None)
                if nxt is not None:
                    pending.append(
                        remote.remote(
                            worker_qualname, installer_qualname, state_ref, nxt,
                        ),
                    )
            while len(buffer) >= chunk_size:
                yield buffer[:chunk_size]
                del buffer[:chunk_size]
            if done and done % (chunk_size * 5) < batch:
                logger.info("%s via ray: %d/%d", log_label, done, total)
        if buffer:
            yield buffer
    finally:
        # Only tear down the cluster if we started it, so a caller that
        # pre-initialised Ray keeps ownership. Runs on error or generator-close
        # too; the graceful shutdown is what triggers Ray's native session cleanup
        # (its /dev/shm plasma file included), so a failed drain doesn't leak.
        if started_here:
            ray.shutdown()
