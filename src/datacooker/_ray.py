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

Admission is **measured, not fixed**. Ray's own memory monitor is only a reactive
last-resort net (it kills + retries a task after it has already blown the threshold,
which oscillates instead of converging when the steady-state footprint is too big). So
this loop actively sizes the in-flight task count to *measured* node memory (AIMD:
grow while memory is comfortable, halve on pressure). Concurrency self-adjusts to what
the node can hold as record sizes vary -- a light schema-A build ramps to full width, a
heavy schema-B attach settles low -- with no per-DB tuning and no env knobs.

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
* **Isolated, fresh-local Ray head.** Each ``ray.init()`` pins a unique per-process
  ``_temp_dir`` and passes ``address="local"``, so a shard task always starts its own
  private instance instead of auto-discovering one via the shared
  ``<temp_dir>/ray_current_cluster`` address file. That file lives at the ``/tmp/ray``
  root (not under the per-session dir), so a stale entry from a dead earlier run would
  otherwise hijack a new init and hang it forever retrying an unreachable GCS. Plasma
  still lives in ``/dev/shm`` per session; a SIGKILL leak is bounded (sized to
  ``object_store_memory``, named per-session) and left to ops (``ray stop`` / SLURM
  epilog) rather than a library-side rmtree.
"""
from __future__ import annotations

import logging
import os
from collections.abc import Iterator, Sequence
from importlib import import_module
from typing import Any

import psutil

# Adaptive-admission watermarks (fraction of node memory actually in use). The controller
# grows concurrency while below _MEM_LOW and backs off above _MEM_HIGH, staying under
# Ray's own last-resort kill threshold (0.9) so the monitor almost never has to fire.
_MEM_HIGH = 0.80
_MEM_LOW = 0.60


def _mem_frac() -> float:
    """Fraction of node memory in use (accounts for reclaimable cache via `available`)."""
    return psutil.virtual_memory().percent / 100.0

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
    module: Any = import_module(module_name)
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
    batch: int = 50,
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
        # Force a FRESH, node-local Ray head -- never attach to a discovered one.
        #
        # ``ray.init()`` with no address auto-discovers an existing cluster via the
        # address file ``<temp_dir>/ray_current_cluster``. Ray isolates each session's
        # *plasma/logs* under ``session_<ts>_<pid>``, but that address file lives at the
        # ``/tmp/ray`` ROOT and is shared across runs (and, on a shared node, across
        # users). A stale file from an earlier run whose GCS is long dead then hijacks a
        # new ``ray.init()``, which hangs forever retrying an unreachable GCS
        # (``Failed to connect to GCS at <ip>:6379``) -- exactly what stranded SLURM
        # shard tasks on nodes carrying a leftover file. Pin a unique per-process
        # ``_temp_dir`` (so we never read someone else's address file) and pass
        # ``address="local"`` (so we never attach to a discovered cluster or a leaked
        # ``RAY_ADDRESS``); each shard task gets its own private, self-cleaning instance.
        token = os.environ.get("SLURM_JOB_ID") or str(os.getpid())
        temp_dir = os.path.join(  # noqa: PTH118
            os.environ.get("TMPDIR", "/tmp"),  # noqa: S108
            f"ray-{token}-{os.getpid()}",
        )
        ray.init(
            address="local",
            _temp_dir=temp_dir,
            ignore_reinit_error=True,
            logging_level=logging.WARNING,
            include_dashboard=False,
        )
    resources = ray.cluster_resources()
    n_cpus = int(resources.get("CPU", os.cpu_count() or 8))
    total = len(items)
    logger.info("%s via ray: cpus=%d batch=%d keys=%d", log_label, n_cpus, batch, total)

    state_ref = ray.put(state)
    remote = ray.remote(max_retries=max_task_retries, num_cpus=1)(_run_batch)
    batch_iter = _chunked(items, max(1, batch))

    def _submit() -> Any:
        nxt = next(batch_iter, None)
        if nxt is None:
            return None
        return remote.remote(worker_qualname, installer_qualname, state_ref, nxt)

    # Measured adaptive admission (AIMD). ``target`` is the in-flight task ceiling; it
    # grows by one whenever measured memory is comfortable and halves the moment it is
    # tight, so concurrency converges to what the node actually holds -- climbing to full
    # width on a light build, settling low while a heavy build's biggest records are in
    # flight (they arrive first, so the target drops early and recovers as they drain).
    # Never below 1 (a single outsized record still runs, alone). Ray's kill threshold is
    # only the backstop.
    target = max(1, n_cpus // 8)
    pending: list[Any] = []
    buffer: list[Any] = []
    done = 0
    exhausted = False
    mem = _mem_frac()
    try:
        while True:
            while not exhausted and len(pending) < min(target, n_cpus):
                fut = _submit()
                if fut is None:
                    exhausted = True
                    break
                pending.append(fut)
            if not pending:
                break
            ready, pending = ray.wait(pending, num_returns=1, timeout=None)
            for ref in ready:
                results = ray.get(ref)
                buffer.extend(results)
                done += len(results)
            mem = _mem_frac()
            if mem > _MEM_HIGH:
                target = max(1, target // 2)
            elif mem < _MEM_LOW and target < n_cpus:
                target += 1
            while len(buffer) >= chunk_size:
                yield buffer[:chunk_size]
                del buffer[:chunk_size]
            if done and done % (chunk_size * 5) < batch:
                logger.info(
                    "%s via ray: %d/%d target=%d mem=%.0f%%",
                    log_label, done, total, target, mem * 100,
                )
        if buffer:
            yield buffer
    finally:
        # Only tear down the cluster if we started it, so a caller that
        # pre-initialised Ray keeps ownership. Runs on error or generator-close
        # too; the graceful shutdown is what triggers Ray's native session cleanup
        # (its /dev/shm plasma file included), so a failed drain doesn't leak.
        if started_here:
            ray.shutdown()
