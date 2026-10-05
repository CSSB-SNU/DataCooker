"""Node-local Ray execution with bounded submission and incremental result draining.

A run installs its recipe state once per worker. Python dictionaries are copied by
Ray deserialization; putting one in the object store does not make it zero-copy.
Domain workflows should therefore supply compact metadata or read-only file-backed
lookups. The scheduler respects the requested CPU limit, starts with one task, and
increases concurrency only after completed work. Memory pressure pauses admission
and reduces concurrency; it is a safety margin, not a proof that any item fits RAM.
"""
from __future__ import annotations

import logging
import os
import shutil
import tempfile
import time
import uuid
from collections.abc import Iterator, Sequence
from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from typing import Any

import psutil

logger = logging.getLogger(__name__)
_RESULT_BYTES = 64 * 1024**2


@dataclass
class Admission:
    """Grow after completed waves; recover gradually after memory pressure."""

    maximum: int
    target: int = 1
    completed: int = 0

    def observe(self, memory_fraction: float, completed: int) -> None:
        """Update the active-task budget from measured memory and completions."""
        if memory_fraction >= 0.75:
            self.target = max(1, self.target // 2)
            self.completed = 0
        elif memory_fraction < 0.70:
            self.completed += completed
            if self.completed >= self.target:
                self.completed = 0
                self.target = min(self.maximum, self.target + 1)

    def admits(self, pending: int, memory_fraction: float) -> bool:
        """Keep progress possible at width one, otherwise leave memory headroom."""
        return pending == 0 or (pending < self.target and memory_fraction < 0.70)


def _resolve(qualname: str) -> Any:
    module_name, attr = qualname.rsplit(".", 1)
    return getattr(import_module(module_name), attr)


def require_ray() -> Any:
    """Import the required execution backend."""
    return import_module("ray")


def _run_batch(
    worker_qualname: str,
    installer_qualname: str,
    state_refs: list[Any],
    run_id: str,
    keys: Sequence[Any],
) -> list[Any]:
    """Install state once per run, even when a caller reuses its Ray cluster.

    Nest the ObjectRef so Ray does not deserialize the entire state for every
    task before entering this function. A plain ObjectRef argument does exactly
    that, even if the installer subsequently skips it.
    """
    module = import_module(worker_qualname.rsplit(".", 1)[0])
    if getattr(module, "_RAY_STATE_RUN", None) != run_id:
        _resolve(installer_qualname)(require_ray().get(state_refs[0]))
        module.__dict__["_RAY_STATE_RUN"] = run_id
    worker = _resolve(worker_qualname)
    return [worker(key) for key in keys]


def _chunked(seq: Sequence[Any], size: int) -> Iterator[Sequence[Any]]:
    for offset in range(0, len(seq), size):
        yield seq[offset : offset + size]


def _result_bytes(result: Any) -> int:
    if isinstance(result, list):
        return sum(_result_bytes(item) for item in result)
    if isinstance(result, tuple) and len(result) >= 2 and isinstance(result[1], bytes):
        return len(result[1])
    return 0


def _cpu_limit(n_jobs: int) -> int:
    available = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    allocated = int(os.environ.get("SLURM_CPUS_PER_TASK", available))
    available = min(available, allocated)
    if n_jobs == 0:
        msg = "n_jobs must not be zero."
        raise ValueError(msg)
    return max(1, min(available, n_jobs if n_jobs > 0 else available + 1 + n_jobs))


def object_store_budget(scratch: str) -> int:
    """Bound each job's object store instead of reserving a fraction of the whole node."""
    allocated_mb = int(os.environ.get("SLURM_MEM_PER_NODE", "0"))
    if not allocated_mb:
        allocated_mb = int(os.environ.get("SLURM_MEM_PER_CPU", "0")) * _cpu_limit(-1)
    budget = allocated_mb * 1024**2 if allocated_mb else psutil.virtual_memory().available
    requested = int(os.environ.get("DATACOOKER_OBJECT_STORE_BYTES", min(8 * 1024**3, budget // 20)))
    minimum = 80 * 1024**2
    if requested < minimum or requested > budget // 4:
        msg = "Object store must be at least 80 MiB and at most a quarter of the job memory budget"
        raise ValueError(msg)
    shared = Path("/dev/shm")  # noqa: S108 - inspect system shared-memory capacity
    shared_free = shutil.disk_usage(shared).free if shared.is_dir() else 0
    if shared_free >= 2 * minimum:
        return min(requested, shared_free // 2)
    if shutil.disk_usage(scratch).free < 2 * requested:
        msg = "Insufficient shared memory and scratch space for bounded Ray object store"
        raise RuntimeError(msg)
    return requested


def iter_ray_chunks(
    *,
    items: Sequence[Any],
    chunk_size: int,
    worker_qualname: str,
    installer_qualname: str,
    state: Any,
    log_label: str,
    n_jobs: int = -1,
    batch: int = 1,
    memory_threshold: float = 0.9,
    max_task_retries: int = 3,
) -> Iterator[list[Any]]:
    """Yield completed results bounded by count, serialized bytes and flush interval."""
    if chunk_size < 1 or batch < 1:
        msg = "chunk_size and batch must be positive."
        raise ValueError(msg)
    if not items:
        return
    ray = require_ray()
    maximum = _cpu_limit(n_jobs)
    started_here = not ray.is_initialized()
    pending: list[Any] = []
    scratch = None
    try:
        if started_here:
            scratch = tempfile.TemporaryDirectory(prefix="dc-ray-")
            ray.init(
                address="local", _temp_dir=scratch.name, num_cpus=maximum,
                include_dashboard=False, logging_level=logging.WARNING,
                object_store_memory=object_store_budget(scratch.name),
                _system_config={
                    "memory_usage_threshold": memory_threshold,
                    "task_oom_retries": max_task_retries,
                },
            )
        maximum = min(maximum, max(1, int(ray.cluster_resources().get("CPU", 1))))
        admission = Admission(maximum)
        state_ref = ray.put(state)
        run_id = uuid.uuid4().hex
        remote = ray.remote(max_retries=max_task_retries, num_cpus=1)(_run_batch)
        batches = iter(_chunked(items, batch))
        exhausted = False
        buffer: list[Any] = []
        buffered_bytes = 0
        done = 0
        flushed = logged = time.monotonic()
        logger.info("%s: ray cpus=%d items=%d batch=%d", log_label, maximum, len(items), batch)
        while pending or not exhausted:
            memory = psutil.virtual_memory().percent / 100.0
            while not exhausted and admission.admits(len(pending), memory):
                keys = next(batches, None)
                if keys is None:
                    exhausted = True
                    break
                pending.append(remote.remote(worker_qualname, installer_qualname, [state_ref], run_id, keys))
            completed = 0
            if pending:
                ready, pending = ray.wait(pending, num_returns=1, timeout=1.0)
                for ref in ready:
                    results = ray.get(ref)
                    completed += len(results)
                    buffered_bytes += _result_bytes(results)
                    buffer.extend(results)
                    del results
                del ready
            done += completed
            memory = psutil.virtual_memory().percent / 100.0
            admission.observe(memory, int(completed > 0))
            now = time.monotonic()
            if buffer and (len(buffer) >= chunk_size or buffered_bytes >= _RESULT_BYTES
                           or now - flushed >= 0.5 or not pending):
                yield buffer
                buffer = []
                buffered_bytes = 0
                flushed = time.monotonic()
            if now - logged >= 10:
                logger.info("%s: completed=%d/%d active=%d limit=%d memory=%.1f%%",
                            log_label, done, len(items), len(pending), admission.target, memory * 100)
                logged = now
    finally:
        for ref in pending:
            ray.cancel(ref)
        if started_here:
            ray.shutdown()
        if scratch is not None:
            scratch.cleanup()
