"""Internal helpers for optional joblib-backed batch execution."""

from __future__ import annotations

import logging
from collections.abc import Iterator, Mapping, Sequence
from functools import lru_cache
from importlib import import_module
from typing import Any

logger = logging.getLogger(__name__)


def require_joblib(*, feature_name: str) -> tuple[Any, Any]:
    """Return ``joblib.Parallel`` and ``joblib.delayed`` or raise a clear error."""
    try:
        joblib_module = _load_joblib_module()
    except ModuleNotFoundError as exc:
        msg = (
            f"{feature_name} requires the optional 'joblib' dependency. "
            "Install DataCooker with joblib support before using this feature."
        )
        raise ModuleNotFoundError(msg) from exc
    return joblib_module.Parallel, joblib_module.delayed


def run_serial_test(
    process_item: Any,
    items: Sequence[Any],
    *,
    enabled: bool,
) -> None:
    """Run a small serial sample before launching parallel work."""
    if not enabled or not items:
        return

    for item in items[:10]:
        result = process_item(item)
        if result is None:
            continue
        # A 1:N producer returns a list of records for one item; every other one
        # returns a single record. Either way the error is each record's last field.
        for record in result if isinstance(result, list) else [result]:
            if record[-1] is not None:
                raise record[-1]


def iter_parallel_chunks(
    *,
    items: Sequence[Any],
    chunk_size: int,
    n_jobs: int,
    process_item: Any,
    prefer: str = "processes",
    log_label: str,
    backend: str | None = None,
    ray_spec: Mapping[str, Any] | None = None,
) -> Iterator[list[Any]]:
    """Yield chunk results from a logged parallel loop (joblib, or Ray).

    Default (``backend`` unset): joblib with the ``prefer`` hint (loky processes /
    threads). Used by ``extract`` and ``processing.batch``.

    ``backend="ray"``: the Ray fan-out (``_ray.iter_ray_chunks``) used by
    ``build_lmdb`` / ``rebuild_lmdb`` -- concurrency self-regulates to measured node
    memory (admission control) instead of a predicted ``n_jobs``, so it needs no
    per-DB tuning. Requires ``ray_spec`` (``{worker_qualname, installer_qualname,
    state, ...}``) because a Ray worker has no fork to inherit module globals from;
    the state is shipped once via the plasma object store. ``process_item`` /
    ``n_jobs`` / ``prefer`` are unused on this path.
    """
    if backend == "ray":
        if ray_spec is None:
            msg = "backend='ray' requires ray_spec (worker/installer/state)."
            raise ValueError(msg)
        from datacooker._ray import iter_ray_chunks

        ray_kwargs = {k: v for k, v in ray_spec.items() if k != "state"}
        yield from iter_ray_chunks(
            items=items,
            chunk_size=chunk_size,
            log_label=log_label,
            state=ray_spec["state"],
            **ray_kwargs,
        )
        return
    parallel, delayed = require_joblib(feature_name=log_label)
    for start in range(0, len(items), chunk_size):
        logger.info(
            "%s %d to %d / %d",
            log_label,
            start,
            min(start + chunk_size, len(items)),
            len(items),
        )
        chunk = items[start : start + chunk_size]
        yield parallel(n_jobs=n_jobs, verbose=10, prefer=prefer)(
            delayed(process_item)(item) for item in chunk
        )


@lru_cache(maxsize=1)
def _load_joblib_module() -> Any:
    return import_module("joblib")
