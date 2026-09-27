"""Internal helpers for optional joblib-backed batch execution."""

from __future__ import annotations

import logging
from collections.abc import Iterator, Sequence
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
        if result is not None and result[-1] is not None:
            raise result[-1]


def iter_parallel_chunks(
    *,
    items: Sequence[Any],
    chunk_size: int,
    n_jobs: int,
    process_item: Any,
    prefer: str,
    log_label: str,
) -> Iterator[list[Any]]:
    """Yield chunk results from a logged joblib-backed parallel loop."""
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
