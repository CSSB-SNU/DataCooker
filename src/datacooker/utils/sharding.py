"""Generic helpers for node-aware sharding."""

from __future__ import annotations

import hashlib
import logging
import os
from collections.abc import Callable, Sequence
from typing import TypeVar

logger = logging.getLogger(__name__)

_NODE_COUNT_ENV_KEYS = ("WORLD_SIZE", "SLURM_NTASKS", "OMPI_COMM_WORLD_SIZE")
_NODE_RANK_ENV_KEYS = ("RANK", "SLURM_PROCID", "OMPI_COMM_WORLD_RANK")

ItemT = TypeVar("ItemT")


def read_first_valid_int_env(keys: tuple[str, ...]) -> int | None:
    """Read the first valid integer from a set of environment variables."""
    for key in keys:
        value = os.environ.get(key)
        if value is None:
            continue
        try:
            return int(value)
        except ValueError:
            logger.warning("Ignoring non-integer env value %s=%s", key, value)
    return None


def resolve_node_config(
    *,
    node_rank: int | None,
    node_count: int | None,
) -> tuple[int, int]:
    """Resolve node rank/count from explicit args first, then environment variables."""
    resolved_node_count = (
        node_count if node_count is not None else read_first_valid_int_env(_NODE_COUNT_ENV_KEYS)
    )
    resolved_node_rank = (
        node_rank if node_rank is not None else read_first_valid_int_env(_NODE_RANK_ENV_KEYS)
    )

    if resolved_node_count is None:
        resolved_node_count = 1
    if resolved_node_rank is None:
        resolved_node_rank = 0

    if resolved_node_count < 1:
        msg = f"Invalid node_count={resolved_node_count}. Expected >= 1."
        raise ValueError(msg)
    if resolved_node_rank < 0 or resolved_node_rank >= resolved_node_count:
        msg = (
            f"Invalid node_rank={resolved_node_rank} for node_count="
            f"{resolved_node_count}."
        )
        raise ValueError(msg)

    return resolved_node_rank, resolved_node_count


def shard_items(
    items: Sequence[ItemT],
    *,
    node_rank: int,
    node_count: int,
) -> list[ItemT]:
    """Return the sub-sequence assigned to the given node by list position.

    WARNING: positional striping is only a stable disjoint partition if every
    build sees the exact same ordered ``items``. Across separate runs (e.g.
    re-running individual shards) the scan order can differ, so the same logical
    record lands in different shards and the merged output ends up with
    duplicate keys. Prefer :func:`shard_items_by_key` for content that is
    sharded across multiple independent jobs.
    """
    return [item for index, item in enumerate(items) if index % node_count == node_rank]


def shard_items_by_key(
    items: Sequence[ItemT],
    *,
    key_func: Callable[[ItemT], object],
    node_rank: int,
    node_count: int,
) -> list[ItemT]:
    """Return the items whose key hashes to ``node_rank``.

    Shards by a stable hash of each item's key (md5 of its string form) rather
    than list position, so a given key always maps to the same shard regardless
    of scan order, file-set changes, or which run built which shard. This makes
    the shard set a true disjoint partition by key even when shards are built or
    re-run independently. ``hashlib`` is used (not the builtin ``hash``, which is
    salted per-process and would not be reproducible across runs).
    """
    selected: list[ItemT] = []
    for item in items:
        digest = hashlib.md5(str(key_func(item)).encode()).hexdigest()  # noqa: S324
        if int(digest, 16) % node_count == node_rank:
            selected.append(item)
    return selected
