"""Backward-compatible exports for the historical core module."""

from __future__ import annotations

from .api import describe, execute, parse, parse_dict, parse_file, rebuild, visualize
from .executor import Cooker
from .loading import load_recipe
from .protocols import (
    ConvertFunc,
    DeserializeFunc,
    LoadFunc,
    SerializeFunc,
    TransformFunc,
)
from .utils.db import (
    build_lmdb,
    count_lmdb_entries,
    extract_lmdb_keys,
    merge_lmdb_shards,
    read_all_lmdb_raw,
    read_lmdb,
    read_lmdb_raw,
    rebuild_lmdb,
)
from .utils.importing import resolve_object
from .utils.paths import scan_paths
from .utils.sharding import resolve_node_config, shard_items

__all__ = [
    "ConvertFunc",
    "Cooker",
    "DeserializeFunc",
    "LoadFunc",
    "SerializeFunc",
    "TransformFunc",
    "build_lmdb",
    "count_lmdb_entries",
    "describe",
    "execute",
    "extract_lmdb_keys",
    "load_recipe",
    "merge_lmdb_shards",
    "parse",
    "parse_dict",
    "parse_file",
    "read_all_lmdb_raw",
    "read_lmdb",
    "read_lmdb_raw",
    "rebuild",
    "rebuild_lmdb",
    "resolve_node_config",
    "resolve_object",
    "scan_paths",
    "shard_items",
    "visualize",
]
