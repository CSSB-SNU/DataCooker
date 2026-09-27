"""Backward-compatible exports for the historical ``datacooker.core`` path."""

from __future__ import annotations

from .api import describe, execute, parse, parse_dict, parse_file, rebuild, visualize
from .executor import Cooker
from .lmdb import (
    build_lmdb,
    count_lmdb_entries,
    default_lmdb_key,
    extract_lmdb_keys,
    extract_lmdb_records,
    filter_pending_lmdb_paths,
    merge_lmdb_shards,
    read_all_lmdb_raw,
    read_lmdb,
    read_lmdb_raw,
    rebuild_lmdb,
)
from .loading import load_recipe
from .processing import parallel_process
from .protocols import (
    ConvertFunc,
    DeserializeFunc,
    LoadFunc,
    SerializeFunc,
    TransformFunc,
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
    "default_lmdb_key",
    "describe",
    "execute",
    "extract_lmdb_keys",
    "extract_lmdb_records",
    "filter_pending_lmdb_paths",
    "load_recipe",
    "merge_lmdb_shards",
    "parallel_process",
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
