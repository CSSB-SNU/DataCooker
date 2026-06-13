"""Core recipe engine for declarative data processing workflows."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

from .api import (
    describe,
    execute,
    parse,
    parse_dict,
    parse_file,
    rebuild,
    visualize,
)
from .cache import ExecutionContext, ParsingCache
from .core import (
    ConvertFunc,
    Cooker,
    DeserializeFunc,
    LoadFunc,
    SerializeFunc,
    TransformFunc,
    load_recipe,
)
from .errors import (
    CycleError,
    DataCookerError,
    DuplicateTargetError,
    InstructionOutputError,
    InvalidRecipeError,
    MissingDependencyError,
    MissingTargetError,
    StepExecutionError,
    UnknownTargetError,
)
from .recipe import Inputs, Recipe, RecipeBook, RecipeError, Variable, variable
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
    "CycleError",
    "DataCookerError",
    "DeserializeFunc",
    "DuplicateTargetError",
    "ExecutionContext",
    "Inputs",
    "InstructionOutputError",
    "InvalidRecipeError",
    "LoadFunc",
    "MissingDependencyError",
    "MissingTargetError",
    "ParsingCache",
    "Recipe",
    "RecipeBook",
    "RecipeError",
    "SerializeFunc",
    "StepExecutionError",
    "TransformFunc",
    "UnknownTargetError",
    "Variable",
    "__version__",
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
    "variable",
    "visualize",
]
