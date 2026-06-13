"""Public core API for DataCooker's static workflow engine."""

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
from .executor import Cooker
from .loading import load_recipe
from .protocols import (
    ConvertFunc,
    DeserializeFunc,
    LoadFunc,
    SerializeFunc,
    TransformFunc,
)
from .recipe import Inputs, Recipe, RecipeBook, RecipeError, Variable, variable
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
    "describe",
    "execute",
    "load_recipe",
    "parse",
    "parse_dict",
    "parse_file",
    "rebuild",
    "resolve_node_config",
    "resolve_object",
    "scan_paths",
    "shard_items",
    "variable",
    "visualize",
]
