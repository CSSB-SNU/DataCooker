"""Core recipe engine for declarative data processing workflows."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

from .api import (
    execute,
    parse,
    parse_dict,
    parse_file,
    rebuild,
)
from .cache import ExecutionContext, ParsingCache
from .core import (
    ConvertFunc,
    Cooker,
    LoadFunc,
    TransformFunc,
    load_recipe,
)
from .errors import (
    CycleError,
    DataCookerError,
    DuplicateTargetError,
    InvalidRecipeError,
    MissingDependencyError,
    MissingTargetError,
    UnknownTargetError,
)
from .recipe import RecipeBook

__all__ = [
    "ConvertFunc",
    "Cooker",
    "CycleError",
    "DataCookerError",
    "DuplicateTargetError",
    "ExecutionContext",
    "InvalidRecipeError",
    "LoadFunc",
    "MissingDependencyError",
    "MissingTargetError",
    "ParsingCache",
    "RecipeBook",
    "TransformFunc",
    "UnknownTargetError",
    "__version__",
    "execute",
    "load_recipe",
    "parse",
    "parse_dict",
    "parse_file",
    "rebuild",
]
