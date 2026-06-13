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
    LoadFunc,
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

__all__ = [
    "ConvertFunc",
    "Cooker",
    "CycleError",
    "DataCookerError",
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
    "variable",
    "visualize",
]
