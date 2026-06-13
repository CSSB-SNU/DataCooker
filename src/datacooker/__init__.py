"""Core recipe engine for declarative data processing workflows."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

from .cache import ParsingCache
from .core import ConvertFunc, Cooker, LoadFunc, TransformFunc, parse, rebuild
from .recipe import RecipeBook

__all__ = [
    "ConvertFunc",
    "Cooker",
    "LoadFunc",
    "ParsingCache",
    "RecipeBook",
    "TransformFunc",
    "__version__",
    "parse",
    "rebuild",
]
