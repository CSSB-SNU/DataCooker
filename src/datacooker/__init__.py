"""ML-ready biomolecular data handling package.

This package provides tools for representing, manipulating, and analyzing biomolecular
structures and sequences. It is designed to integrate seamlessly with machine learning
workflows by offering standardized data representations (e.g., arrays) and utilities for
scientific computing.
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

from .cache import ParsingCache
from .core import Cooker, LoadFunc, TransformFunc, parse
from .recipe import RecipeBook

__all__ = [
    "Cooker",
    "LoadFunc",
    "ParsingCache",
    "RecipeBook",
    "TransformFunc",
    "__version__",
    "parse",
]
