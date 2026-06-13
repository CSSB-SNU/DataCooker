"""Backward-compatible exports for the historical core module."""

from __future__ import annotations

from .api import execute, parse, parse_dict, parse_file, rebuild
from .executor import Cooker
from .loading import load_recipe
from .protocols import ConvertFunc, LoadFunc, TransformFunc

__all__ = [
    "ConvertFunc",
    "Cooker",
    "LoadFunc",
    "TransformFunc",
    "execute",
    "load_recipe",
    "parse",
    "parse_dict",
    "parse_file",
    "rebuild",
]
