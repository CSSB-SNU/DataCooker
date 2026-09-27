"""Utilities for resolving dotted Python object paths."""

from __future__ import annotations

from importlib import import_module
from typing import Any


def resolve_object(path: str) -> Any:
    """Resolve a dotted ``module.attr`` path into a Python object."""
    if "." not in path:
        msg = f"Expected a dotted path like 'package.module.object', got '{path}'."
        raise ValueError(msg)
    module_name, attr_name = path.rsplit(".", 1)
    module = import_module(module_name)
    return getattr(module, attr_name)
