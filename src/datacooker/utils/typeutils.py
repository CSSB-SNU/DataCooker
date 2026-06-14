"""Shared helpers for interpreting Python type annotations.

These were previously duplicated between :mod:`datacooker.recipe` and
:mod:`datacooker.executor` (the latter inlined ``type(None) in get_args(...)``).
Centralizing them keeps annotation semantics consistent across the engine.
"""

from __future__ import annotations

from typing import Any, get_args, get_origin

NoneType = type(None)


def allows_none(annotation: object) -> bool:
    """Return whether ``annotation`` admits ``None`` (e.g. ``str | None``)."""
    if annotation is None or annotation is NoneType:
        return True
    return NoneType in get_args(annotation)


def is_supported_annotation(annotation: object) -> bool:
    """Return whether ``annotation`` is a usable variable type declaration.

    Accepted forms:
      * ``Any``
      * a concrete class (``str``, ``MyModel``)
      * a typing construct with an origin (``list[int]``, ``str | None``,
        ``Optional[X]``, ``dict[str, int]``)

    Everything else (bare ``None``, arbitrary instances, malformed objects) is
    rejected so that invalid declarations fail fast instead of silently being
    treated as a type.
    """
    if annotation is Any:
        return True
    if isinstance(annotation, type):
        return True
    return get_origin(annotation) is not None
