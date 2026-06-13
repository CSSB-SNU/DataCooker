"""Public execution APIs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .cache import ExecutionContext
from .executor import Cooker
from .protocols import LoadFunc, TransformFunc
from .recipe import RecipeBook


def execute(
    recipebook: RecipeBook | str | Path,
    inputs: Mapping[str, Any],
    *,
    transform_func: TransformFunc | None = None,
    targets: Sequence[str] | str | None = None,
    validate: bool = True,
) -> dict[str, Any]:
    """Execute a static workflow graph and return the requested outputs."""
    parse_cache = ExecutionContext(transform_func)
    cooker = Cooker(parse_cache=parse_cache, recipebook=recipebook)
    cooker.prep(dict(inputs))
    cooker.cook(targets=targets, validate=validate)
    return cooker.serve(targets=targets)


def parse_file(
    recipe_path: RecipeBook | str | Path,
    file_path: Path,
    load_func: LoadFunc | None,
    inputs: Mapping[str, Any] | None = None,
    transform_func: TransformFunc | None = None,
    targets: Sequence[str] | str | None = None,
    validate: bool = True,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Load a file, merge additional inputs, and execute the requested recipe."""
    data_dict = dict(inputs or {})
    if load_func is not None:
        data_dict.update(load_func(file_path))
    else:
        data_dict["file_path"] = file_path
    data_dict.update(extra_kwargs)
    return execute(
        recipe_path,
        data_dict,
        transform_func=transform_func,
        targets=targets,
        validate=validate,
    )


def parse_dict(
    recipe_path: RecipeBook | str | Path,
    datadict: Mapping[str, Any],
    transform_func: TransformFunc | None = None,
    targets: Sequence[str] | str | None = None,
    validate: bool = True,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Execute a recipe against an already prepared input dictionary."""
    data_dict = dict(datadict)
    data_dict.update(extra_kwargs)
    return execute(
        recipe_path,
        data_dict,
        transform_func=transform_func,
        targets=targets,
        validate=validate,
    )


def parse(
    recipe_path: RecipeBook | str | Path,
    file_path: Path,
    load_func: LoadFunc,
    transform_func: TransformFunc | None = None,
    targets: Sequence[str] | str | None = None,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Backward-compatible alias for :func:`parse_file`."""
    return parse_file(
        recipe_path,
        file_path,
        load_func,
        transform_func=transform_func,
        targets=targets,
        **extra_kwargs,
    )


def rebuild(
    recipe_path: RecipeBook | str | Path,
    datadict: Mapping[str, Any],
    transform_func: TransformFunc | None = None,
    targets: Sequence[str] | str | None = None,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Backward-compatible alias for :func:`parse_dict`."""
    return parse_dict(
        recipe_path,
        datadict,
        transform_func=transform_func,
        targets=targets,
        **extra_kwargs,
    )
