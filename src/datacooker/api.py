"""Public execution APIs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .cache import ExecutionContext
from .executor import Cooker
from .loading import load_recipe
from .protocols import LoadFunc, TransformFunc
from .recipe import RecipeBook


def execute(
    recipebook: RecipeBook | str | Path,
    inputs: Mapping[str, Any],
    *,
    transform_func: TransformFunc | None = None,
    targets: Sequence[str] | str | None = None,
    step_names: Sequence[str] | str | None = None,
    tags: Sequence[str] | str | None = None,
    namespaces: Sequence[str] | str | None = None,
    validate: bool = True,
) -> dict[str, Any]:
    """Execute a static workflow graph and return the requested outputs."""
    selected_recipe = _select_recipebook(
        recipebook,
        targets=targets,
        step_names=step_names,
        tags=tags,
        namespaces=namespaces,
    )
    parse_cache = ExecutionContext(transform_func)
    cooker = Cooker(parse_cache=parse_cache, recipebook=selected_recipe)
    cooker.prep(dict(inputs))
    cooker.cook(targets=targets, validate=validate)
    return cooker.serve(targets=targets)


def describe(
    recipebook: RecipeBook | str | Path,
    *,
    targets: Sequence[str] | str | None = None,
    step_names: Sequence[str] | str | None = None,
    tags: Sequence[str] | str | None = None,
    namespaces: Sequence[str] | str | None = None,
    available_inputs: Sequence[str] | set[str] | None = None,
) -> str:
    """Describe the reachable workflow graph for the requested targets."""
    loaded_recipe = _select_recipebook(
        recipebook,
        targets=targets,
        step_names=step_names,
        tags=tags,
        namespaces=namespaces,
    )
    resolved_targets = loaded_recipe.default_targets if targets is None else targets
    return loaded_recipe.describe(
        targets=resolved_targets,
        available_inputs=available_inputs,
    )


def visualize(
    recipebook: RecipeBook | str | Path,
    *,
    output_format: str = "mermaid",
    targets: Sequence[str] | str | None = None,
    step_names: Sequence[str] | str | None = None,
    tags: Sequence[str] | str | None = None,
    namespaces: Sequence[str] | str | None = None,
    available_inputs: Sequence[str] | set[str] | None = None,
) -> str:
    """Render a workflow graph in Mermaid or DOT format."""
    loaded_recipe = _select_recipebook(
        recipebook,
        targets=targets,
        step_names=step_names,
        tags=tags,
        namespaces=namespaces,
    )
    resolved_targets = loaded_recipe.default_targets if targets is None else targets
    return loaded_recipe.visualize(
        output_format=output_format,
        targets=resolved_targets,
        available_inputs=available_inputs,
    )


def parse_file(
    recipe_path: RecipeBook | str | Path,
    file_path: Path,
    load_func: LoadFunc | None,
    inputs: Mapping[str, Any] | None = None,
    transform_func: TransformFunc | None = None,
    targets: Sequence[str] | str | None = None,
    step_names: Sequence[str] | str | None = None,
    tags: Sequence[str] | str | None = None,
    namespaces: Sequence[str] | str | None = None,
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
        step_names=step_names,
        tags=tags,
        namespaces=namespaces,
        validate=validate,
    )


def parse_dict(
    recipe_path: RecipeBook | str | Path,
    datadict: Mapping[str, Any],
    transform_func: TransformFunc | None = None,
    targets: Sequence[str] | str | None = None,
    step_names: Sequence[str] | str | None = None,
    tags: Sequence[str] | str | None = None,
    namespaces: Sequence[str] | str | None = None,
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
        step_names=step_names,
        tags=tags,
        namespaces=namespaces,
        validate=validate,
    )


def parse(
    recipe_path: RecipeBook | str | Path,
    file_path: Path,
    load_func: LoadFunc,
    transform_func: TransformFunc | None = None,
    targets: Sequence[str] | str | None = None,
    step_names: Sequence[str] | str | None = None,
    tags: Sequence[str] | str | None = None,
    namespaces: Sequence[str] | str | None = None,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Backward-compatible alias for :func:`parse_file`."""
    return parse_file(
        recipe_path,
        file_path,
        load_func,
        transform_func=transform_func,
        targets=targets,
        step_names=step_names,
        tags=tags,
        namespaces=namespaces,
        **extra_kwargs,
    )


def rebuild(
    recipe_path: RecipeBook | str | Path,
    datadict: Mapping[str, Any],
    transform_func: TransformFunc | None = None,
    targets: Sequence[str] | str | None = None,
    step_names: Sequence[str] | str | None = None,
    tags: Sequence[str] | str | None = None,
    namespaces: Sequence[str] | str | None = None,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Backward-compatible alias for :func:`parse_dict`."""
    return parse_dict(
        recipe_path,
        datadict,
        transform_func=transform_func,
        targets=targets,
        step_names=step_names,
        tags=tags,
        namespaces=namespaces,
        **extra_kwargs,
    )


def _select_recipebook(
    recipebook: RecipeBook | str | Path,
    *,
    targets: Sequence[str] | str | None = None,
    step_names: Sequence[str] | str | None = None,
    tags: Sequence[str] | str | None = None,
    namespaces: Sequence[str] | str | None = None,
) -> RecipeBook:
    has_filters = any(selection is not None for selection in (step_names, tags, namespaces))
    loaded_recipe, default_targets = load_recipe(recipebook)
    if targets is None and not has_filters:
        if default_targets is None:
            return loaded_recipe
        return loaded_recipe.subset(targets=default_targets)

    return loaded_recipe.subset(
        targets=targets,
        step_names=step_names,
        tags=tags,
        namespaces=namespaces,
    )
