"""Recipe loading helpers."""

from __future__ import annotations

import importlib.util
from collections.abc import Sequence
from pathlib import Path

from .recipe import RecipeBook


def normalize_targets(targets: Sequence[str] | str | None) -> list[str] | None:
    """Normalize target specifications into a list."""
    if targets is None:
        return None
    if isinstance(targets, str):
        return [targets]
    return list(targets)


def load_recipe(recipebook: RecipeBook | str | Path) -> tuple[RecipeBook, list[str] | None]:
    """Load a recipe from a file path or return the recipe directly."""
    if isinstance(recipebook, RecipeBook):
        return recipebook, recipebook.default_targets

    recipe_path = Path(recipebook).resolve()
    if not recipe_path.exists():
        msg = f"Recipe file '{recipe_path}' does not exist."
        raise FileNotFoundError(msg)

    module_name = recipe_path.stem
    spec = importlib.util.spec_from_file_location(module_name, recipe_path)
    if spec is None or spec.loader is None:
        msg = f"Could not load module from '{recipe_path}'."
        raise ImportError(msg)

    recipe_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(recipe_module)
    loaded_recipe = getattr(recipe_module, "RECIPE", None)
    if loaded_recipe is None or not isinstance(loaded_recipe, RecipeBook):
        msg = f"'RECIPE' not found or invalid in '{recipe_path}'."
        raise AttributeError(msg)
    default_targets = getattr(recipe_module, "TARGETS", loaded_recipe.default_targets)

    return loaded_recipe, normalize_targets(default_targets)
