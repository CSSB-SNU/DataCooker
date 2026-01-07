# API Reference

This section documents the public surface of `datacooker`. Module paths are relative to `src/datacooker`.

## `core.py`
- `ParsingCache(key_transform: Callable[[str], tuple[str, ...]] | None = None)`: nested, transformable key-value store for intermediate data. Methods:
  - `add_data(name, data)`: insert if not already present.
  - `__contains__(name)`: membership check.
  - `__getitem__(name)`: retrieve stored value.
  - `keys()`: flattened list of stored keys.
- `RecipeBook`: builder for ordered recipe steps.
  - `add(targets, instruction, inputs)`: register a step. `targets` can be a tuple of `(name, type)` pairs or a list of such tuples. `inputs` supports `args`, `kwargs`, and static `params`.
  - `__contains__(target_name)`: whether a target is defined.
  - `__getitem__(target_name)`: retrieve the `Recipe` for a given target.
  - `targets()`: list declared targets.
- `Cooker(parse_cache, recipebook, targets=None)`: orchestrator.
  - `prep(data_dict, fields=None)`: seed the cache.
  - `cook()`: execute all targets in dependency order.
  - `serve(targets=None)`: return selected targets (or defaults from the recipe book).
- `parse(recipe_path, file_path, load_func, transform_func=None, targets=None, **extra_kwargs)`: helper to load a file, run the recipe on it, and return results.
- `rebuild(recipe_path, datadict, transform_func=None, targets=None, **extra_kwargs)`: same as `parse` but starts from an existing `datadict`.

## `recipe.py`
- `Variable(name, type)`: typed target definition.
- `Inputs(args=(), kwargs={}, params={})`: inputs container used inside `Recipe`.
- `Recipe(targets, instruction, inputs)`: immutable step description.
- `RecipeError`: raised when a target is missing.

## Patterns
- **Optional targets**: declare a target type as `type | None` so the cooker returns `None` when missing instead of raising.
- **Wildcard args**: include glob characters (e.g., `"input_*"`) in `args` to gather all matching cache entries.
- **Multiple outputs**: when a step defines multiple targets, the instruction must return a tuple with the same arity.
