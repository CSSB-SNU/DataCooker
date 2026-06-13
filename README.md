# DataCooker

`DataCooker` is the core recipe engine extracted from the structure-oriented data
workflows that now live in downstream projects such as `StructCooker`.

## Scope

This repository intentionally keeps only the reusable core abstractions:

- `ParsingCache` for storing intermediate values
- `RecipeBook` for declaring recipe steps
- `Cooker` for resolving and executing dependency graphs
- `parse` and `rebuild` helpers for loading recipe modules

Domain-specific pipelines, configs, submit scripts, and structure-processing
recipes are expected to live in downstream repositories.

## Installation

```bash
pip install -e .
```

## Example

```python
from datacooker import Cooker, ParsingCache, RecipeBook


def add(a: int, b: int) -> int:
    return a + b


recipe = RecipeBook().add(
    (("sum", int),),
    add,
    {"args": (("a", int), ("b", int))},
)

cache = ParsingCache()
cooker = Cooker(cache, recipe, targets=["sum"])
cooker.prep({"a": 1, "b": 2})
cooker.cook()

result, _ = cooker.serve()
assert result["sum"] == 3
```
