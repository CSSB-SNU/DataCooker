# DataCooker

`DataCooker` is an API library for declaring, validating, and executing static
data workflows over a fixed dependency graph.

## Scope

This repository intentionally keeps only reusable workflow primitives:

- `RecipeBook` for declaring workflow steps
- `Cooker` for resolving requested targets over a static graph
- `ExecutionContext` / `ParsingCache` for runtime values
- `execute()` for running workflows from in-memory inputs
- `parse_dict()` and `parse_file()` as convenience wrappers

Domain-specific pipelines, IO loaders, model calls, and data products are
expected to live in downstream repositories.

## What It Is

DataCooker models a workflow as:

- named target artifacts
- deterministic dependency edges between artifacts
- step functions that consume declared inputs and produce declared outputs

The graph shape is static. The runtime values are dynamic.

## What It Is Not

DataCooker does not try to be:

- a scheduler
- a distributed task runner
- a dynamic workflow engine that mutates the graph at runtime
- a general orchestration platform like Airflow

It is a small execution core for static data workflows.

## Installation

```bash
pip install -e .
```

## Quick Start

```python
from datacooker import RecipeBook, execute


def add(a: int, b: int) -> int:
    return a + b


recipe = RecipeBook().add(
    (("sum", int),),
    add,
    {"args": (("a", int), ("b", int))},
)

result = execute(recipe, {"a": 1, "b": 2}, targets="sum")
assert result["sum"] == 3
```

## Multi-Step Example

```python
from datacooker import RecipeBook, execute


recipe = RecipeBook()
recipe.add(
    (("sum", int),),
    lambda a, b: a + b,
    {"args": (("a", int), ("b", int))},
)
recipe.add(
    (("label", str),),
    lambda total: f"total={total}",
    {"args": (("sum", int),)},
)

result = execute(recipe, {"a": 2, "b": 5}, targets="label")
assert result == {"label": "total=7"}
```

## File-Based Recipes

Recipe modules can expose:

- `RECIPE`: a `RecipeBook`
- `TARGETS`: optional default targets

Then run them with:

```python
from pathlib import Path

from datacooker import parse_dict

result = parse_dict(Path("my_recipe.py"), {"a": 1, "b": 2})
```

## Public API

- `RecipeBook`: declare workflow steps
- `Cooker`: stateful executor for advanced usage
- `execute`: primary execution entrypoint
- `parse_dict`, `parse_file`: convenience wrappers
- `load_recipe`: load `RECIPE` / `TARGETS` from a module path
- `ExecutionContext`: runtime key-value store

## Current Guarantees

- returned outputs are always `dict[str, Any]`
- target selection supports full-graph or subset execution
- recipe validation checks missing dependencies and cycles
- wildcard args resolve against keys already present in the execution context
