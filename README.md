# DataCooker

`DataCooker` is a typed API library for declaring, validating, and executing
static data workflows over a fixed dependency graph.

## Scope

This repository intentionally keeps only reusable workflow primitives:

- `RecipeBook` for declaring workflow steps
- `Cooker` for resolving requested targets over a static graph
- `ExecutionContext` / `ParsingCache` for runtime values
- `execute()` for running workflows from in-memory inputs
- `parse_dict()` and `parse_file()` as convenience wrappers
- `describe()` for workflow introspection and debugging
- `visualize()` for Mermaid and Graphviz graph rendering
- `datacooker.utils.db` for recipe-driven LMDB build/rebuild/read flows

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
from datacooker import RecipeBook, execute, variable


def add(a: int, b: int) -> int:
    return a + b


recipe = RecipeBook().step(
    outputs=variable("sum", int),
    instruction=add,
    args=[variable("a", int), variable("b", int)],
)

result = execute(recipe, {"a": 1, "b": 2}, targets="sum")
assert result["sum"] == 3
```

## Multi-Step Example

```python
from datacooker import RecipeBook, execute, variable


recipe = RecipeBook()
recipe.step(
    outputs=variable("sum", int),
    instruction=lambda a, b: a + b,
    args=[variable("a", int), variable("b", int)],
)
recipe.step(
    outputs=variable("label", str),
    instruction=lambda total: f"total={total}",
    args=[variable("sum", int)],
)
recipe.set_default_targets("label")

result = execute(recipe, {"a": 2, "b": 5})
assert result == {"label": "total=7"}
```

## Workflow Introspection

```python
from datacooker import describe, visualize

print(
    describe(
        recipe,
        available_inputs={"a", "b"},
    )
)
```

Example output:

```text
Targets: label
Required inputs: a, b
Missing inputs: <none>
Execution order:
1. sum <- a, b
2. label <- sum
```

Mermaid output is available directly:

```python
print(visualize(recipe, available_inputs={"a", "b"}))
```

Graphviz DOT is also supported:

```python
print(visualize(recipe, output_format="dot", available_inputs={"a", "b"}))
```

## Database Utilities

`DataCooker` also ships reusable LMDB workflow helpers under
`datacooker.utils.db`.

They cover:

- building an LMDB from raw files plus a recipe
- rebuilding an LMDB through another recipe
- reading decoded entries
- merging shard databases
- listing keys

Serialization stays downstream-specific through callbacks, so domain projects
can keep their own byte format while reusing the workflow orchestration layer.

## File-Based Recipes

Recipe modules can expose:

- `RECIPE`: a `RecipeBook`
- `TARGETS`: optional default targets

If `TARGETS` is omitted, `RecipeBook.set_default_targets(...)` can define the
same behavior inside the module.

Then run them with:

```python
from pathlib import Path

from datacooker import parse_dict

result = parse_dict(Path("my_recipe.py"), {"a": 1, "b": 2})
```

## Public API

- `RecipeBook`: declare workflow steps, defaults, and graph metadata
- `Cooker`: stateful executor for advanced usage
- `execute`: primary execution entrypoint
- `describe`: graph introspection entrypoint
- `visualize`: graph rendering entrypoint
- `parse_dict`, `parse_file`: convenience wrappers
- `load_recipe`: load `RECIPE` / `TARGETS` from a module path
- `ExecutionContext`: runtime key-value store
- `datacooker.utils.db`: generic LMDB workflow helpers
- `datacooker.utils.sharding`: node-rank aware sharding helpers
- `Variable`, `Inputs`, `Recipe`, `variable`: typed declaration helpers
- `MissingDependencyError`, `StepExecutionError`, `InstructionOutputError`: structured runtime diagnostics

## Current Guarantees

- returned outputs are always `dict[str, Any]`
- target selection supports full-graph or subset execution
- recipe validation checks missing dependencies and cycles
- wildcard args resolve against keys already present in the execution context
- recipe books can expose default targets directly
- Mermaid and DOT graph rendering use the same static dependency model as execution
- runtime errors carry structured context such as target names, dependency chains, and step descriptions
- the package ships inline type information via `py.typed`

## Compatibility

The legacy `RecipeBook.add(...)`, `parse(...)`, and `rebuild(...)` APIs remain
available for existing downstream code. New code should prefer
`RecipeBook.step(...)` plus `execute(...)`.
