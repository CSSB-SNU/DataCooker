# Write a recipe

A recipe is a `RecipeBook` of steps. Each step declares **targets** (named
outputs with a type), an **instruction** (a callable), and **inputs** that bind
the instruction's parameters to other targets in the cache.

## A two-step recipe

```python
from datacooker import RecipeBook, execute

book = RecipeBook()

# Step 1: sequence -> features
book.add(
    targets=(("features", dict),),
    instruction=lambda sequence: {"length": len(sequence)},
    inputs={"args": (("sequence", str),)},
)

# Step 2: features -> summary
book.add(
    targets=(("summary", str),),
    instruction=lambda features: f"length={features['length']}",
    inputs={"args": (("features", dict),)},
)

print(execute(book, {"sequence": "MKTAYIAK"}, targets=["summary"]))
# {'summary': 'length=8'}
```

## Binding inputs

The `inputs` mapping supports three keys:

- `args` — positional inputs, a tuple of `(target_name, type)` pairs.
- `kwargs` — keyword inputs, a mapping of `param_name -> (target_name, type)`.
- `params` — literal constants passed straight to the instruction.

```python
def trim(sequence: str, *, max_len: int) -> str:
    return sequence[:max_len]

book.add(
    targets=(("trimmed", str),),
    instruction=trim,
    inputs={
        "args": (("sequence", str),),
        "params": {"max_len": 4},
    },
)
```

## Instruction factories

For configurable instructions, return the worker from a factory so the recipe
stays declarative:

```python
def build_summary(prefix: str):
    def _worker(features: dict) -> str:
        return f"{prefix}={features['length']}"
    return _worker

book.add(
    targets=(("summary", str),),
    instruction=build_summary("length"),
    inputs={"args": (("features", dict),)},
)
```

## Multiple targets and inspection

A single step may return several targets. A `RecipeBook` can also describe and
validate itself before running:

```python
book.describe(["summary"])     # reachable graph for a target
book.execution_order(["summary"])
book.validate()                # raises on cycles / duplicates / missing inputs
```

To run against a prepared dict use `parse_dict`; to run only the default
targets, call `execute(book, inputs)` without `targets`.
