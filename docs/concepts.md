# Concepts & Architecture

## Core pieces
- `RecipeBook`: ordered collection of steps. Each step declares `targets`, an `instruction` callable, and `inputs` that point to other targets (args/kwargs) plus static `params`.
- `ParsingCache`: a keyed store for intermediate and final values. It can flatten keys (`foo.bar`) or apply a custom `key_transform` for nested access.
- `Cooker`: orchestrates execution. It resolves dependencies, handles wildcard inputs (glob patterns like `input*`), detects cycles, runs instructions, and stores outputs back into the cache.
- Helpers: `parse()` loads an input file, seeds the cache, runs a recipe book from disk, and returns selected targets. `rebuild()` is similar but starts from an in-memory `datadict`.

## Execution model
1. **Prep**: seed the cache with initial fields via `cooker.prep(data_dict)`.
2. **Resolve**: for each requested target, resolve all upstream inputs:
   - Args/kwargs in `Inputs` are mapped by name to other targets.
   - Wildcard args are expanded against existing cache keys.
   - Optional targets (`type | None`) return `None` if not declared.
3. **Cook**: run instructions in dependency order, caching outputs as they are produced.
4. **Serve**: retrieve specific targets or all defaults defined by the recipe book.

## Declaring steps
```python
from datacooker import RecipeBook

book = RecipeBook()
book.add(
    targets=(("features", dict),),
    instruction=lambda seq: {"len": len(seq)},
    inputs={"args": (("sequence", str),)},
)
book.add(
    targets=(("summary", str),),
    instruction=lambda feats: f"length={feats['len']}",
    inputs={"args": (("features", dict),)},
)
```

## Working with files vs memory
- Use `datacooker.core.parse()` when your starting point is a file path and you have a loader function that returns a data dictionary.
- Use `datacooker.core.rebuild()` when you already have a populated `datadict` (e.g., from a previous pipeline stage or cached LMDB read).

## Error handling
- Cycles raise `RuntimeError`.
- Missing targets raise `RecipeError`.
- Type mismatches in returned tuples raise `ValueError` with expected target counts.
- Duplicate targets in a recipe book raise `ValueError`.
