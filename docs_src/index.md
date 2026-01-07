# DataCooker (Concepts)

DataCooker is a small recipe engine for biomolecular data. You describe *what* you want to build (targets) and *how* to compute them (instructions), and DataCooker resolves dependencies, executes steps, and keeps intermediate results in a cache.

## Core concepts
- **RecipeBook**: ordered steps that declare targets, an instruction callable, and inputs (args/kwargs/params) pointing to other targets.
- **ParsingCache**: keyed store for intermediate and final values; supports custom key transforms and flattening.
- **Cooker**: orchestrates dependency resolution and execution, expanding wildcards, detecting cycles, and writing outputs back to the cache.

## Execution flow
1. **Prep** the cache with initial fields (`cooker.prep(data_dict)`).
2. **Resolve** dependencies for requested targets (including wildcard inputs and optional targets).
3. **Cook** instructions in order.
4. **Serve** requested targets or defaults from the recipe book.

## What you can build
- Parse A3M/FASTA alignments, mmCIF/CCD files, and package graph data into LMDBs using the provided recipe books.
- Combine built-in instructions with your own to assemble custom pipelines tailored to new datasets.

## Declaring steps (example)
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

## Error handling
- Cycles raise `RuntimeError`; missing targets raise `RecipeError`.
- Type mismatches in multi-target returns raise `ValueError` with expected arity.
- Duplicate target declarations also raise `ValueError`.

## Next steps
- [Getting Started](getting-started.md): install from GitHub and run a first recipe.
- [Tutorials](tutorials.md): minimal custom recipe and packaged pipeline examples.
- [API](api.md): Python API reference for `datacooker`.
