# DataCooker

DataCooker is a small recipe engine for biomolecular data. You describe *what* you want to build (targets) and *how* to compute them (instructions), and DataCooker resolves dependencies, executes steps, and keeps intermediate results in a cache. The repository also ships ready-made pipelines for A3M/FASTA parsing, mmCIF/CCD processing, LMDB packaging, and graph splits used by downstream training jobs.

## Concepts at a glance
- Recipes are declared with a small DSL (`RecipeBook`) that names targets and their inputs.
- Execution is handled by the `Cooker`, which walks dependencies and persists intermediate results in a `ParsingCache`.
- You can mix built-in instructions from `pipelines/instructions/` with your own functions to assemble new pipelines quickly.

## What you can build
- Parse alignments, mmCIF/CCD files, and package graph data into LMDBs using the provided recipe books.
- Create lightweight custom pipelines by combining targets, instructions, and transforms tailored to your data.

## Where to go next
- [Getting Started](getting-started.md): install from GitHub and run a first recipe.
- [Tutorials](tutorials.md): minimal custom recipe and packaged pipeline examples.
- [Concepts](concepts.md): dependency resolution and cache design details.
- [API](api.md): Python API reference for `datacooker`.
