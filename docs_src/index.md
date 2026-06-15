# DataCooker (Concepts)

DataCooker is a small static-workflow engine for biomolecular data. You declare
*what* you want to build (**targets**) and *how* to compute it
(**instructions**) with **inputs** that reference other targets; the engine
resolves dependencies, executes the steps in order, and returns the requested
outputs.

## Core concepts

- **RecipeBook** — an ordered set of steps. Each step declares its output
  targets, an instruction callable, and the inputs (args / kwargs / params)
  that feed it. Add steps with `add()` (canonical) or `step()`.
- **Recipe / execution graph** — a RecipeBook resolves to a dependency graph
  (DAG) over named targets. Cycles and missing targets are reported as errors.
- **ExecutionContext** — the keyed value store that holds seeded inputs and
  intermediate/final results during a run.
- **Cooker** — the executor that resolves dependencies, expands wildcard inputs,
  and runs instructions.
- **Reader / Writer hooks** — pluggable I/O around a recipe:
  `ReaderHooks` (loader, key transform) turn a file into seed inputs;
  `WriterHooks` (serializer, materializer) turn the result into stored bytes.

## Running a recipe

The public entry points all live on the top-level `datacooker` package:

| Function | Use it when |
| --- | --- |
| `execute(book, inputs, targets=...)` | you already have an in-memory input dict |
| `parse_dict(book, datadict, targets=...)` | same, with the parsing/IO conventions |
| `parse_file(recipe, file_path, loader=..., targets=...)` | a file should be loaded into inputs first |
| `run_recipe(recipe, inputs=..., reader=..., writer=...)` | you want reader/writer hooks applied |

## LMDB workflows

For large datasets DataCooker ships a config-driven LMDB pipeline
(`python -m datacooker.cli.lmdb`) that scans a directory, runs a recipe per
file, and writes one record per entry. See
[Build an LMDB](tutorials/build_lmdb.md).

## Next steps

- [Getting Started](getting-started.md) — install and run a first recipe.
- [Tutorials](tutorials.md) — write a recipe, then build an LMDB.
- [API](api.md) — the public `datacooker` surface.
