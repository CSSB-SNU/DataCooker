# API

The public surface is re-exported from the top-level `datacooker` package
(`import datacooker`). Items are grouped by role below.

## Recipes

- `RecipeBook` — build a workflow with `add(targets, instruction, inputs)` or
  `step(outputs, instruction, *, args, kwargs, params, ...)`; introspect with
  `describe`, `execution_order`, `validate`, `to_dot`, `to_mermaid`.
- `Recipe`, `Inputs`, `Variable`, `variable` — recipe / step value model.

## Execution

- `execute(book, inputs, *, targets=None, ...)` — run a recipe over an
  in-memory input dict and return the requested targets.
- `parse_dict(recipe, datadict, *, targets=None, ...)` — run against a prepared
  dict using the parsing conventions.
- `parse_file(recipe, file_path, loader=None, inputs=None, *, targets=None, ...)`
  — load a file into inputs (via `loader`) and run the recipe.
- `run_recipe(recipe, *, inputs, reader=None, writer=None, ...)` /
  `run_recipe_batch(...)` — run with reader/writer hooks applied (single / batch).
- `Cooker`, `ExecutionContext` — the executor and its keyed value store.
- `describe`, `visualize`, `load_recipe`, `load_inputs` — graph inspection and
  loading helpers.

## I/O hooks and payloads

- `ReaderHooks`, `WriterHooks`, `InputHooks`, `PayloadHooks` — pluggable I/O
  around a recipe (loader / key transform / serializer / materializer).
- `encode_output`, `decode_payload`, `write_output` — encode results, decode
  payloads, and materialise outputs.

## LMDB workflows

- `run_lmdb_extract` — extract records from an LMDB through a recipe.
- The `datacooker.lmdb` subpackage provides `build_lmdb`, `rebuild_lmdb`,
  `merge_lmdb_shards`, `count_lmdb_entries`, `extract_lmdb_keys`,
  `read_lmdb_raw`, `default_lmdb_key`, and `filter_pending_lmdb_paths`.
- The `datacooker.cli.lmdb` CLI exposes `build`, `rebuild`, `merge`, and
  `count` (see [Build an LMDB](tutorials/build_lmdb.md)).

## Generic instructions and helpers

- `extract_float_single`, `key_stack`, `merge_dict`, `single_value_instruction`
  — dependency-free transform instructions for recipes.
- `scan_paths`, `shard_items`, `resolve_node_config` — file discovery and
  sharding utilities.
- `resolve_object`, `dot_path` — resolve dotted import paths used in configs.
- `parallel_process_report`, `BatchProcessReport`, `BatchProcessResult` —
  parallel processing and its reports.

## Type aliases

`LoadFunc`, `SerializeFunc`, `DeserializeFunc`, `TransformFunc`, `ConvertFunc`,
`ProjectFunc`, `KeyFunc`.

## Errors

`DataCookerError` (base), `RecipeError`, `InvalidRecipeError`, `CycleError`,
`DuplicateTargetError`, `MissingDependencyError`, `MissingTargetError`,
`UnknownTargetError`, `InstructionOutputError`, `StepExecutionError`.
