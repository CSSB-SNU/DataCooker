"""High-level runners composed from the lower-level DataCooker APIs."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from datacooker.api import execute, parse_dict
from datacooker.lmdb import extract_lmdb_records
from datacooker.processing import BatchProcessReport, parallel_process_report
from datacooker.protocols import TransformFunc
from datacooker.readers import ReaderHooks, normalize_input_sequence
from datacooker.recipe import RecipeBook
from datacooker.writers import WriterHooks, write_output


def run_recipe(
    recipe: RecipeBook | str | Path,
    *,
    inputs: Mapping[str, Any],
    reader: ReaderHooks | None = None,
    writer: WriterHooks | None = None,
    key_transform: TransformFunc | None = None,
    targets: list[str] | str | None = None,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Execute a recipe from prepared inputs and optionally materialize the result."""
    output_path = extra_kwargs.pop("output_path", None)
    results = execute(
        recipe,
        dict(inputs),
        key_transform=key_transform,
        transform_func=None if reader is None else reader.key_transform,
        targets=targets,
        **extra_kwargs,
    )
    write_output(results, writer=writer, output_path=output_path)
    return results


def run_recipe_batch(
    split_recipe: RecipeBook | str | Path,
    recipe: RecipeBook | str | Path,
    *,
    inputs: Mapping[str, Any],
    reader: ReaderHooks | None = None,
    split_output_name: str = "data_list",
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    test_run: bool = True,
    node_rank: int | None = None,
    node_count: int | None = None,
    **extra_kwargs: Any,
) -> BatchProcessReport:
    """Expand work items with one recipe and process each item with another."""
    split_results = parse_dict(
        split_recipe,
        dict(inputs),
        key_transform=None if reader is None else reader.key_transform,
    )
    data_list = normalize_input_sequence(
        split_results.get(split_output_name),
        output_name=split_output_name,
    )
    return parallel_process_report(
        data_list,
        inputs=inputs,
        recipe=Path(recipe) if isinstance(recipe, str) else recipe,
        key_transform=None if reader is None else reader.key_transform,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        test_run=test_run,
        node_rank=node_rank,
        node_count=node_count,
        **extra_kwargs,
    )


def run_lmdb_extract(
    env_path: Path,
    recipe: RecipeBook | str | Path,
    *,
    reader: ReaderHooks,
    writer: WriterHooks | None = None,
    inputs: Mapping[str, Any] | None = None,
    metadata_recipe: Path | None = None,
    metadata_input: Mapping[str, Any] | None = None,
    merge_recipe: Path | None = None,
    merge_inputs: Mapping[str, Any] | None = None,
    merge_input_name: str = "data_dict",
    chunk_size: int = 100,
    n_jobs: int = -1,
    test_run: bool = True,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Extract recipe outputs from an LMDB and optionally materialize the result."""
    output_path = extra_kwargs.pop("output_path", None)
    if reader.deserializer is None:
        msg = "run_lmdb_extract requires reader.deserializer."
        raise ValueError(msg)
    results = extract_lmdb_records(
        env_path,
        Path(recipe) if isinstance(recipe, str) else recipe,
        reader=reader,
        inputs=inputs,
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        merge_recipe=merge_recipe,
        merge_inputs=merge_inputs,
        merge_input_name=merge_input_name,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        test_run=test_run,
        **extra_kwargs,
    )
    write_output(results, writer=writer, output_path=output_path)
    return results
