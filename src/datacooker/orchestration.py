"""Backward-compatible orchestration wrappers built on top of v2 runners."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from datacooker.processing import BatchProcessReport
from datacooker.protocols import (
    ConvertFunc,
    DeserializeFunc,
    ProjectFunc,
    TransformFunc,
)
from datacooker.readers import ReaderHooks
from datacooker.recipe import RecipeBook
from datacooker.runners import run_lmdb_extract, run_recipe, run_recipe_batch
from datacooker.writers import WriterHooks


def run_workflow(
    recipe: RecipeBook | str | Path,
    *,
    inputs: Mapping[str, Any],
    transform_func: TransformFunc | None = None,
    key_transform: TransformFunc | None = None,
    project_func: ProjectFunc | None = None,
    writer: WriterHooks | None = None,
    output_path: Path | None = None,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Compatibility wrapper for the v2 ``run_recipe`` runner."""
    reader_hooks = ReaderHooks(key_transform=key_transform or transform_func)
    writer_hooks = WriterHooks(materializer=project_func)
    if writer is not None:
        writer_hooks = WriterHooks(
            serializer=writer.serializer,
            materializer=writer.materializer or writer_hooks.materializer,
        )
    return run_recipe(
        recipe,
        inputs=inputs,
        reader=reader_hooks,
        writer=writer_hooks,
        output_path=output_path,
        **extra_kwargs,
    )


def run_parallel_workflow(
    split_recipe: RecipeBook | str | Path,
    recipe: RecipeBook | str | Path,
    *,
    inputs: Mapping[str, Any],
    split_output_name: str = "data_list",
    transform_func: TransformFunc | None = None,
    key_transform: TransformFunc | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    test_run: bool = True,
    node_rank: int | None = None,
    node_count: int | None = None,
    **extra_kwargs: Any,
) -> BatchProcessReport:
    """Compatibility wrapper for the v2 ``run_recipe_batch`` runner."""
    return run_recipe_batch(
        split_recipe,
        recipe,
        inputs=inputs,
        reader=ReaderHooks(key_transform=key_transform or transform_func),
        split_output_name=split_output_name,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        test_run=test_run,
        node_rank=node_rank,
        node_count=node_count,
        **extra_kwargs,
    )


def extract_lmdb_workflow(
    env_path: Path,
    recipe: RecipeBook | str | Path,
    *,
    deserialize: DeserializeFunc,
    inputs: Mapping[str, Any] | None = None,
    metadata_recipe: Path | None = None,
    metadata_input: Mapping[str, Any] | None = None,
    convert_func: ConvertFunc | None = None,
    transform_func: TransformFunc | None = None,
    key_transform: TransformFunc | None = None,
    merge_recipe: Path | None = None,
    merge_inputs: Mapping[str, Any] | None = None,
    merge_input_name: str = "data_dict",
    chunk_size: int = 100,
    n_jobs: int = -1,
    test_run: bool = True,
    project_func: ProjectFunc | None = None,
    writer: WriterHooks | None = None,
    output_path: Path | None = None,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Compatibility wrapper for the v2 ``run_lmdb_extract`` runner."""
    reader_hooks = ReaderHooks(
        deserializer=deserialize,
        adapter=convert_func,
        key_transform=key_transform or transform_func,
    )
    writer_hooks = WriterHooks(materializer=project_func)
    if writer is not None:
        writer_hooks = WriterHooks(
            serializer=writer.serializer,
            materializer=writer.materializer or writer_hooks.materializer,
        )
    return run_lmdb_extract(
        env_path,
        recipe,
        reader=reader_hooks,
        writer=writer_hooks,
        inputs=inputs,
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        merge_recipe=merge_recipe,
        merge_inputs=merge_inputs,
        merge_input_name=merge_input_name,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        test_run=test_run,
        output_path=output_path,
        **extra_kwargs,
    )
