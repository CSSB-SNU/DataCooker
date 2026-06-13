"""High-level workflow runners composed from the lower-level DataCooker APIs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from datacooker.api import parse_dict
from datacooker.lmdb import extract_lmdb_records
from datacooker.processing import BatchProcessReport, parallel_process_report
from datacooker.protocols import (
    ConvertFunc,
    DeserializeFunc,
    ProjectFunc,
    TransformFunc,
)
from datacooker.recipe import RecipeBook


def run_workflow(
    recipe: RecipeBook | str | Path,
    *,
    inputs: Mapping[str, Any],
    transform_func: TransformFunc | None = None,
    project_func: ProjectFunc | None = None,
    output_path: Path | None = None,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Execute a recipe from prepared inputs and optionally project the result."""
    results = parse_dict(
        recipe,
        dict(inputs),
        transform_func=transform_func,
        **extra_kwargs,
    )
    _project_results(results, project_func=project_func, output_path=output_path)
    return results


def run_parallel_workflow(
    split_recipe: RecipeBook | str | Path,
    recipe: RecipeBook | str | Path,
    *,
    inputs: Mapping[str, Any],
    split_output_name: str = "data_list",
    transform_func: TransformFunc | None = None,
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
        transform_func=transform_func,
    )
    data_list = split_results.get(split_output_name)
    normalized_data_list = _normalize_data_list(
        data_list,
        split_output_name=split_output_name,
    )
    return parallel_process_report(
        normalized_data_list,
        inputs=inputs,
        recipe=Path(recipe) if isinstance(recipe, str) else recipe,
        transform_func=transform_func,
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
    merge_recipe: Path | None = None,
    merge_inputs: Mapping[str, Any] | None = None,
    merge_input_name: str = "data_dict",
    chunk_size: int = 100,
    n_jobs: int = -1,
    test_run: bool = True,
    project_func: ProjectFunc | None = None,
    output_path: Path | None = None,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Extract recipe outputs from an LMDB and optionally project the result."""
    results = extract_lmdb_records(
        env_path,
        Path(recipe) if isinstance(recipe, str) else recipe,
        deserialize=deserialize,
        inputs=inputs,
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        convert_func=convert_func,
        transform_func=transform_func,
        merge_recipe=merge_recipe,
        merge_inputs=merge_inputs,
        merge_input_name=merge_input_name,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        test_run=test_run,
        **extra_kwargs,
    )
    _project_results(results, project_func=project_func, output_path=output_path)
    return results


def _normalize_data_list(
    data_list: Any,
    *,
    split_output_name: str,
) -> list[dict[str, Any]]:
    if data_list is None:
        msg = f"Expected '{split_output_name}' in split workflow output."
        raise KeyError(msg)
    if isinstance(data_list, (str, bytes, bytearray)) or not isinstance(data_list, Sequence):
        msg = (
            f"Expected '{split_output_name}' to be a sequence of mappings, "
            f"got {type(data_list).__name__}."
        )
        raise TypeError(msg)

    normalized: list[dict[str, Any]] = []
    for index, item in enumerate(data_list):
        if not isinstance(item, Mapping):
            msg = (
                f"Expected '{split_output_name}[{index}]' to be a mapping, "
                f"got {type(item).__name__}."
            )
            raise TypeError(msg)
        normalized.append(dict(item))
    return normalized


def _project_results(
    results: dict[str, Any],
    *,
    project_func: ProjectFunc | None,
    output_path: Path | None,
) -> None:
    if project_func is None:
        return
    if output_path is None:
        msg = "output_path is required when project_func is provided."
        raise ValueError(msg)
    project_func(data=results, output_path=output_path)
