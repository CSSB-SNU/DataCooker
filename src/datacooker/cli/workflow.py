"""Config-driven workflow execution entrypoints."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click

from datacooker.cli._common import resolve_optional_int
from datacooker.config import load_config
from datacooker.orchestration import (
    extract_lmdb_workflow,
    run_parallel_workflow,
    run_workflow,
)


@click.group()
def cli() -> None:
    """Run config-driven DataCooker workflows."""


@cli.command("run")
@click.argument("config_path", type=click.Path(exists=True, path_type=Path))
def run_command(config_path: Path) -> None:
    """Run a recipe against prepared inputs and optionally project the result."""
    config = load_config(config_path)
    recipe = _pop_recipe_path(config, "recipe", "recipe_path")
    inputs = _pop_inputs(config)
    output_path = _pop_optional_path(config, "output_path", "output_data_path")
    project_func = config.pop("project_func", None)

    results = run_workflow(
        recipe,
        inputs=inputs,
        transform_func=config.pop("transform_func", None),
        project_func=project_func,
        output_path=output_path,
        **config,
    )
    click.echo(f"produced {len(results)} outputs")


@cli.command("parallel-run")
@click.argument("config_path", type=click.Path(exists=True, path_type=Path))
@click.option("--node-rank", type=int, default=None)
@click.option("--node-count", type=int, default=None)
def parallel_run_command(
    config_path: Path,
    node_rank: int | None,
    node_count: int | None,
) -> None:
    """Expand work items with one recipe and process them with another."""
    config = load_config(config_path)
    split_recipe = _pop_recipe_path(config, "split_recipe", "split_recipe_path")
    recipe = _pop_recipe_path(config, "recipe", "recipe_path")
    inputs = _pop_inputs(config)

    resolved_node_rank = resolve_optional_int(
        node_rank,
        config,
        ("node_rank", "rank", "shard_idx"),
    )
    resolved_node_count = resolve_optional_int(
        node_count,
        config,
        ("node_count", "world_size", "n_shards"),
    )

    results = run_parallel_workflow(
        split_recipe,
        recipe,
        inputs=inputs,
        split_output_name=str(config.pop("split_output_name", "data_list")),
        transform_func=config.pop("transform_func", None),
        chunk_size=int(config.pop("chunk_size", 10_000)),
        n_jobs=int(config.pop("n_jobs", -1)),
        test_run=bool(config.pop("test_run", True)),
        node_rank=resolved_node_rank,
        node_count=resolved_node_count,
        **config,
    )
    click.echo(f"processed {results.attempted} items with {results.failed} errors")


@cli.command("extract-lmdb")
@click.argument("config_path", type=click.Path(exists=True, path_type=Path))
def extract_lmdb_command(config_path: Path) -> None:
    """Extract recipe outputs from all LMDB records."""
    config = load_config(config_path)
    env_path = _pop_path(config, "env_path", "db_path")
    recipe = _pop_recipe_path(config, "extract_recipe", "extract_recipe_path", "recipe")
    output_path = _pop_optional_path(config, "output_path", "output_data_path")
    project_func = config.pop("project_func", None)
    deserialize = config.pop("deserialize", None)
    if deserialize is None:
        msg = "extract-lmdb config requires a 'deserialize' callable."
        raise click.ClickException(msg)

    inputs = dict(config.pop("inputs", {}))
    additional_inputs = config.pop("additional_inputs", None)
    if additional_inputs is not None:
        inputs.update(additional_inputs)

    merge_recipe = _pop_optional_path(config, "merge_recipe", "merge_recipe_path")
    metadata_recipe = _pop_optional_path(config, "metadata_recipe")

    results = extract_lmdb_workflow(
        env_path,
        recipe,
        deserialize=deserialize,
        inputs=inputs,
        metadata_recipe=metadata_recipe,
        metadata_input=config.pop("metadata_input", None),
        convert_func=config.pop("convert_func", None),
        transform_func=config.pop("transform_func", None),
        merge_recipe=merge_recipe,
        merge_inputs=config.pop("merge_inputs", None),
        merge_input_name=str(config.pop("merge_input_name", "data_dict")),
        chunk_size=int(config.pop("chunk_size", 100)),
        n_jobs=int(config.pop("n_jobs", -1)),
        test_run=bool(config.pop("test_run", True)),
        project_func=project_func,
        output_path=output_path,
        **config,
    )
    click.echo(f"extracted {len(results)} records")


def _pop_inputs(config: dict[str, Any]) -> dict[str, Any]:
    inputs = config.pop("inputs", None)
    if inputs is None:
        msg = "workflow config requires an 'inputs' mapping."
        raise click.ClickException(msg)
    if not isinstance(inputs, dict):
        msg = f"'inputs' must be a mapping, got {type(inputs).__name__}."
        raise click.ClickException(msg)
    return inputs


def _pop_path(config: dict[str, Any], *names: str) -> Path:
    for name in names:
        value = config.pop(name, None)
        if value is not None:
            return Path(value) if isinstance(value, str) else value
    joined = ", ".join(names)
    msg = f"workflow config requires one of: {joined}."
    raise click.ClickException(msg)


def _pop_optional_path(config: dict[str, Any], *names: str) -> Path | None:
    for name in names:
        if name not in config:
            continue
        value = config.pop(name)
        if value is None:
            return None
        return Path(value) if isinstance(value, str) else value
    return None


def _pop_recipe_path(config: dict[str, Any], *names: str) -> Path:
    return _pop_path(config, *names)


def main() -> None:
    """Entry point for ``python -m datacooker.cli.workflow``."""
    cli()


if __name__ == "__main__":
    main()
