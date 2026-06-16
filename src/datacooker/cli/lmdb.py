"""Config-driven LMDB workflow entrypoints."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import click

from datacooker import scan_paths, shard_items
from datacooker.cli._common import coerce_path, pop_config_value, with_shard_suffix
from datacooker.config import load_config
from datacooker.lmdb import (
    build_lmdb,
    count_lmdb_entries,
    merge_lmdb_shards,
    rebuild_lmdb,
)
from datacooker.readers import ReaderHooks
from datacooker.writers import WriterHooks


@click.group()
def cli() -> None:
    """Run reusable LMDB workflows from declarative config files."""


@cli.command("build")
@click.argument("config_path", type=click.Path(exists=True, path_type=Path))
@click.option("--map-size", type=int, default=None, help="Override LMDB map size.")
@click.option("--shard-idx", type=int, default=None, help="0-based shard index.")
@click.option("--n-shards", type=int, default=1, show_default=True)
def build_command(
    config_path: Path,
    map_size: int | None,
    shard_idx: int | None,
    n_shards: int,
) -> None:
    """Build an LMDB from scanned files plus a DataCooker recipe."""
    config = load_config(config_path)
    file_list = pop_config_value(config, "file_list", default=None)
    if file_list is None:
        data_dir = Path(pop_config_value(config, "data_dir"))
        file_pattern = str(pop_config_value(config, "file_pattern", default="*"))
    else:
        pop_config_value(config, "data_dir", default=None)
        pop_config_value(config, "file_pattern", default=None)
    coerce_path(config, "env_path")
    _normalize_recipe_alias(config)
    coerce_path(config, "metadata_recipe")
    config["reader"] = _pop_reader_hooks(config)
    config["writer"] = _pop_writer_hooks(config)

    if file_list is not None:
        with Path(file_list).open() as handle:
            data_list = [Path(line.strip()) for line in handle if line.strip()]
    else:
        data_list = scan_paths(data_dir, pattern=file_pattern)
    if shard_idx is not None:
        if shard_idx < 0 or shard_idx >= n_shards:
            msg = f"Invalid shard index {shard_idx} for {n_shards} shards."
            raise click.ClickException(msg)
        data_list = shard_items(data_list, node_rank=shard_idx, node_count=n_shards)
        env_path = coerce_path(config, "env_path")
        if env_path is None:
            msg = "env_path is required in the LMDB build config."
            raise click.ClickException(msg)
        config["env_path"] = with_shard_suffix(env_path, shard_idx)

    if map_size is not None:
        config["map_size"] = map_size

    report = build_lmdb(*data_list, **config)
    env_path = coerce_path(config, "env_path")
    if env_path is None:
        msg = "env_path is required in the LMDB build config."
        raise click.ClickException(msg)
    click.echo(
        f"built {env_path} "
        f"(attempted={report.attempted}, written={report.written}, "
        f"skipped_existing={report.skipped_existing}, failed={report.failed})"
    )
    click.echo(f"entries={count_lmdb_entries(env_path)}")


@cli.command("rebuild")
@click.argument("config_path", type=click.Path(exists=True, path_type=Path))
@click.option("--map-size", type=int, default=None, help="Override LMDB map size.")
def rebuild_command(config_path: Path, map_size: int | None) -> None:
    """Rebuild an LMDB through another DataCooker recipe."""
    config = load_config(config_path)
    coerce_path(config, "old_env_path")
    coerce_path(config, "new_env_path")
    _normalize_recipe_alias(config)
    coerce_path(config, "metadata_recipe")
    config["reader"] = _pop_reader_hooks(config)
    config["writer"] = _pop_writer_hooks(config)
    if map_size is not None:
        config["map_size"] = map_size

    report = rebuild_lmdb(**config)
    new_env_path = coerce_path(config, "new_env_path")
    if new_env_path is None:
        msg = "new_env_path is required in the LMDB rebuild config."
        raise click.ClickException(msg)
    click.echo(
        f"rebuilt {new_env_path} "
        f"(attempted={report.attempted}, written={report.written}, "
        f"skipped_existing={report.skipped_existing}, failed={report.failed})"
    )
    click.echo(f"entries={count_lmdb_entries(new_env_path)}")


@cli.command("merge")
@click.argument("shard_pattern", type=str)
@click.option(
    "--output",
    required=True,
    type=click.Path(path_type=Path),
    help="Output merged LMDB path.",
)
@click.option("--map-size", type=int, default=int(1e12), show_default=True)
@click.option("--overwrite", is_flag=True, help="Overwrite the output if it exists.")
@click.option(
    "--on-duplicate",
    type=click.Choice(["overwrite", "skip", "error"]),
    default="overwrite",
    show_default=True,
)
def merge_command(
    shard_pattern: str,
    output: Path,
    map_size: int,
    overwrite: bool,
    on_duplicate: str,
) -> None:
    """Merge shard LMDBs matched by a glob pattern."""
    pattern_path = Path(shard_pattern)
    shard_paths = sorted(pattern_path.parent.glob(pattern_path.name))
    if not shard_paths:
        msg = f"No LMDB files matched: {shard_pattern}"
        raise click.ClickException(msg)

    merge_lmdb_shards(
        shard_paths,
        output,
        map_size=map_size,
        overwrite=overwrite,
        on_duplicate=on_duplicate,
    )
    click.echo(f"merged {len(shard_paths)} shards into {output}")
    click.echo(f"entries={count_lmdb_entries(output)}")


@cli.command("count")
@click.argument("env_path", type=click.Path(exists=True, path_type=Path))
def count_command(env_path: Path) -> None:
    """Count the number of records in an LMDB."""
    click.echo(f"{env_path}: {count_lmdb_entries(env_path)}")


def main() -> None:
    """Entry point for ``python -m datacooker.cli.lmdb``."""
    cli()


def _normalize_recipe_alias(config: dict[str, Any]) -> None:
    if "recipe" in config:
        coerce_path(config, "recipe")
        return
    recipe_path = config.pop("recipe_path", None)
    if recipe_path is None:
        return
    config["recipe"] = Path(recipe_path) if isinstance(recipe_path, str) else recipe_path


def _pop_reader_hooks(config: dict[str, Any]) -> ReaderHooks:
    reader_config = config.pop("reader", None)
    return ReaderHooks.from_mapping(
        reader_config if isinstance(reader_config, dict) else None,
        loader=config.pop("loader", None),
        adapter=config.pop("adapter", None),
        deserializer=config.pop("deserializer", None),
        key_transform=config.pop("key_transform", None),
    )


def _pop_writer_hooks(config: dict[str, Any]) -> WriterHooks:
    writer_config = config.pop("writer", None)
    return WriterHooks.from_mapping(
        writer_config if isinstance(writer_config, dict) else None,
        serializer=config.pop("serializer", None),
        materializer=config.pop("materializer", None),
    )


if __name__ == "__main__":
    main()
