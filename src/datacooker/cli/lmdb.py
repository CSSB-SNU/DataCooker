"""Config-driven LMDB workflow entrypoints."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import click

from datacooker import scan_paths
from datacooker.cli._common import coerce_path, pop_config_value, with_shard_suffix
from datacooker.config import load_config
from datacooker.config.runtime import activate_environment, canonical_inputs
from datacooker.lmdb import (
    build_lmdb,
    count_lmdb_entries,
    default_lmdb_key,
    merge_lmdb_shards,
    rebuild_lmdb,
    rekey_lmdb,
    suffix_key_map,
)
from datacooker.readers import ReaderHooks
from datacooker.utils.importing import resolve_object
from datacooker.utils.sharding import shard_items_by_key
from datacooker.writers import WriterHooks

if TYPE_CHECKING:
    from datacooker.lmdb.core import DuplicatePolicy


@click.group()
def cli() -> None:
    """Run reusable LMDB workflows from declarative config files."""


@cli.command("build")
@click.argument("config_path", type=click.Path(exists=True, path_type=Path))
@click.option("--map-size", type=int, default=None, help="Override LMDB map size.")
@click.option("--shard-idx", type=int, default=None, help="0-based shard index.")
@click.option("--n-shards", type=int, default=1, show_default=True)
@click.option("--file-list", "file_list_override",
              type=click.Path(exists=True, path_type=Path), default=None,
              help="Process only the files in this list (e.g. a planner size tier).")
@click.option("--output", type=click.Path(path_type=Path), default=None,
              help="Override env_path (e.g. per tier-shard output).")
@click.option("--n-jobs", "n_jobs", type=int, default=None,
              help="Maximum concurrent workers, bounded by allocated CPUs.")
@click.option("--chunk-size", "chunk_size", type=int, default=None,
              help="Results buffered before writing (default 10k lets it hoard tens of "
                   "GB and OOM); planner sets ~a few x n_jobs.")
@click.option("--skip-existing", is_flag=True, default=False,
              help="Resume: skip keys already in the output (idempotent re-run).")
@click.option("--ray-batch", "ray_batch", type=int, default=None,
              help="Items per Ray task (default 1). Rarely needed; the memory gate bounds concurrency.")
def build_command(
    config_path: Path,
    map_size: int | None,
    shard_idx: int | None,
    n_shards: int,
    file_list_override: Path | None,
    output: Path | None,
    n_jobs: int | None,
    chunk_size: int | None,
    skip_existing: bool,
    ray_batch: int | None,
) -> None:
    """Build an LMDB from scanned files plus a DataCooker recipe."""
    config = load_config(config_path)
    activate_environment(config)
    if output is not None:
        config["env_path"] = output
    if n_jobs is not None:
        config["n_jobs"] = n_jobs
    if chunk_size is not None:
        config["chunk_size"] = chunk_size
    if ray_batch is not None:
        config["ray_batch"] = ray_batch
    if skip_existing:
        config["skip_existing"] = True
    # `keyed` is a planner-only hint (items are LMDB keys, not files); it never
    # reaches build_lmdb.
    pop_config_value(config, "keyed", default=None)
    # Always pop the config's file_list (avoid leaking it into build_lmdb kwargs);
    # an explicit --file-list overrides it (e.g. a planner size tier).
    config_file_list = pop_config_value(config, "file_list", default=None)
    file_list = file_list_override or config_file_list
    data_dir: Path | None = None
    file_pattern = "*"
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
        data_list = scan_paths(cast("Path", data_dir), pattern=file_pattern)
    key_func = config.get("key_builder", default_lmdb_key)
    if isinstance(key_func, str):
        key_func = resolve_object(key_func)
    data_list, _duplicates = canonical_inputs(data_list, key_func, config.get("duplicate_input_policy", "error"))
    if shard_idx is not None:
        if shard_idx < 0 or shard_idx >= n_shards:
            msg = f"Invalid shard index {shard_idx} for {n_shards} shards."
            raise click.ClickException(msg)
        # Shard by a stable hash of each file's LMDB key (not list position), so
        # the same record always lands in the same shard regardless of scan
        # order / which run builds which shard. Positional striping over an
        # unsorted scan let the same key land in different shards across re-runs,
        # producing cross-shard duplicate keys in the merged output.
        key_builder = config.get("key_builder", default_lmdb_key)
        key_func = (
            resolve_object(key_builder) if isinstance(key_builder, str) else key_builder
        )
        data_list = shard_items_by_key(
            data_list,
            key_func=key_func,
            node_rank=shard_idx,
            node_count=n_shards,
        )
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
    report_path = Path(str(env_path) + ".build-report.json")
    temporary = report_path.with_suffix(".tmp")
    temporary.write_text(json.dumps(asdict(report), indent=2))
    temporary.replace(report_path)


@cli.command("rebuild")
@click.argument("config_path", type=click.Path(exists=True, path_type=Path))
@click.option("--map-size", type=int, default=None, help="Override LMDB map size.")
@click.option("--shard-idx", type=int, default=None, help="0-based shard index.")
@click.option("--n-shards", type=int, default=1, show_default=True)
@click.option("--key-list", type=click.Path(exists=True, path_type=Path), default=None,
              help="Process only the keys in this file (one per line), e.g. a size tier.")
@click.option("--output", type=click.Path(path_type=Path), default=None,
              help="Override new_env_path (e.g. per tier-shard output).")
@click.option("--n-jobs", "n_jobs", type=int, default=None,
              help="Maximum concurrent workers, bounded by allocated CPUs.")
@click.option("--chunk-size", "chunk_size", type=int, default=None,
              help="Results buffered before writing; the 10k default lets the drain "
                   "hoard tens of GB and OOM. Planner sets ~a few x n_jobs.")
@click.option("--ray-batch", "ray_batch", type=int, default=None,
              help="Items per Ray task (default 1). Rarely needed; the memory gate bounds concurrency.")
def rebuild_command(
    config_path: Path,
    map_size: int | None,
    shard_idx: int | None,
    n_shards: int,
    key_list: Path | None,
    output: Path | None,
    n_jobs: int | None,
    chunk_size: int | None,
    ray_batch: int | None,
) -> None:
    """Rebuild an LMDB through another DataCooker recipe."""
    config = load_config(config_path)
    activate_environment(config)
    coerce_path(config, "old_env_path")
    coerce_path(config, "new_env_path")
    _normalize_recipe_alias(config)
    coerce_path(config, "metadata_recipe")
    config["reader"] = _pop_reader_hooks(config)
    config["writer"] = _pop_writer_hooks(config)
    if map_size is not None:
        config["map_size"] = map_size
    if output is not None:
        config["new_env_path"] = output
    if key_list is not None:
        config["key_list"] = key_list
    if n_jobs is not None:
        config["n_jobs"] = n_jobs
    if chunk_size is not None:
        config["chunk_size"] = chunk_size
    if ray_batch is not None:
        config["ray_batch"] = ray_batch

    if shard_idx is not None:
        if shard_idx < 0 or shard_idx >= n_shards:
            msg = f"Invalid shard index {shard_idx} for {n_shards} shards."
            raise click.ClickException(msg)
        new_env_path = coerce_path(config, "new_env_path")
        if new_env_path is None:
            msg = "new_env_path is required in the LMDB rebuild config."
            raise click.ClickException(msg)
        config["new_env_path"] = with_shard_suffix(new_env_path, shard_idx)
        config["shard_idx"] = shard_idx
        config["n_shards"] = n_shards

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
    report_path = Path(str(new_env_path) + ".build-report.json")
    temporary = report_path.with_suffix(".tmp")
    temporary.write_text(json.dumps(asdict(report), indent=2))
    temporary.replace(report_path)


@cli.command("rekey")
@click.argument("config_path", type=click.Path(exists=True, path_type=Path))
@click.option("--map-size", type=int, default=None, help="Override LMDB map size.")
def rekey_command(config_path: Path, map_size: int | None) -> None:
    """Copy an LMDB to a new database, transforming keys but copying values as-is.

    The config needs ``old_env_path`` and ``new_env_path`` plus exactly one key
    rule: ``suffix`` (append a constant string) or ``key_map`` (dotted import
    path to a ``(str) -> str | None`` callable). Optional keys: ``map_size``,
    ``commit_every``, ``skip_existing``, ``on_duplicate``, ``append``.
    """
    config = load_config(config_path)
    activate_environment(config)
    old_env_path = coerce_path(config, "old_env_path")
    new_env_path = coerce_path(config, "new_env_path")
    if old_env_path is None or new_env_path is None:
        msg = "old_env_path and new_env_path are required in the rekey config."
        raise click.ClickException(msg)

    kwargs: dict[str, Any] = {
        "old_env_path": old_env_path,
        "new_env_path": new_env_path,
        "key_map": _resolve_key_map(config),
    }
    for option in (
        "map_size",
        "commit_every",
        "skip_existing",
        "on_duplicate",
        "append",
    ):
        if option in config:
            kwargs[option] = config[option]
    if map_size is not None:
        kwargs["map_size"] = map_size

    report = rekey_lmdb(**kwargs)
    click.echo(
        f"rekeyed {new_env_path} "
        f"(attempted={report.attempted}, written={report.written}, "
        f"skipped_existing={report.skipped_existing}, "
        f"duplicate_keys={report.duplicate_keys}, "
        f"skipped_empty={report.skipped_empty})"
    )
    click.echo(f"entries={count_lmdb_entries(new_env_path)}")


def required_map_size(shard_paths: list[Path], headroom: float = 1.25) -> int:
    """Map size a merge of ``shard_paths`` needs: their data files' total size plus headroom."""
    total = sum((Path(p) / "data.mdb").stat().st_size for p in shard_paths if (Path(p) / "data.mdb").is_file())
    return int(total * headroom)


@cli.command("merge")
@click.argument("shard_pattern", type=str, required=False)
@click.option("--shard-list", type=click.Path(exists=True, path_type=Path),
              help="Exact shard paths, one per line; mutually exclusive with the glob.")
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
    type=click.Choice(["overwrite", "skip", "error", "union"]),
    default="overwrite",
    show_default=True,
)
@click.option(
    "--value-deserializer",
    type=str,
    default=None,
    help="Dotted path to a bytes->dict decoder (required for --on-duplicate union).",
)
@click.option(
    "--value-serializer",
    type=str,
    default=None,
    help="Dotted path to a dict->bytes encoder (required for --on-duplicate union).",
)
def merge_command(
    shard_pattern: str | None,
    shard_list: Path | None,
    output: Path,
    map_size: int,
    overwrite: bool,
    on_duplicate: str,
    value_deserializer: str | None,
    value_serializer: str | None,
) -> None:
    """Merge shard LMDBs matched by a glob pattern."""
    if (shard_pattern is None) == (shard_list is None):
        msg = "Provide exactly one shard glob or --shard-list."
        raise click.ClickException(msg)
    if shard_list is not None:
        shard_paths = [Path(line) for line in shard_list.read_text().splitlines() if line]
        missing = [str(path) for path in shard_paths if not (path / "data.mdb").is_file()]
        if missing:
            msg = "Missing planned shards: " + ", ".join(missing)
            raise click.ClickException(msg)
    else:
        pattern_path = Path(cast("str", shard_pattern))
        shard_paths = sorted(pattern_path.parent.glob(pattern_path.name))
        if not shard_paths:
            msg = f"No LMDB files matched: {shard_pattern}"
            raise click.ClickException(msg)

    deserializer = (
        resolve_object(value_deserializer) if value_deserializer else None
    )
    serializer = resolve_object(value_serializer) if value_serializer else None

    # The map is a ceiling, not an allocation: never below what the shards already hold.
    map_size = max(map_size, required_map_size(shard_paths))
    merge_lmdb_shards(
        shard_paths,
        output,
        map_size=map_size,
        overwrite=overwrite,
        on_duplicate=cast("DuplicatePolicy", on_duplicate),
        value_deserializer=deserializer,
        value_serializer=serializer,
    )
    storage = "shard_lmdb" if (output / "shards.json").is_file() else "merged LMDB"
    click.echo(f"published {len(shard_paths)} shards as {storage}: {output}")
    click.echo(f"entries={count_lmdb_entries(output)}")


@cli.command("count")
@click.argument("env_path", type=click.Path(exists=True, path_type=Path))
def count_command(env_path: Path) -> None:
    """Count the number of records in an LMDB."""
    click.echo(f"{env_path}: {count_lmdb_entries(env_path)}")


@cli.command("index")
@click.argument("db_path", type=click.Path(exists=True, path_type=Path))
@click.option("--schema", default="unknown", help="schema class A-F for meta.json")
@click.option("--key-convention", default="unknown")
@click.option("--codec", default="biomol_zstd")
@click.option("--source", default=None, help="lineage: DB this was built from")
@click.option("--recipe", "recipe", default=None)
@click.option("--expansion", type=float, default=None, help="measured E for the planner")
def index_command(
    db_path: Path,
    schema: str,
    key_convention: str,
    codec: str,
    source: str | None,
    recipe: str | None,
    expansion: float | None,
) -> None:
    """Write a DB's per-item size sidecar + metadata descriptor.

    Mandatory post-step of every build.
    """
    from datacooker.lmdb.index import DbMeta, write_index, write_meta

    stats = write_index(db_path)
    write_meta(db_path, DbMeta(
        schema=schema, key_convention=key_convention, codec=codec,
        entries=stats.entries, total_bytes=stats.total_bytes,
        max_bytes=stats.max_bytes, expansion=expansion, source=source, recipe=recipe,
    ))
    click.echo(f"indexed {db_path}: entries={stats.entries} "
               f"max_bytes={stats.max_bytes} n_gt_1GB={stats.n_gt_1gb}")


@cli.command("pipeline")
@click.argument("engine_config", type=click.Path(exists=True, path_type=Path))
@click.option("--workdir", type=click.Path(path_type=Path), required=True,
              help="scratch for item-lists, sbatch scripts, logs")
@click.option("--expansion", type=float, default=100.0)
@click.option("--include-keys", type=click.Path(exists=True, path_type=Path), default=None,
              help="rebuild only: restrict to this key subset (e.g. a date-scoped set)")
@click.option("--schema", default="unknown")
@click.option("--node-mem-gb", type=click.FloatRange(min=1), default=490.0, envvar="DATACOOKER_NODE_MEM_GB")
@click.option("--cores", type=click.IntRange(min=1), default=112, envvar="DATACOOKER_CPUS_PER_NODE")
@click.option("--n-nodes", type=click.IntRange(min=1), default=8, envvar="DATACOOKER_MAX_NODES")
@click.option("--repo", type=click.Path(path_type=Path), default=None)
@click.option("--depends-on", "depends_on", default="",
              help="comma-separated upstream job ids the tier arrays afterok-wait "
                   "on (multi-DB chaining: this DB starts after upstream DBs finish)")
@click.option("--dry-run", is_flag=True, help="plan + write scripts, do not sbatch")
def pipeline_command(
    engine_config: Path,
    workdir: Path,
    expansion: float,
    include_keys: Path | None,
    schema: str,
    node_mem_gb: float,
    cores: int,
    n_nodes: int,
    repo: Path | None,
    depends_on: str,
    dry_run: bool,
) -> None:
    """Plan a build/rebuild and submit chained tier arrays + merge + index.

    Size the items, plan size tiers, submit tier arrays + merge + index -- all
    chained. Op is inferred from the engine config (old_env_path => rebuild;
    env_path + file_list/data_dir => build).
    """
    from datacooker.executors.slurm import SlurmExecutor, wait_for_jobs
    from datacooker.lmdb.index import ItemSize, measure_files, read_index, write_index
    from datacooker.pipeline import PipelineConfig, run
    from datacooker.planning import NodeSpec

    deps = tuple(d.strip() for d in (depends_on or "").split(",") if d.strip())
    if deps and not dry_run:
        wait_for_jobs(deps)
        deps = ()  # Completed jobs can expire from the scheduler before submission.
    raw = load_config(engine_config)
    activate_environment(raw)
    coerce_path(raw, "old_env_path")
    coerce_path(raw, "new_env_path")
    coerce_path(raw, "env_path")

    if raw.get("old_env_path") is not None:            # ---- rebuild ----
        op, source, target = "rebuild", raw["old_env_path"], raw["new_env_path"]
        sizes = read_index(source)
        if not sizes:
            click.echo(f"source index missing -> building {source}.index.tsv")
            write_index(source)
            sizes = read_index(source)
        if include_keys:
            keep = {ln.strip() for ln in Path(include_keys).read_text().splitlines() if ln.strip()}
            sizes = [s for s in sizes if s.key in keep]
    elif raw.get("env_path") is not None:              # ---- build ----
        op, source, target = "build", None, raw["env_path"]
        file_list = raw.get("file_list")
        if file_list is not None:
            items = [ln.strip() for ln in Path(file_list).read_text().splitlines() if ln.strip()]
            if raw.get("keyed"):
                # Keyed build: file_list holds LMDB keys (seqids / chain-ids), not
                # files on disk -- the recipe resolves each item by lookup, so there
                # is no input st_size. Size uniformly and let the planner balance by
                # count (template Phase 3/4). Closes the key-list sizing gap.
                sizes = [ItemSize(k, 1) for k in items]
            else:
                paths = [Path(k) for k in items]
                sizes = measure_files(paths)
        else:
            files = list(scan_paths(Path(raw["data_dir"]), pattern=str(raw.get("file_pattern", "*"))))
            click.echo(f"[planner] discovered {len(files)} files; collecting sizes")
            sizes = measure_files(files)
    else:
        msg = "engine config needs old_env_path (rebuild) or env_path+file source (build)."
        raise click.ClickException(msg)

    if op == "build":
        key_func = raw.get("key_builder", default_lmdb_key)
        if isinstance(key_func, str):
            key_func = resolve_object(key_func)
        selected, duplicates = canonical_inputs([Path(item.key) for item in sizes], key_func,
                                                 raw.get("duplicate_input_policy", "error"))
        keep = {str(p) for p in selected}
        sizes = [item for item in sizes if item.key in keep]
        Path(workdir).mkdir(parents=True, exist_ok=True)
        (Path(workdir) / "input-selection.json").write_text(json.dumps(duplicates, indent=2))
    pcfg = PipelineConfig(
        op=op, source=source, target=target, engine_config=engine_config,
        expansion=expansion, schema=schema,
        node=NodeSpec(mem_gb=node_mem_gb, cores=cores, n_nodes=n_nodes),
        map_size=int(raw.get("map_size", 300_000_000_000)),
    )
    run_workdir = Path(workdir) / pcfg.run_id
    execu = SlurmExecutor(
        workdir=run_workdir, repo=Path(repo) if repo else Path.cwd(),
        submit=not dry_run,
    )
    handle = run(pcfg, execu, keyfiles_dir=run_workdir / "items", sizes=sizes,
                 depends_on=deps)
    click.echo(f"pipeline submitted; final stage = {handle.name} ({handle.job_id})")
    # Machine-parseable terminal job id so a multi-DB driver can afterok-chain the
    # next config's pipeline onto this one's final (index) stage.
    click.echo(f"PIPELINE_TERMINAL_JOB={handle.job_id or ''}")


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
    config["recipe"] = (
        Path(recipe_path) if isinstance(recipe_path, str) else recipe_path
    )


def _resolve_key_map(config: dict[str, Any]) -> Any:
    """Resolve the rekey ``key_map`` from a config's ``suffix``/``key_map`` field."""
    suffix = config.pop("suffix", None)
    key_map = config.pop("key_map", None)
    if suffix is not None:
        if key_map is not None:
            msg = "Specify only one of 'suffix' or 'key_map' in the rekey config."
            raise click.ClickException(msg)
        return suffix_key_map(str(suffix))
    if isinstance(key_map, str):
        return resolve_object(key_map)
    return key_map


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
