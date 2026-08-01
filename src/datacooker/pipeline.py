"""Planning-first pipeline driver.

The single entry point that makes every LMDB build/rebuild *plan from per-item
sizes before it submits*. Given a declarative config it:

  1. sizes the items -- stored value bytes (rebuild, from ``<db>.index.tsv``) or
     source-file bytes (build) -- building the source index if missing;
  2. plans size tiers (small @ full cores; giants isolated) from those sizes;
  3. submits one job array per tier via an Executor, each shard running the
     engine (``datacooker.cli.lmdb {build,rebuild}``) with the tier's item list
     (the engine fans out over Ray -- no per-shard worker knobs);
  4. merges the tier-shard outputs into the target (``afterok``);
  5. writes the target's own ``.index.tsv`` + ``.meta.json`` (mandatory post-step).

Hash-sharding is gone: items are partitioned by size, so giants never cluster
onto a worker and OOM. No human ever tunes n_jobs / mem / chunk / shards.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from datacooker.executors.base import Executor, TaskHandle
from datacooker.lmdb.index import ItemSize
from datacooker.planning import NodeSpec, plan, shard_keys

CLI = ["python", "-u", "-m", "datacooker.cli.lmdb"]

# per-op flag that hands the engine its item subset for a tier-shard
LIST_FLAG = {"rebuild": "--key-list", "build": "--file-list"}


@dataclass
class PipelineConfig:
    op: str                                  # "rebuild" | "build"
    target: Path                             # output LMDB
    engine_config: Path                      # datacooker yaml: recipe/adapter/params/...
    source: Path | None = None               # lineage note (source LMDB / data dir)
    expansion: float = 100.0                 # E, per-recipe (measured)
    node: NodeSpec = field(default_factory=NodeSpec)
    schema: str = "unknown"
    map_size: int = 300_000_000_000
    merge_on_duplicate: str = "overwrite"

    def __post_init__(self) -> None:
        self.target = Path(self.target)
        self.engine_config = Path(self.engine_config)
        if self.source is not None:
            self.source = Path(self.source)
        if self.op not in LIST_FLAG:
            msg = f"unsupported op {self.op!r}; expected one of {sorted(LIST_FLAG)}"
            raise ValueError(msg)

    @property
    def base(self) -> str:
        """Target path without ``.lmdb`` -- tier-shard outputs hang off this."""
        s = str(self.target)
        return s[:-5] if s.endswith(".lmdb") else s


def run(
    cfg: PipelineConfig,
    executor: Executor,
    keyfiles_dir: Path,
    sizes: list[ItemSize],
) -> TaskHandle:
    """Plan ``sizes`` into tiers and submit tier arrays -> merge -> index, chained.

    ``sizes`` are (item_id, bytes): item_id is an LMDB key (rebuild) or a source
    file path (build); ``bytes`` is stored value size or file size. The caller
    supplies them so this driver stays op-agnostic.
    """
    keyfiles_dir = Path(keyfiles_dir)
    keyfiles_dir.mkdir(parents=True, exist_ok=True)
    list_flag = LIST_FLAG[cfg.op]

    print(f"[pipeline] {cfg.op}: planning {len(sizes)} items -> {cfg.target.name}")
    p = plan(sizes, cfg.node, cfg.expansion)
    print(p.describe())

    # one array per tier
    tier_handles: list[TaskHandle] = []
    for tier in p.nonempty():
        shards = shard_keys(tier)
        for i, ids in enumerate(shards):
            (keyfiles_dir / f"items_{tier.name}_shard{i}.txt").write_text(
                "\n".join(ids) + "\n", encoding="utf-8")

        def argv_for_shard(i: int, tier=tier) -> list[str]:
            lst = keyfiles_dir / f"items_{tier.name}_shard{i}.txt"
            out = f"{cfg.base}_{tier.name}_shard{i}.lmdb"
            argv = [*CLI, cfg.op, str(cfg.engine_config),
                    list_flag, str(lst), "--output", out,
                    "--n-jobs", str(tier.n_jobs)]
            if cfg.op == "rebuild":
                argv += ["--n-shards", "1"]
            else:  # build: resume-safe re-runs (rebuild already skip_existing=True)
                argv += ["--skip-existing"]
            # Bound the drain buffer to ~a few worker-waves: --chunk-size caps how
            # many results are regrouped before a write, keeping the parent at
            # worker-memory scale (the 10k default let it hoard and OOM).
            argv += ["--chunk-size", str(max(64, 4 * tier.n_jobs))]
            return argv

        tier_handles.append(executor.run_tier_shards(
            name=f"{cfg.target.stem}_{tier.name}",
            argv_for_shard=argv_for_shard,
            n_shards=len(shards),
            mem_gb=tier.mem_gb,
            cores=cfg.node.cores,
        ))

    # merge (afterok all tiers)
    merge = executor.run_once(
        name=f"{cfg.target.stem}_merge",
        argv=[*CLI, "merge", f"{cfg.base}_*_shard*.lmdb",
              "--output", str(cfg.target), "--map-size", str(cfg.map_size),
              "--overwrite", "--on-duplicate", cfg.merge_on_duplicate],
        mem_gb=200, cores=16,
        depends_on=tuple(h.job_id for h in tier_handles if h.job_id),
    )

    # index + meta for the NEW db (afterok merge) -- mandatory post-step
    return executor.run_once(
        name=f"{cfg.target.stem}_index",
        argv=[*CLI, "index", str(cfg.target),
              "--schema", cfg.schema,
              *(["--source", str(cfg.source)] if cfg.source else []),
              "--expansion", str(cfg.expansion)],
        mem_gb=64, cores=8,
        depends_on=(merge.job_id,) if merge.job_id else (),
    )
