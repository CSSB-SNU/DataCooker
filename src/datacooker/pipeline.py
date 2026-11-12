"""Size-balanced pipeline submission with exact shard manifests.

The planner distributes descending-size inputs across nodes once. Each shard runs
bounded Ray workers, then an afterok merge publishes a fresh database and indexing
writes its sidecars. Input size balances I/O; it does not prove peak memory fits.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from uuid import uuid4

from datacooker.executors.base import Executor, TaskHandle
from datacooker.lmdb.index import ItemSize
from datacooker.planning import NodeSpec, Tier, plan, shard_keys

CLI = ["python", "-u", "-m", "datacooker.cli.lmdb"]

# per-op flag that hands the engine its item subset for a tier-shard
LIST_FLAG = {"rebuild": "--key-list", "build": "--file-list"}


@dataclass
class PipelineConfig:
    """Declarative inputs for one planning-first LMDB build/rebuild pipeline."""

    op: str                                  # "rebuild" | "build"
    target: Path                             # output LMDB
    engine_config: Path                      # datacooker yaml: recipe/adapter/params/...
    source: Path | None = None               # lineage note (source LMDB / data dir)
    expansion: float = 100.0                 # E, per-recipe (measured)
    node: NodeSpec = field(default_factory=NodeSpec)
    schema: str = "unknown"
    codec: str = "biomol_zstd"               # value codec recorded in meta.json (the index step)
    map_size: int = 300_000_000_000
    merge_on_duplicate: str = "error"
    run_id: str = field(default_factory=lambda: uuid4().hex, init=False)

    def __post_init__(self) -> None:
        """Normalise paths and validate ``op``."""
        self.target = Path(self.target)
        self.engine_config = Path(self.engine_config)
        if self.source is not None:
            self.source = Path(self.source)
        if self.op not in LIST_FLAG:
            msg = f"unsupported op {self.op!r}; expected one of {sorted(LIST_FLAG)}"
            raise ValueError(msg)

    @property
    def base(self) -> str:
        """Unique per-run shard prefix, on the destination filesystem."""
        return str(self.target.parent / f".{self.target.stem}.build" / self.run_id / self.target.stem)


def run(
    cfg: PipelineConfig,
    executor: Executor,
    keyfiles_dir: Path,
    sizes: list[ItemSize],
    depends_on: tuple[str, ...] = (),
) -> TaskHandle:
    """Plan ``sizes`` into tiers and submit tier arrays -> merge -> index, chained.

    ``sizes`` are (item_id, bytes): item_id is an LMDB key (rebuild) or a source
    file path (build); ``bytes`` is stored value size or file size. The caller
    supplies them so this driver stays op-agnostic.

    ``depends_on`` are upstream job ids the FIRST stage (the tier arrays) must wait
    for -- this is how a multi-DB build chains one config's whole pipeline after
    another's (e.g. cif_attached's tiers afterok cif's index).
    """
    keyfiles_dir = Path(keyfiles_dir)
    keyfiles_dir.mkdir(parents=True, exist_ok=True)
    list_flag = LIST_FLAG[cfg.op]

    print(f"[pipeline] {cfg.op}: planning {len(sizes)} items -> {cfg.target.name}")  # noqa: T201
    p = plan(sizes, cfg.node, cfg.expansion)
    print(p.describe())  # noqa: T201

    # one array per tier
    tier_handles: list[TaskHandle] = []
    shard_outputs: list[str] = []
    for tier in p.nonempty():
        shards = shard_keys(tier)
        shard_outputs.extend(f"{cfg.base}_{tier.name}_shard{i}.lmdb"
                             for i in range(len(shards)))
        for i, ids in enumerate(shards):
            (keyfiles_dir / f"items_{tier.name}_shard{i}.txt").write_text(
                "\n".join(ids) + "\n", encoding="utf-8")

        def argv_for_shard(i: int, tier: Tier = tier) -> list[str]:
            lst = keyfiles_dir / f"items_{tier.name}_shard{i}.txt"
            out = f"{cfg.base}_{tier.name}_shard{i}.lmdb"
            argv = [*CLI, cfg.op, str(cfg.engine_config),
                    list_flag, str(lst), "--output", out,
                    "--n-jobs", str(tier.n_jobs)]
            if cfg.op == "rebuild":
                argv += ["--n-shards", "1"]
            else:  # build: resume-safe re-runs (rebuild already skip_existing=True)
                argv += ["--skip-existing"]
            # Bound the drain buffer to ~a few worker-waves: --chunk-size caps how many
            # results are regrouped before a write, so the parent holds at most that many
            # (small, serialized) records at once (the 10k default let it hoard and OOM).
            # One item per Ray task (the engine default) plus the measured-memory gate own
            # the within-node concurrency; there is nothing per-shard to size here.
            argv += ["--chunk-size", str(max(64, 4 * tier.n_jobs))]
            return argv

        tier_handles.append(executor.run_tier_shards(
            name=f"{cfg.target.stem}_{tier.name}",
            argv_for_shard=argv_for_shard,
            n_shards=len(shards),
            mem_gb=tier.mem_gb,
            cores=cfg.node.cores,
            depends_on=depends_on,
        ))

    shard_manifest = keyfiles_dir / "merge_shards.txt"
    shard_manifest.write_text("".join(path + "\n" for path in shard_outputs), encoding="utf-8")

    # Finalize after all tiers: <=5M merges, >5M publishes immutable shards.
    merge = executor.run_once(
        name=f"{cfg.target.stem}_publish",
        argv=[*CLI, "merge", "--shard-list", str(shard_manifest),
              "--output", str(cfg.target), "--map-size", str(cfg.map_size),
              "--overwrite", "--on-duplicate", cfg.merge_on_duplicate],
        mem_gb=min(200, int(cfg.node.mem_gb)), cores=min(16, cfg.node.cores),
        depends_on=tuple(h.job_id for h in tier_handles if h.job_id),
    )

    # index + meta for the NEW db (afterok merge) -- mandatory post-step
    return executor.run_once(
        name=f"{cfg.target.stem}_index",
        argv=[*CLI, "index", str(cfg.target),
              "--schema", cfg.schema, "--codec", cfg.codec,
              *(["--source", str(cfg.source)] if cfg.source else []),
              "--expansion", str(cfg.expansion)],
        mem_gb=min(64, int(cfg.node.mem_gb)), cores=min(8, cfg.node.cores),
        depends_on=(merge.job_id,) if merge.job_id else (),
    )
