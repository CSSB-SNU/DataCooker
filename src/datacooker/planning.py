"""Size-aware execution planner.

Turns per-item sizes (from ``lmdb.index``) + a node spec into a memory-safe,
throughput-optimal PLAN: items are partitioned into size tiers, and each tier is
assigned the worker count / memory / shard count it can safely afford. Giants run
in their own low-concurrency lane while the small majority runs at full cores --
so a heavy per-item build (CIFMol reconstruction, ~E x value_bytes live) never
exhausts node RAM the way a memory-blind hash-shard does.

Memory model (one item in flight per worker):

    peak_item_live  ~= expansion * value_bytes           # E, measured per recipe
    safe n_jobs      = clamp(floor(budget / (E * tier_max_bytes)), 1, cores)
    budget           = node_mem * safety - base

Under the Ray fan-out these ``n_jobs`` values are advisory (Ray sizes concurrency
to measured node memory); the planner is kept for the parts Ray does not do --
the cross-node shard count and the per-tier SLURM ``--mem`` request.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from datacooker.lmdb.index import ItemSize, human

GB = 1024**3


@dataclass(frozen=True, slots=True)
class NodeSpec:
    mem_gb: float = 490.0
    cores: int = 112
    n_nodes: int = 8
    mem_safety: float = 0.85
    base_gb: float = 130.0            # source mmap page cache + parent + slack
    giant_njobs_ceiling: int = 8      # items forcing n_jobs below this are 'giant'

    @property
    def budget_bytes(self) -> float:
        return (self.mem_gb * self.mem_safety - self.base_gb) * GB


@dataclass(frozen=True, slots=True)
class Tier:
    name: str
    keys: tuple[str, ...]
    n_jobs: int
    mem_gb: int
    n_shards: int
    max_bytes: int

    @property
    def total_bytes(self) -> int:
        return 0  # filled by the planner via replace-free construction below


@dataclass(frozen=True, slots=True)
class Plan:
    tiers: tuple[Tier, ...]
    node: NodeSpec
    expansion: float

    def nonempty(self) -> list[Tier]:
        return [t for t in self.tiers if t.keys]

    def describe(self) -> str:
        lines = [
            f"budget/worker-lane={human(self.node.budget_bytes)}  "
            f"E={self.expansion} cores={self.node.cores} nodes={self.node.n_nodes}",
        ]
        for t in self.tiers:
            if not t.keys:
                lines.append(f"[{t.name}] 0 items -- skip")
                continue
            lines.append(
                f"[{t.name}] items={len(t.keys)} max={human(t.max_bytes)} "
                f"n_jobs={t.n_jobs} mem={t.mem_gb}g n_shards={t.n_shards}"
            )
        return "\n".join(lines)


def safe_njobs(max_bytes: int, budget_bytes: float, expansion: float, cores: int) -> int:
    if max_bytes <= 0:
        return cores
    n = math.floor(budget_bytes / (expansion * max_bytes))
    return max(1, min(cores, n))


def _snake_shards(keys: list[str], n_shards: int) -> list[list[str]]:
    """Distribute size-sorted keys across shards in a boustrophedon (snake)
    order so each shard gets a balanced mix of large and small -- not the
    largest all in one shard."""
    buckets: list[list[str]] = [[] for _ in range(n_shards)]
    for i, k in enumerate(keys):
        row = i // n_shards
        idx = i % n_shards if row % 2 == 0 else n_shards - 1 - (i % n_shards)
        buckets[idx].append(k)
    return buckets


def plan(
    sizes: list[ItemSize],
    node: NodeSpec = NodeSpec(),
    expansion: float = 100.0,
) -> Plan:
    """Partition ``sizes`` into small/medium/giant tiers with per-tier resources.

    * small  -- fits at full ``cores`` -> n_jobs=cores.
    * giant  -- forces n_jobs below ``giant_njobs_ceiling`` -> its own low-mem lane.
    * medium -- everything between.

    Tiers still set each shard's ``--mem`` and the cross-node shard count; the
    ``n_jobs`` they derive is advisory under the Ray fan-out.
    """
    budget = node.budget_bytes
    if budget <= 0:
        msg = "memory budget <= 0: check node.base_gb / node.mem_gb"
        raise ValueError(msg)
    E = expansion
    cores = node.cores

    small_max = budget / (E * cores)                        # <= -> full cores
    giant_min = budget / (E * node.giant_njobs_ceiling)     # >  -> giant lane

    small: list[ItemSize] = []
    medium: list[ItemSize] = []
    giant: list[ItemSize] = []
    for it in sizes:
        if it.value_bytes <= small_max:
            small.append(it)
        elif it.value_bytes >= giant_min:
            giant.append(it)
        else:
            medium.append(it)

    def _tier(name: str, items: list[ItemSize]) -> Tier:
        if not items:
            return Tier(name, (), 0, 0, 0, 0)
        items = sorted(items, key=lambda s: s.value_bytes, reverse=True)
        tmax = items[0].value_bytes
        njobs = safe_njobs(tmax, budget, E, cores)
        mem_gb = min(int(node.mem_gb), math.ceil(node.base_gb + njobs * E * tmax / GB))
        n_shards = min(node.n_nodes, len(items))
        keys = tuple(k for shard in _snake_shards([i.key for i in items], n_shards) for k in shard)
        # keys kept flat here; the driver re-shards via _snake_shards when it
        # materializes per-shard key-lists. Store max/njobs so the driver can.
        return Tier(name, keys, njobs, mem_gb, n_shards, tmax)

    tiers = (
        _tier("small", small),
        _tier("medium", medium),
        _tier("giant", giant),
    )
    return Plan(tiers=tiers, node=node, expansion=E)


def shard_keys(tier: Tier) -> list[list[str]]:
    """Re-derive the per-shard key-lists for a tier (snake-balanced)."""
    if not tier.keys or tier.n_shards <= 1:
        return [list(tier.keys)]
    return _snake_shards(list(tier.keys), tier.n_shards)
