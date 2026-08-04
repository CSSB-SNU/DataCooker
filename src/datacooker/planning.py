"""Execution planner for the Ray fan-out.

Within-node memory is **Ray's** job now (measured admission control -- it throttles
and retries tasks by real node memory), so the planner no longer predicts per-item
memory (the old ``peak = n_jobs x E x max_bytes`` model that OOM'd on its blind
spots). It does only what Ray does not: decide **how many nodes** to spread a build
across, and **balance** heavy items across them so no shard is lopsided.

Every shard is given the **whole node** (``mem_gb`` / ``cores``); Ray self-regulates
inside that. ``expansion`` (E) survives only as a recorded stat, not a knob.
"""
from __future__ import annotations

from dataclasses import dataclass

from datacooker.lmdb.index import ItemSize, human


@dataclass(frozen=True, slots=True)
class NodeSpec:
    """Whole-node resource profile a single shard is pinned to."""

    mem_gb: float = 490.0
    cores: int = 112
    n_nodes: int = 8


@dataclass(frozen=True, slots=True)
class Tier:
    """One size lane of items to be split across snake-balanced shards."""

    name: str
    keys: tuple[str, ...]
    n_jobs: int          # advisory under Ray (it sizes to node CPUs)
    mem_gb: int          # SLURM --mem per shard (whole node)
    n_shards: int
    max_bytes: int


@dataclass(frozen=True, slots=True)
class Plan:
    """Full fan-out plan: the tiers to submit plus the node profile and stats."""

    tiers: tuple[Tier, ...]
    node: NodeSpec
    expansion: float

    def nonempty(self) -> list[Tier]:
        """Return the tiers that actually have items to build."""
        return [t for t in self.tiers if t.keys]

    def describe(self) -> str:
        """Render a one-line-per-tier human summary of the plan."""
        lines = [
            f"ray fan-out: cores={self.node.cores} nodes<={self.node.n_nodes} "
            f"(within-node memory = Ray admission control)",
        ]
        for t in self.tiers:
            if not t.keys:
                lines.append(f"[{t.name}] 0 items -- skip")
                continue
            lines.append(
                f"[{t.name}] items={len(t.keys)} max={human(t.max_bytes)} "
                f"mem={t.mem_gb}g n_shards={t.n_shards}"
            )
        return "\n".join(lines)


def _snake_shards(keys: list[str], n_shards: int) -> list[list[str]]:
    """Distribute size-sorted keys across shards in boustrophedon (snake) order.

    Each shard gets a balanced mix of large and small keys -- not the largest all
    in one shard.
    """
    buckets: list[list[str]] = [[] for _ in range(n_shards)]
    for i, k in enumerate(keys):
        row = i // n_shards
        idx = i % n_shards if row % 2 == 0 else n_shards - 1 - (i % n_shards)
        buckets[idx].append(k)
    return buckets


_DEFAULT_NODE = NodeSpec()


def plan(
    sizes: list[ItemSize],
    node: NodeSpec = _DEFAULT_NODE,
    expansion: float = 100.0,
) -> Plan:
    """Split items into <= ``n_nodes`` snake-balanced shards, one full-node lane.

    One tier ("all"): each shard is a SLURM job pinned to a whole node running Ray;
    the planner sets the node count and balances heavy items across shards, and Ray
    owns memory within each node. ``expansion`` is recorded, not used to size.
    """
    items = sorted(sizes, key=lambda s: s.value_bytes, reverse=True)
    if not items:
        return Plan((Tier("all", (), node.cores, int(node.mem_gb), 0, 0),), node, expansion)
    tmax = items[0].value_bytes
    n_shards = min(node.n_nodes, len(items))
    keys = tuple(
        k for shard in _snake_shards([i.key for i in items], n_shards) for k in shard
    )
    tier = Tier("all", keys, node.cores, int(node.mem_gb), n_shards, tmax)
    return Plan((tier,), node, expansion)


def shard_keys(tier: Tier) -> list[list[str]]:
    """Re-derive the per-shard key-lists for a tier (snake-balanced)."""
    if not tier.keys or tier.n_shards <= 1:
        return [list(tier.keys)]
    return _snake_shards(list(tier.keys), tier.n_shards)
