"""Execution planner for the Ray fan-out.

The planner does exactly one thing: decide how many whole-node shards to spread a
build across, and snake-balance the items so no shard is lopsided. It does NOT predict
per-item memory or split items into size tiers -- that predict-with-E model had blind
spots (a record that is cheap on disk but heavy once materialised) and needed per-DB
tuning. Within-node concurrency is instead bounded by the engine's *measured* memory
admission gate (``datacooker._ray``): a light record runs at full width, a heavy one
simply slows intake, so one lane is enough and there is nothing to size.
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
    """One lane of items to be split across snake-balanced whole-node shards."""

    name: str
    keys: tuple[str, ...]
    n_jobs: int          # advisory under Ray (it sizes to node CPUs)
    mem_gb: int          # SLURM --mem per shard (whole node)
    n_shards: int
    max_bytes: int


@dataclass(frozen=True, slots=True)
class Plan:
    """Full fan-out plan: the lane(s) to submit plus the node profile and stats."""

    tiers: tuple[Tier, ...]
    node: NodeSpec
    expansion: float

    def nonempty(self) -> list[Tier]:
        """Return the tiers that actually have items to build."""
        return [t for t in self.tiers if t.keys]

    def describe(self) -> str:
        """Render a one-line-per-lane human summary of the plan."""
        lines = [
            f"ray fan-out: cores={self.node.cores} nodes<={self.node.n_nodes} "
            f"(within-node memory = Ray measured admission)",
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
    """Snake-balance items across ``<= n_nodes`` whole-node shards, one lane.

    Each shard is a SLURM job pinned to a whole node running Ray; the planner only sets
    the node count and spreads heavy items evenly, and Ray's measured admission gate owns
    memory within each node. ``expansion`` is recorded as a stat, not used to size.
    """
    items = sorted(sizes, key=lambda s: s.value_bytes, reverse=True)
    if not items:
        return Plan((Tier("all", (), node.cores, int(node.mem_gb), 0, 0),), node, expansion)
    tmax = items[0].value_bytes
    n_shards = min(node.n_nodes, len(items))
    keys = tuple(i.key for i in items)
    tier = Tier("all", keys, node.cores, int(node.mem_gb), n_shards, tmax)
    return Plan((tier,), node, expansion)


def shard_keys(tier: Tier) -> list[list[str]]:
    """Re-derive the per-shard key-lists for a lane (snake-balanced)."""
    if not tier.keys or tier.n_shards <= 1:
        return [list(tier.keys)]
    return _snake_shards(list(tier.keys), tier.n_shards)
