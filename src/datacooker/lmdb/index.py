"""Per-item size index + per-DB metadata sidecars.

Every LMDB a workflow produces gets two sibling files:

* ``<db>.index.tsv``  -- one row per key with its stored value size in bytes,
  largest first, plus a ``#``-comment header of summary stats. This is the
  per-item information the planner reads to size a build BEFORE submitting it.
  Sizes are read with a zero-copy ``buffers=True`` cursor, so ``len(value)``
  comes from the mmap page header without materializing the payload -- safe and
  fast even on multi-hundred-GB DBs.

* ``<db>.meta.json`` -- per-DB descriptor: schema class, key convention, codec
  framing, entry count, source lineage, and the measured memory-expansion factor
  the planner needs. Written as a mandatory post-step of every build so coverage
  is verifiable and the next build can plan against it.

Human-readable by design (TSV + JSON) and dependency-light (only ``lmdb``).
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from types import ModuleType

INDEX_SUFFIX = ".index.tsv"
META_SUFFIX = ".meta.json"


def index_path(db_path: str | Path) -> Path:
    """Return the ``<db>.index.tsv`` sidecar path for a DB."""
    return Path(str(db_path) + INDEX_SUFFIX)


def meta_path(db_path: str | Path) -> Path:
    """Return the ``<db>.meta.json`` sidecar path for a DB."""
    return Path(str(db_path) + META_SUFFIX)


def human(n: float) -> str:
    """Format a byte count as a human-readable IEC size string."""
    for unit in ("B", "K", "M", "G", "T"):
        if abs(n) < 1024.0:
            return f"{n:.1f}{unit}"
        n /= 1024.0
    return f"{n:.1f}P"


@dataclass(frozen=True, slots=True)
class ItemSize:
    """One key and its stored value size in bytes."""

    key: str
    value_bytes: int


@dataclass(frozen=True, slots=True)
class IndexStats:
    """Summary statistics over a DB's per-item value sizes."""

    entries: int
    total_bytes: int
    max_bytes: int
    max_key: str
    median_bytes: int
    n_gt_100mb: int
    n_gt_1gb: int


@dataclass
class DbMeta:
    """Per-DB descriptor written to ``<db>.meta.json``."""

    schema: str = "unknown"          # A|B|C|D|E|F (see the schema catalogue)
    key_convention: str = "unknown"  # pdbid | pdbid_chain | seqid | mgyp | ...
    codec: str = "biomol_zstd"       # biomol_zstd | biomol_zstd_streaming | raw
    entries: int = 0
    total_bytes: int = 0
    max_bytes: int = 0
    expansion: float | None = None   # measured peak-live-mem / value_bytes (E)
    source: str | None = None        # lineage: the DB/dir this was built from
    recipe: str | None = None
    extra: dict = field(default_factory=dict)


def _require_lmdb() -> ModuleType:
    try:
        import lmdb
    except ImportError as exc:  # pragma: no cover
        msg = "python-lmdb is required for size indexing."
        raise RuntimeError(msg) from exc
    return lmdb


def scan_sizes(db_path: str | Path) -> list[ItemSize]:
    """Return every key's stored value size via a zero-copy cursor."""
    lmdb = _require_lmdb()
    env = lmdb.open(str(db_path), readonly=True, lock=False, subdir=True, max_readers=4096)
    try:
        out: list[ItemSize] = []
        with env.begin(buffers=True) as txn:
            for key, value in txn.cursor().iternext(keys=True, values=True):
                out.append(ItemSize(bytes(key).decode("utf-8", "replace"), len(value)))
        return out
    finally:
        env.close()


def _stats(sizes: list[ItemSize]) -> IndexStats:
    if not sizes:
        return IndexStats(0, 0, 0, "", 0, 0, 0)
    ordered = sorted(sizes, key=lambda s: s.value_bytes, reverse=True)
    vals = [s.value_bytes for s in ordered]
    return IndexStats(
        entries=len(ordered),
        total_bytes=sum(vals),
        max_bytes=vals[0],
        max_key=ordered[0].key,
        median_bytes=vals[len(vals) // 2],
        n_gt_100mb=sum(1 for v in vals if v > 100 * 1024**2),
        n_gt_1gb=sum(1 for v in vals if v > 1024**3),
    )


def write_index(db_path: str | Path, sizes: list[ItemSize] | None = None) -> IndexStats:
    """Write ``<db>.index.tsv`` (largest first, with stats header).

    Scans the DB if ``sizes`` is not supplied. Returns the summary stats.
    """
    if sizes is None:
        sizes = scan_sizes(db_path)
    ordered = sorted(sizes, key=lambda s: s.value_bytes, reverse=True)
    stats = _stats(ordered)
    out = index_path(db_path)
    with out.open("w", encoding="utf-8") as f:
        f.write(f"# db\t{Path(db_path).name}\n")
        f.write(f"# entries\t{stats.entries}\n")
        f.write(f"# total_bytes\t{stats.total_bytes}\t{human(stats.total_bytes)}\n")
        f.write(f"# max_bytes\t{stats.max_bytes}\t{human(stats.max_bytes)}\t{stats.max_key}\n")
        f.write(f"# median_bytes\t{stats.median_bytes}\t{human(stats.median_bytes)}\n")
        f.write(f"# n_value_gt_100MB\t{stats.n_gt_100mb}\n")
        f.write(f"# n_value_gt_1GB\t{stats.n_gt_1gb}\n")
        f.write("key\tvalue_bytes\tvalue_human\n")
        for s in ordered:
            f.write(f"{s.key}\t{s.value_bytes}\t{human(s.value_bytes)}\n")
    return stats


def read_index(db_path: str | Path) -> list[ItemSize]:
    """Read per-key sizes from the sidecar (no DB open). Empty if absent."""
    p = index_path(db_path)
    if not p.exists():
        return []
    out: list[ItemSize] = []
    with p.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith(("#", "key\t")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    out.append(ItemSize(parts[0], int(parts[1])))
                except ValueError:
                    continue
    return out


def has_index(db_path: str | Path) -> bool:
    """Return whether the ``<db>.index.tsv`` sidecar exists."""
    return index_path(db_path).exists()


def write_meta(db_path: str | Path, meta: DbMeta) -> Path:
    """Write ``<db>.meta.json``."""
    p = meta_path(db_path)
    p.write_text(json.dumps(asdict(meta), indent=2), encoding="utf-8")
    return p


def read_meta(db_path: str | Path) -> DbMeta | None:
    """Read ``<db>.meta.json`` into a :class:`DbMeta`, or ``None`` if absent."""
    p = meta_path(db_path)
    if not p.exists():
        return None
    return DbMeta(**json.loads(p.read_text(encoding="utf-8")))
