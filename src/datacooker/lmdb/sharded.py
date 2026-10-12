"""Immutable shard collections with the read API used by DataCooker.

The manifest is the logical database. Its referenced shards must remain immutable
and must not be deleted as scratch files. No payload copy is made at publication.
"""
from __future__ import annotations

import heapq
import json
import os
import tempfile
from contextlib import ExitStack
from pathlib import Path
from typing import Any

import lmdb
from typing_extensions import Self

MERGE_LIMIT = 5_000_000
MANIFEST = "shards.json"


def should_merge(entries: int) -> bool:
    """Exactly five million entries still use a merged database."""
    return entries <= MERGE_LIMIT


def shard_paths(path: str | Path) -> list[Path]:
    """Resolve and validate a published collection's physical shard paths."""
    root = Path(path)
    data = json.loads((root / MANIFEST).read_text())
    if data.get("format") != "shard_lmdb" or data.get("version") != 1:
        msg = f"Unsupported shard manifest: {root}"
        raise ValueError(msg)
    paths = [(root / item["path"]).resolve() for item in data["shards"]]
    if not paths or len(set(paths)) != len(paths):
        msg = "Empty or repeated shard paths"
        raise ValueError(msg)
    if any(not (p / "data.mdb").is_file() for p in paths):
        msg = f"Missing physical shard in {root}"
        raise FileNotFoundError(msg)
    return paths


class ShardCursor:
    """Globally sorted, forward-only cursor over disjoint physical shards."""

    def __init__(self, txns: list[Any]) -> None:
        self.txns = txns

    def __iter__(self) -> Any:
        """Iterate key/value pairs."""
        return self.iternext()

    def iternext(self, *, keys: bool = True, values: bool = True) -> Any:
        """Iterate without materializing values when only keys are requested."""
        streams = [t.cursor().iternext(keys=True, values=values) for t in self.txns]
        previous = None
        for item in heapq.merge(*streams, key=(lambda kv: bytes(kv[0])) if values else bytes):
            key = bytes(item[0]) if values else bytes(item)
            if key == previous:
                msg = f"Duplicate key across shards: {key!r}"
                raise ValueError(msg)
            previous = key
            yield item if keys else item[1]


class ShardTransaction:
    """Read-only snapshot transactions; get probes the component shards."""

    def __init__(self, envs: list[Any], *, buffers: bool = False) -> None:
        self.stack = ExitStack()
        try:
            self.txns = [self.stack.enter_context(e.begin(buffers=buffers)) for e in envs]
        except BaseException:
            self.stack.close()
            raise

    def __enter__(self) -> Self:
        """Enter the read transaction."""
        return self

    def __exit__(self, *_: object) -> None:
        """Release readers."""
        self.abort()

    def abort(self) -> None:
        """Release all reader transactions."""
        self.stack.close()

    def cursor(self) -> ShardCursor:
        """Return a forward cursor in global key order."""
        return ShardCursor(self.txns)

    def get(self, key: bytes, default: Any = None) -> Any:
        """Return a key's payload, or the supplied default if absent."""
        for txn in self.txns:
            value = txn.get(key)
            if value is not None:
                return value
        return default

    def stat(self) -> dict[str, int]:
        """Return the logical entry count."""
        return {"entries": sum(t.stat()["entries"] for t in self.txns)}


class ShardEnvironment:
    """Read-only environment for the published shard manifest."""

    def __init__(self, path: str | Path, **kwargs: Any) -> None:
        self.stack = ExitStack()
        try:
            self.envs = [self.stack.enter_context(lmdb.open(str(p), **kwargs)) for p in shard_paths(path)]
        except BaseException:
            self.stack.close()
            raise

    def __enter__(self) -> Self:
        """Enter the read environment."""
        return self

    def __exit__(self, *_: object) -> None:
        """Release readers."""
        self.close()

    def close(self) -> None:
        """Release physical environments."""
        self.stack.close()

    def begin(self, *, buffers: bool = False, write: bool = False) -> ShardTransaction:
        """Start a read-only logical transaction."""
        if write:
            msg = "Published shard_lmdb collections are immutable"
            raise ValueError(msg)
        return ShardTransaction(self.envs, buffers=buffers)


def open_env(path: str | Path, **kwargs: Any) -> Any:
    """Open an ordinary LMDB or an immutable shard collection."""
    if (Path(path) / MANIFEST).is_file():
        if not kwargs.get("readonly", False):
            msg = "Published shard_lmdb collections are immutable"
            raise ValueError(msg)
        return ShardEnvironment(path, **kwargs)
    return lmdb.open(str(path), **kwargs)


def publish_shards(paths: list[Path], target: Path, *, overwrite: bool = False) -> int:
    """Check unique keys and atomically publish references, never copy payloads."""
    paths = [p.resolve() for p in paths]
    target = target.resolve()
    if not paths or len(set(paths)) != len(paths) or target in paths:
        msg = "Require unique physical shards distinct from the output"
        raise ValueError(msg)
    if target.exists() and ((target / "data.mdb").exists() or not overwrite):
        msg = f"Refusing to replace existing database: {target}"
        raise FileExistsError(msg)
    with ExitStack() as stack:
        envs = [stack.enter_context(lmdb.open(str(p), readonly=True, lock=False, create=False)) for p in paths]
        txns = [stack.enter_context(e.begin()) for e in envs]
        counts = [t.stat()["entries"] for t in txns]
        # Only B-tree keys are visited: no decompression or copying of values.
        checked = sum(1 for _ in ShardCursor(txns).iternext(values=False))
        if checked != sum(counts):
            msg = "Shard count changed during publication"
            raise ValueError(msg)
    _write_manifest(target, paths, counts)
    return checked


def _write_manifest(target: Path, paths: list[Path], counts: list[int]) -> None:
    target.mkdir(parents=True, exist_ok=True)
    descriptor = {"format": "shard_lmdb", "version": 1, "entries": sum(counts),
                  "merge_limit": MERGE_LIMIT,
                  "shards": [{"path": os.path.relpath(p, target), "entries": n}
                             for p, n in zip(paths, counts, strict=True)]}
    with tempfile.NamedTemporaryFile(mode="w", dir=target, delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(descriptor, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(target / MANIFEST)


def patch_lmdb(target: str | Path, patches: str | Path | list[str | Path]) -> dict[str, int]:
    """Replace or add the records of ``patches`` (one LMDB or several) in ``target``.

    For records that were built wrong, or not at all, when rebuilding the whole DB is not
    worth it. A shard collection drops the patched keys from its shards and gains the
    patches as more immutable shards; a single LMDB takes the patch values in place.
    """
    target = Path(target).resolve()
    plist = [Path(p).resolve() for p in (patches if isinstance(patches, list) else [patches])]
    with ExitStack() as stack:
        ptxns = [stack.enter_context(stack.enter_context(
            lmdb.open(str(p), readonly=True, lock=False)).begin()) for p in plist]
        keys = [bytes(k) for t in ptxns for k in t.cursor().iternext(values=False)]
        if len(set(keys)) != len(keys):
            msg = "Patch LMDBs share keys"
            raise ValueError(msg)
        if (target / MANIFEST).is_file():
            shards = [p.resolve() for p in shard_paths(target)]
            if set(plist) & set(shards):
                msg = "A patch is already one of the target's shards"
                raise ValueError(msg)
            removed = 0
            for path in shards:
                with lmdb.open(str(path), lock=True, map_size=1 << 40) as env, env.begin(write=True) as txn:
                    removed += sum(bool(txn.delete(k)) for k in keys)
            # The published shards already hold unique keys and none of the patch keys any
            # more, so republishing needs counts and a look-up of the patch keys -- not a
            # walk of every key (16M keys, an hour on GPFS, to patch 14 records).
            counts = []
            for path in shards:
                with lmdb.open(str(path), readonly=True, lock=False) as env, env.begin() as txn:
                    if any(txn.get(k) is not None for k in keys):
                        msg = f"Patched key still present in {path}"
                        raise ValueError(msg)
                    counts.append(txn.stat()["entries"])
            counts += [t.stat()["entries"] for t in ptxns]
            _write_manifest(target, [*shards, *plist], counts)
            return {"patched": len(keys), "replaced": removed, "added": len(keys) - removed, "entries": sum(counts)}
        size = (target / "data.mdb").stat().st_size + sum((p / "data.mdb").stat().st_size for p in plist)
        with lmdb.open(str(target), map_size=max(1 << 30, int(size * 1.25))) as env, env.begin(write=True) as txn:
            replaced = 0
            for ptxn in ptxns:
                for k, v in ptxn.cursor():
                    replaced += txn.get(k) is not None
                    txn.put(k, v)
            entries = txn.stat()["entries"]
    return {"patched": len(keys), "replaced": replaced, "added": len(keys) - replaced, "entries": entries}


def update_index_after_patch(target: str | Path, patches: list[str | Path]) -> bool:
    """Rewrite ``target``'s index from its old rows plus the patches' sizes.

    Re-indexing a multi-TB DB reads every value (hours, and it ran a node out of memory);
    only the patched keys changed. Returns False when there is no old index to update.
    """
    from datacooker.lmdb.index import (
        has_index,
        read_index,
        scan_sizes,
        write_index,
    )

    if not has_index(target):
        return False
    new = {s.key: s for p in patches for s in scan_sizes(p)}
    old = read_index(target)
    seen = {s.key for s in old}
    write_index(target, [new.get(s.key, s) for s in old] + [s for k, s in new.items() if k not in seen])
    return True

