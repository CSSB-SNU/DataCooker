"""Core LMDB workflow helpers for recipe-driven database pipelines."""

from __future__ import annotations

import gc
import heapq
import logging
import os
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import import_module
from pathlib import Path
from typing import Any, Literal

from datacooker._parallel import iter_parallel_chunks, run_serial_test
from datacooker.api import parse_file
from datacooker.lmdb._shared import load_metadata, rebuild_entry
from datacooker.protocols import (
    ConvertFunc,
    DeserializeFunc,
    KeyFunc,
    LoadFunc,
    SerializeFunc,
    TransformFunc,
)
from datacooker.readers import ReaderHooks, decode_payload, resolve_key_transform
from datacooker.utils.sharding import shard_items_by_key
from datacooker.writers import WriterHooks, encode_output

logger = logging.getLogger(__name__)

_ENV_CACHE: dict[tuple[str, bool, bool, int], Any] = {}
ErrorMode = Literal["skip", "raise"]
DuplicatePolicy = Literal["overwrite", "skip", "error", "union"]

# Maps an existing LMDB key to its new key, or ``None`` to drop the entry.
KeyMapFunc = Callable[[str], "str | None"]


@dataclass(frozen=True, slots=True)
class LmdbWriteReport:
    """Summarize a write-oriented LMDB workflow run."""

    attempted: int
    written: int
    skipped_existing: int = 0
    skipped_empty: int = 0
    failed: int = 0
    duplicate_keys: int = 0
    failed_keys: tuple[str, ...] = field(default_factory=tuple)


def default_lmdb_key(path: Path) -> str:
    """Derive the default LMDB key from a file path."""
    return path.name.split(".")[0]


def extract_lmdb_keys(env_path: Path) -> list[str]:
    """Return all keys in an LMDB database."""
    lmdb_module = _require_lmdb()
    if not env_path.exists():
        return []

    env = lmdb_module.open(str(env_path), readonly=True, lock=False)
    try:
        with env.begin() as txn:
            return [
                key.decode() for key in txn.cursor().iternext(keys=True, values=False)
            ]
    finally:
        env.close()


def count_lmdb_entries(env_path: Path) -> int:
    """Count entries in an LMDB database."""
    lmdb_module = _require_lmdb()
    if not env_path.exists():
        return 0

    env = lmdb_module.open(str(env_path), readonly=True, lock=False)
    try:
        with env.begin() as txn:
            return sum(1 for _ in txn.cursor())
    finally:
        env.close()


def filter_pending_lmdb_paths(
    paths: Sequence[Path],
    env_path: Path,
    *,
    key_func: Callable[[Path], str] = default_lmdb_key,
) -> list[Path]:
    """Filter out paths whose derived keys already exist in an LMDB database."""
    existing_keys = set(extract_lmdb_keys(env_path))
    return [path for path in paths if key_func(path) not in existing_keys]


def read_lmdb_raw(env_path: Path, key: str) -> bytes | None:
    """Read a raw LMDB payload without applying deserialization."""
    env = _get_env(env_path, readonly=True, lock=False)
    with env.begin(buffers=True) as txn:
        value = txn.get(key.encode())
    return None if value is None else bytes(value)


def read_all_lmdb_raw(env_path: Path) -> dict[str, bytes]:
    """Read every raw entry from an LMDB database."""
    env = _get_env(env_path, readonly=True, lock=False)
    data: dict[str, bytes] = {}
    with env.begin(buffers=True) as txn:
        for key, value in txn.cursor():
            data[bytes(key).decode()] = bytes(value)
    return data


def read_lmdb(
    env_path: Path,
    key: str,
    *,
    deserializer: DeserializeFunc,
) -> dict[str, Any]:
    """Read a decoded entry from an LMDB database."""
    value = read_lmdb_raw(env_path, key)
    if value is None:
        msg = f"Key '{key}' not found in LMDB database at '{env_path}'."
        raise KeyError(msg)
    return deserializer(value)


# Build worker state at module scope: with the "multiprocessing" backend the
# top-level worker is picklable by reference and forked workers inherit this
# dict (metadata map, recipe, reader) via copy-on-write instead of pickling it
# per task.
_BUILD_STATE: dict[str, Any] = {}


_GC_ENABLED_PIDS: set[int] = set()


def _ensure_worker_gc() -> None:
    """Re-enable cyclic GC once per forked worker: the parent disables it before
    the fork (to keep the inherited state copy-on-write), but a worker that keeps
    GC off never collects the circular references in its per-task objects, so its
    RssAnon climbs unbounded -> node OOM. Frozen parent objects are in the
    permanent generation (unscanned), so COW is preserved."""
    pid = os.getpid()
    if pid not in _GC_ENABLED_PIDS:
        gc.enable()
        _GC_ENABLED_PIDS.add(pid)


def _build_worker(data_file: Any) -> tuple[bytes, bytes, Exception | None]:
    state = _BUILD_STATE
    _ensure_worker_gc()
    key = state["key_func"](data_file)
    try:
        data_dict = dict(state["base_inputs"])
        data_dict.update(state["metadata_dict"])
        output = parse_file(
            recipe_path=state["recipe"],
            file_path=data_file,
            loader=state["loader"],
            inputs=data_dict,
            key_transform=state["key_transform"],
            **state["extra_kwargs"],
        )
        return key.encode(), encode_output(output, writer=state["writer"]), None
    except Exception as error:  # noqa: BLE001
        return key.encode(), b"", error


# Rebuild worker state at module scope, mirroring _BUILD_STATE. With the
# "multiprocessing" backend the top-level worker is picklable by reference and
# forked workers inherit the metadata map / recipe / reader via copy-on-write.
# The source LMDB handle is NOT inherited across fork (its fds/locks are unsafe
# to share); each worker lazily reopens its own read-only handle in _REBUILD_ENV.
_REBUILD_STATE: dict[str, Any] = {}
_REBUILD_ENV: dict[str, Any] = {}


def install_build_state(state: Mapping[str, Any]) -> None:
    """Populate ``_BUILD_STATE`` in a fresh worker process.

    The fork backend inherits ``_BUILD_STATE`` via copy-on-write; the ``ray``
    backend has no fork, so each Ray worker calls this once (with the state
    shipped from the driver's plasma object store) to make ``_build_worker``
    read the same recipe / metadata / hooks as the fork path.
    """
    _BUILD_STATE.clear()
    _BUILD_STATE.update(state)


def install_rebuild_state(state: Mapping[str, Any]) -> None:
    """Populate ``_REBUILD_STATE`` in a fresh worker process (ray backend).

    Mirror of :func:`install_build_state` for the rebuild path. Clears the
    per-process source-handle cache too, so a restarted Ray worker reopens its
    own read-only handle instead of reusing a dead one.
    """
    _REBUILD_STATE.clear()
    _REBUILD_ENV.clear()
    _REBUILD_STATE.update(state)

_RebuildResult = tuple[bytes, bytes, Exception | None]


def explode_source_key(key: str, sep: str) -> str:
    """Recover the source key an exploded destination key came from."""
    return key.split(sep, 1)[0]


def _validate_explode(*, explode_entries: bool, split_entries: bool) -> None:
    if explode_entries and not split_entries:
        msg = "explode_entries=True requires split_entries=True."
        raise ValueError(msg)


def _pending_source_keys(
    old_keys: Sequence[str],
    new_env_path: Path,
    *,
    skip_existing: bool,
    explode_entries: bool,
    explode_key_sep: str,
) -> tuple[list[str], int]:
    """Return the source keys still to process, and how many are already done.

    An exploded rebuild wrote its records under ``<source key><sep><sub-entry
    key>``, so resumability is decided on the source key each of them came from
    rather than on the destination key itself.
    """
    parsed = set(extract_lmdb_keys(new_env_path))
    if explode_entries:
        parsed = {explode_source_key(key, explode_key_sep) for key in parsed}
    pending = [key for key in old_keys if not skip_existing or key not in parsed]
    return pending, len(parsed)


def _emit_rebuilt(
    key: str,
    rebuilt: Mapping[str, Any],
    *,
    writer: WriterHooks,
    explode_entries: bool,
    explode_key_sep: str,
    explode_value_key: str | None,
) -> _RebuildResult | list[_RebuildResult] | None:
    """Turn one recipe result into the LMDB record(s) it should be written as.

    Without ``explode_entries`` this is the 1:1 rebuild: one source key in, one
    destination record out. With it, a ``split_entries`` result -- a mapping of
    sub-entry key -> recipe outputs -- fans out into one record per sub-entry,
    keyed ``<source key><sep><sub-entry key>``. ``explode_value_key`` unwraps a
    single recipe target so the stored value is the payload itself rather than a
    one-key envelope around it.
    """
    if not explode_entries:
        return key.encode(), encode_output(dict(rebuilt), writer=writer), None
    records: list[_RebuildResult] = []
    for sub_key, sub_value in rebuilt.items():
        payload = sub_value[explode_value_key] if explode_value_key else sub_value
        if payload is None:
            continue
        records.append(
            (
                f"{key}{explode_key_sep}{sub_key}".encode(),
                encode_output(payload, writer=writer),
                None,
            ),
        )
    return records or None


def _rebuild_worker(key: str) -> _RebuildResult | list[_RebuildResult] | None:
    state = _REBUILD_STATE
    # Key the handle cache by PID: an LMDB env handle is NOT safe to share across
    # a fork (workers would collide on the same reader-lock-table slot and stall).
    # The parent (serial test) and every forked worker each open their own.
    _ensure_worker_gc()
    pid = os.getpid()
    env = _REBUILD_ENV.get(pid)
    if env is None:
        env = _require_lmdb().open(
            str(state["old_env_path"]),
            readonly=True,
            lock=False,
            readahead=False,
        )
        _REBUILD_ENV[pid] = env
    try:
        with env.begin() as txn:
            raw_value = txn.get(key.encode())
        if raw_value is None:
            return None
        data = decode_payload(bytes(raw_value), reader=state["reader"])
        rebuilt = rebuild_entry(
            recipe=state["recipe"],
            data=data,
            metadata_dict=state["metadata_dict"],
            parameter_dict=state["parameter_dict"],
            key_transform=state["key_transform"],
            split_entries=state["split_entries"],
            extra_kwargs=state["extra_kwargs"],
        )
        if rebuilt is None:
            return None
        return _emit_rebuilt(
            key,
            rebuilt,
            writer=state["writer"],
            explode_entries=state["explode_entries"],
            explode_key_sep=state["explode_key_sep"],
            explode_value_key=state["explode_value_key"],
        )
    except Exception as error:  # noqa: BLE001
        return key.encode(), b"", error


def build_lmdb(
    *data_list: Path,
    env_path: Path,
    recipe: Path,
    inputs: Mapping[str, Any] | None = None,
    metadata_recipe: Path | None = None,
    metadata_input: Mapping[str, Any] | None = None,
    reader: ReaderHooks | None = None,
    writer: WriterHooks | None = None,
    loader: LoadFunc | None = None,
    key_transform: TransformFunc | None = None,
    serializer: SerializeFunc | None = None,
    key_builder: KeyFunc = default_lmdb_key,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    map_size: int = int(1e12),
    skip_existing: bool = False,
    test_run: bool = True,
    error_mode: ErrorMode = "skip",
    ray_batch: int | None = None,
    **extra_kwargs: Any,
) -> LmdbWriteReport:
    """Build an LMDB database from raw files and a DataCooker recipe.

    The fan-out backend is **Ray, always** -- there is no backend choice. Ray's
    measured admission control replaced the fork/joblib predictive memory model
    that OOM'd on its blind spots; ``n_jobs`` is advisory (Ray sizes to node CPUs)
    and ``ray_batch`` sets keys-per-task. See ``datacooker._ray``.
    """
    resolved_reader = ReaderHooks.from_mapping(
        None,
        loader=loader,
        key_transform=key_transform,
    ).merge(reader)
    resolved_writer = WriterHooks.from_mapping(
        None,
        serializer=serializer,
    ).merge(writer)
    if resolved_writer.serializer is None:
        msg = "build_lmdb requires a serializer or writer.serializer."
        raise ValueError(msg)
    resolved_key_func = key_builder
    _validate_error_mode(error_mode)
    lmdb_module = _require_lmdb()
    Path(env_path).parent.mkdir(parents=True, exist_ok=True)
    env = lmdb_module.open(str(env_path), map_size=int(map_size))
    metadata_dict = load_metadata(
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        extra_kwargs=extra_kwargs,
    )
    base_inputs = dict(inputs or {})
    source_paths = list(data_list)
    pending_paths = (
        filter_pending_lmdb_paths(source_paths, env_path, key_func=resolved_key_func)
        if skip_existing
        else source_paths
    )
    skipped_existing = len(source_paths) - len(pending_paths)

    # Ray is DataCooker's only build fan-out backend: the fork/joblib predictive
    # memory model (n_jobs x E x max_bytes) OOM'd on its blind spots (a deep MSA's
    # dense expansion, the parent chunk buffer), so within-node concurrency is now
    # Ray's MEASURED admission control -- no per-DB n_jobs/E/shard tuning. The
    # top-level ``_build_worker`` reads its recipe/metadata/loader from
    # ``_BUILD_STATE``, shipped to each Ray worker once via the plasma object store.
    build_state = {
        "key_func": resolved_key_func,
        "base_inputs": base_inputs,
        "metadata_dict": metadata_dict,
        "recipe": recipe,
        "loader": resolved_reader.loader,
        "key_transform": resolved_reader.key_transform,
        "writer": resolved_writer,
        "extra_kwargs": dict(extra_kwargs),
    }
    install_build_state(build_state)
    run_serial_test(_build_worker, pending_paths, enabled=test_run)
    stats = _parallel_write(
        env=env,
        items=pending_paths,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        process_item=_build_worker,
        error_mode=error_mode,
        ray_spec={
            "worker_qualname": "datacooker.lmdb.core._build_worker",
            "installer_qualname": "datacooker.lmdb.core.install_build_state",
            "state": build_state,
            **({"batch": ray_batch} if ray_batch else {}),
        },
    )
    env.close()
    return LmdbWriteReport(
        attempted=len(source_paths),
        written=stats.written,
        skipped_existing=skipped_existing,
        skipped_empty=stats.skipped_empty,
        failed=len(stats.failed_keys),
        failed_keys=stats.failed_keys,
    )


def rebuild_lmdb(
    old_env_path: Path,
    new_env_path: Path,
    recipe: Path,
    parameters: Mapping[str, Any] | None = None,
    metadata_recipe: Path | None = None,
    metadata_input: Mapping[str, Any] | None = None,
    reader: ReaderHooks | None = None,
    writer: WriterHooks | None = None,
    adapter: ConvertFunc | None = None,
    key_transform: TransformFunc | None = None,
    serializer: SerializeFunc | None = None,
    deserializer: DeserializeFunc | None = None,
    split_entries: bool = False,
    explode_entries: bool = False,
    explode_key_sep: str = "_",
    explode_value_key: str | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    map_size: int = int(1e12),
    skip_existing: bool = True,
    test_run: bool = True,
    error_mode: ErrorMode = "skip",
    shard_idx: int | None = None,
    n_shards: int = 1,
    key_list: Path | str | list | None = None,
    ray_batch: int | None = None,
    **extra_kwargs: Any,
) -> LmdbWriteReport:
    """Rebuild an LMDB database from an existing one through a DataCooker recipe.

    When ``shard_idx`` is set, only the source keys whose stable hash maps to
    that shard are processed (a disjoint partition into ``n_shards``), so the
    rebuild can be fanned out across nodes and the shard outputs merged. Point
    each shard at its own ``new_env_path`` (e.g. via a ``_shard{idx}`` suffix).

    The fan-out backend is **Ray, always** -- no backend choice. Ray self-regulates
    concurrency by measured node memory (admission control), which replaced the
    fork/joblib predictive model that OOM'd on deep MSAs; ``n_jobs`` is advisory.
    Cross-node sharding (``shard_idx`` / ``n_shards``) still applies -- Ray governs
    only the within-node fan-out. See ``datacooker._ray``.

    ``explode_entries`` turns the rebuild from 1:1 into 1:N -- it requires
    ``split_entries``, and writes one destination record per sub-entry keyed
    ``<source key><explode_key_sep><sub-entry key>`` instead of one record
    holding every sub-entry. Use it when the sub-entry, not the source entry, is
    the unit a consumer wants to read (each read then decompresses only its own
    sub-entry). ``explode_value_key`` names a single recipe target to unwrap, so
    the stored value is the payload rather than a one-key envelope. Source keys
    must not contain ``explode_key_sep``, since ``skip_existing`` recovers a
    destination record's source key by splitting on its first occurrence.
    """
    _validate_explode(explode_entries=explode_entries, split_entries=split_entries)
    resolved_reader = ReaderHooks.from_mapping(
        None,
        adapter=adapter,
        deserializer=deserializer,
        key_transform=key_transform,
    ).merge(reader)
    resolved_writer = WriterHooks.from_mapping(
        None,
        serializer=serializer,
    ).merge(writer)
    if resolved_reader.deserializer is None:
        msg = "rebuild_lmdb requires a deserializer or reader.deserializer."
        raise ValueError(msg)
    if resolved_writer.serializer is None:
        msg = "rebuild_lmdb requires a serializer or writer.serializer."
        raise ValueError(msg)
    _validate_error_mode(error_mode)
    lmdb_module = _require_lmdb()
    Path(new_env_path).parent.mkdir(parents=True, exist_ok=True)
    new_env = lmdb_module.open(str(new_env_path), map_size=int(map_size))
    metadata_dict = load_metadata(
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        extra_kwargs=extra_kwargs,
    )
    parameter_dict = dict(parameters or {})

    old_keys = extract_lmdb_keys(old_env_path)
    if key_list is not None:
        # Process only an explicit key set (e.g. a size tier's key-list), instead
        # of a hash slice -- lets giants run in their own low-concurrency lane.
        if isinstance(key_list, (str, Path)):
            with Path(key_list).open(encoding="utf-8") as f:
                wanted = {ln.strip() for ln in f if ln.strip()}
        else:
            wanted = {str(k) for k in key_list}
        old_keys = [k for k in old_keys if k in wanted]
    if shard_idx is not None:
        old_keys = shard_items_by_key(
            old_keys,
            key_func=lambda key: key,
            node_rank=shard_idx,
            node_count=n_shards,
        )
    pending_keys, n_parsed = _pending_source_keys(
        old_keys,
        new_env_path,
        skip_existing=skip_existing,
        explode_entries=explode_entries,
        explode_key_sep=explode_key_sep,
    )
    skipped_existing = len(old_keys) - len(pending_keys)
    logger.info("Already parsed %d entries into %s.", n_parsed, new_env_path)
    logger.info("To be parsed %d entries.", len(pending_keys))

    # Ray-only fan-out (see build_lmdb): the top-level ``_rebuild_worker`` reads its
    # recipe/metadata/reader from ``_REBUILD_STATE``, shipped to each Ray worker once
    # via the plasma object store. Concurrency is Ray's measured admission control,
    # not a predicted n_jobs -- this is the path that stopped the deep-MSA OOMs.
    rebuild_state = {
        "old_env_path": old_env_path,
        "reader": resolved_reader,
        "recipe": recipe,
        "metadata_dict": metadata_dict,
        "parameter_dict": parameter_dict,
        "key_transform": resolve_key_transform(reader=resolved_reader),
        "split_entries": split_entries,
        "explode_entries": explode_entries,
        "explode_key_sep": explode_key_sep,
        "explode_value_key": explode_value_key,
        "writer": resolved_writer,
        "extra_kwargs": dict(extra_kwargs),
    }
    install_rebuild_state(rebuild_state)
    run_serial_test(_rebuild_worker, pending_keys, enabled=test_run)
    stats = _parallel_write(
        env=new_env,
        items=pending_keys,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        process_item=_rebuild_worker,
        error_mode=error_mode,
        ray_spec={
            "worker_qualname": "datacooker.lmdb.core._rebuild_worker",
            "installer_qualname": "datacooker.lmdb.core.install_rebuild_state",
            "state": rebuild_state,
            **({"batch": ray_batch} if ray_batch else {}),
        },
    )
    new_env.close()
    return LmdbWriteReport(
        attempted=len(old_keys),
        written=stats.written,
        skipped_existing=skipped_existing,
        skipped_empty=stats.skipped_empty,
        failed=len(stats.failed_keys),
        failed_keys=stats.failed_keys,
    )


def suffix_key_map(suffix: str) -> KeyMapFunc:
    """Return a key map that appends a constant suffix to every key."""

    def _map(key: str) -> str:
        return key + suffix

    return _map


def rekey_lmdb(
    old_env_path: Path,
    new_env_path: Path,
    *,
    key_map: KeyMapFunc | None = None,
    map_size: int = int(1e12),
    commit_every: int = 20_000,
    skip_existing: bool = True,
    on_duplicate: DuplicatePolicy = "error",
    append: bool | None = None,
) -> LmdbWriteReport:
    """Copy an LMDB to a new database, transforming keys but copying values as-is.

    Unlike :func:`rebuild_lmdb`, no deserialize/recipe/serialize round-trip is
    performed: each raw value is written byte-for-byte unchanged, making this a
    fast, lossless, key-only rewrite. ``key_map`` maps every existing key to its
    new key; return ``None`` to drop the entry. ``key_map=None`` copies keys
    unchanged (a plain copy).

    When several source keys map to the same new key, ``on_duplicate`` decides:
    ``"overwrite"`` (last value wins), ``"skip"`` (first value wins) or
    ``"error"`` (raise).

    The copy streams through a single writer -- it is write-I/O bound, so worker
    parallelism would not help. When the transformed keys are produced in
    strictly ascending order into a fresh database, the fast ``MDB_APPEND`` path
    is used automatically; pass ``append=False`` to force ordinary puts or
    ``append=True`` to require the fast path (raising if a key is out of order).
    ``skip_existing`` makes the run resumable by skipping keys already present in
    the destination (which also disables the auto ``MDB_APPEND`` path).

    Parameters
    ----------
    old_env_path
        Source LMDB directory.
    new_env_path
        Destination LMDB directory (created if absent, appended to on resume).
    key_map
        Callable mapping an existing key to its new key, or ``None`` to drop the
        entry. ``None`` (the default) copies keys unchanged.
    map_size
        Destination map size in bytes.
    commit_every
        Number of writes per write transaction.
    skip_existing
        Skip source entries whose new key already exists in the destination.
    on_duplicate
        Policy when two source keys collide on one new key.
    append
        Force (``True``) or disable (``False``) the ``MDB_APPEND`` fast path.
        ``None`` selects it automatically for a fresh, ascending-key run.

    Returns
    -------
    LmdbWriteReport
        Counts of attempted/written/skipped/duplicate entries.
    """
    _validate_duplicate_policy(on_duplicate)
    lmdb_module = _require_lmdb()

    existing_keys = set(extract_lmdb_keys(new_env_path)) if skip_existing else set()
    # MDB_APPEND needs an empty database and strictly ascending keys, so it is
    # disabled whenever we resume into a non-empty destination.
    can_append = append if append is not None else not existing_keys

    old_env = lmdb_module.open(
        str(old_env_path),
        readonly=True,
        lock=False,
        readahead=True,
    )
    Path(new_env_path).parent.mkdir(parents=True, exist_ok=True)
    new_env = lmdb_module.open(str(new_env_path), map_size=int(map_size))
    try:
        stats = _stream_rekey(
            old_env=old_env,
            new_env=new_env,
            new_env_path=new_env_path,
            key_map=key_map,
            existing_keys=existing_keys,
            commit_every=commit_every,
            on_duplicate=on_duplicate,
            can_append=can_append,
        )
    finally:
        new_env.sync()
        new_env.close()
        old_env.close()

    return LmdbWriteReport(
        attempted=stats.attempted,
        written=stats.written,
        skipped_existing=stats.skipped_existing,
        skipped_empty=stats.skipped_empty,
        duplicate_keys=stats.duplicate_keys,
    )


def _raise_duplicate_key(new_key: str, old_key: str, dst: Path) -> None:
    msg = (
        f"Duplicate new key '{new_key}' (from '{old_key}') "
        f"while re-keying into '{dst}'."
    )
    raise KeyError(msg)


def _stream_rekey(
    *,
    old_env: Any,
    new_env: Any,
    new_env_path: Path,
    key_map: KeyMapFunc | None,
    existing_keys: set[str],
    commit_every: int,
    on_duplicate: DuplicatePolicy,
    can_append: bool,
) -> _RekeyStats:
    """Stream one LMDB into another with transformed keys and verbatim values."""
    track_duplicates = key_map is not None
    stats = _RekeyStats()
    seen: set[str] = set()
    last_key: bytes | None = None
    pending = 0
    txn = new_env.begin(write=True)
    try:
        with old_env.begin(buffers=True) as read_txn:
            for raw_key, raw_value in read_txn.cursor():
                stats.attempted += 1
                old_key = bytes(raw_key).decode()
                new_key = old_key if key_map is None else key_map(old_key)
                if new_key is None:
                    stats.skipped_empty += 1
                    continue
                if new_key in existing_keys:
                    stats.skipped_existing += 1
                    continue
                if track_duplicates and new_key in seen:
                    stats.duplicate_keys += 1
                    if on_duplicate == "error":
                        _raise_duplicate_key(new_key, old_key, new_env_path)
                    if on_duplicate == "skip":
                        continue
                    # "overwrite": fall through; a plain put replaces the value.
                if track_duplicates:
                    seen.add(new_key)
                encoded = new_key.encode()
                use_append = can_append and (last_key is None or encoded > last_key)
                if can_append and not use_append:
                    # Keys are not globally ascending; drop to ordinary puts.
                    can_append = False
                txn.put(encoded, bytes(raw_value), append=use_append)
                stats.written += 1
                if last_key is None or encoded > last_key:
                    last_key = encoded
                pending += 1
                if pending >= commit_every:
                    txn.commit()
                    txn = new_env.begin(write=True)
                    pending = 0
        txn.commit()
    except BaseException:
        txn.abort()
        raise
    return stats


def merge_lmdb_shards(
    shard_paths: Sequence[Path],
    merged_env_path: Path,
    *,
    map_size: int = int(1e12),
    overwrite: bool = False,
    on_duplicate: DuplicatePolicy = "overwrite",
    value_deserializer: Callable[[bytes], dict] | None = None,
    value_serializer: Callable[[dict], bytes] | None = None,
) -> LmdbWriteReport:
    """Merge multiple LMDB shard databases into one database.

    ``on_duplicate="union"`` combines colliding records by deep-merging their
    deserialized ``dict`` values (``{**existing, **incoming}``) instead of
    overwriting. It requires ``value_deserializer`` and ``value_serializer``
    so records can be decoded, merged, and re-encoded; later shards win on
    per-subkey conflicts.
    """
    _validate_duplicate_policy(on_duplicate)
    if on_duplicate == "union" and (
        value_deserializer is None or value_serializer is None
    ):
        msg = (
            "on_duplicate='union' requires both value_deserializer and "
            "value_serializer to decode/encode records for merging."
        )
        raise ValueError(msg)
    lmdb_module = _require_lmdb()
    if merged_env_path.exists() and not overwrite:
        msg = f"{merged_env_path} already exists. Use overwrite=True to replace it."
        raise FileExistsError(msg)

    merged_env = lmdb_module.open(str(merged_env_path), map_size=map_size)
    total_keys = 0
    duplicate_keys = 0

    # K-way merge of the shards' key-sorted cursors -> one globally key-sorted
    # stream. New keys go in with MDB_APPEND (sequential bulk fill of the B-tree:
    # far faster than random puts, no page churn). Duplicates land CONSECUTIVELY
    # in the sorted stream, so they are detected by comparing to the last key --
    # no O(n) ``seen_keys`` set (which itself OOMs on 16M-entry merges). The write
    # txn is flushed every ``flush_bytes`` so dirty pages stay bounded regardless
    # of value size; MDB_APPEND stays valid across commits since keys only ascend.
    flush_bytes = 2 * 1024**3
    pending_bytes = 0

    shard_envs = [
        lmdb_module.open(str(p), readonly=True, lock=False) for p in shard_paths
    ]
    shard_txns = [e.begin() for e in shard_envs]
    streams = [t.cursor().iternext(keys=True, values=True) for t in shard_txns]
    merged_stream = heapq.merge(*streams, key=lambda kv: bytes(kv[0]))

    merged_txn = merged_env.begin(write=True)
    last_key: bytes | None = None
    try:
        for key, value in merged_stream:
            k = bytes(key)
            v = bytes(value)
            if k == last_key:  # duplicate (consecutive in the sorted stream)
                duplicate_keys += 1
                if on_duplicate == "skip":
                    continue
                if on_duplicate == "error":
                    msg = (
                        f"Duplicate LMDB key '{k.decode()}' encountered while "
                        f"merging shards into '{merged_env_path}'."
                    )
                    raise KeyError(msg)
                if on_duplicate == "union":
                    existing = merged_txn.get(k)
                    if existing is not None:
                        v = value_serializer({
                            **value_deserializer(bytes(existing)),
                            **value_deserializer(v),
                        })
                merged_txn.put(k, v)  # append=False: overwrite the just-written key
            else:
                merged_txn.put(k, v, append=True)  # sequential bulk load
                total_keys += 1
                last_key = k
            pending_bytes += len(v)
            if pending_bytes >= flush_bytes:
                merged_txn.commit()
                merged_txn = merged_env.begin(write=True)
                pending_bytes = 0
        merged_txn.commit()
    except BaseException:
        merged_txn.abort()
        raise
    finally:
        for txn in shard_txns:
            txn.abort()
        for env in shard_envs:
            env.close()
        merged_env.sync()
        merged_env.close()

    logger.info("[Done] Merged %d shards into %s", len(shard_paths), merged_env_path)
    logger.info("Total keys merged: %d", total_keys)
    return LmdbWriteReport(
        attempted=total_keys + duplicate_keys,
        written=total_keys,
        duplicate_keys=duplicate_keys,
    )


def _get_env(
    path: Path,
    *,
    readonly: bool,
    lock: bool,
) -> Any:
    lmdb_module = _require_lmdb()
    # Key by PID: an LMDB env must not be shared across a fork (its internal
    # mutexes deadlock). Under the "multiprocessing" backend, forked workers
    # then miss the parent's cached env and open their own handle.
    key = (str(path), readonly, lock, os.getpid())
    env = _ENV_CACHE.get(key)
    if env is None:
        env = lmdb_module.open(str(path), readonly=readonly, lock=lock)
        _ENV_CACHE[key] = env
    return env


def _parallel_write(
    *,
    env: Any,
    items: Sequence[Any],
    chunk_size: int,
    n_jobs: int,
    process_item: Any,
    error_mode: ErrorMode = "skip",
    ray_spec: Mapping[str, Any] | None = None,
) -> _WriteStats:
    # build/rebuild fan out over Ray only (see build_lmdb); ``process_item`` is the
    # worker Ray runs per key, ``n_jobs`` is advisory (Ray sizes to node CPUs).
    _validate_error_mode(error_mode)
    written = 0
    skipped_empty = 0
    failed_keys: list[str] = []
    for chunk_results in iter_parallel_chunks(
        items=items,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        process_item=process_item,
        backend="ray",
        ray_spec=ray_spec,
        log_label="Processing files",
    ):
        pending_writes: list[tuple[bytes, bytes]] = []
        with env.begin(write=True) as txn:
            for result in chunk_results:
                if result is None:
                    skipped_empty += 1
                    continue
                # A 1:N (exploded) rebuild returns a list of records for one item;
                # every other producer returns a single record.
                records = result if isinstance(result, list) else [result]
                for key, payload, error in records:
                    if error is not None:
                        key_text = key.decode()
                        failed_keys.append(key_text)
                        if error_mode == "raise":
                            msg = f"Error processing LMDB key '{key_text}': {error}"
                            raise RuntimeError(msg) from error
                        logger.error("Error processing %s: %s", key_text, error)
                        continue
                    pending_writes.append((key, payload))
            for key, payload in pending_writes:
                txn.put(key, payload)
                written += 1
    return _WriteStats(
        written=written,
        skipped_empty=skipped_empty,
        failed_keys=tuple(failed_keys),
    )


def _require_lmdb() -> Any:
    try:
        return _load_lmdb_module()
    except ModuleNotFoundError as exc:
        msg = (
            "LMDB support requires the optional 'lmdb' dependency. "
            "Install DataCooker with LMDB support before using datacooker.lmdb."
        )
        raise ModuleNotFoundError(msg) from exc


@lru_cache(maxsize=1)
def _load_lmdb_module() -> Any:
    return import_module("lmdb")


@dataclass(frozen=True, slots=True)
class _WriteStats:
    written: int
    skipped_empty: int
    failed_keys: tuple[str, ...]


@dataclass(slots=True)
class _RekeyStats:
    attempted: int = 0
    written: int = 0
    skipped_existing: int = 0
    skipped_empty: int = 0
    duplicate_keys: int = 0


def _validate_error_mode(error_mode: str) -> None:
    if error_mode in {"skip", "raise"}:
        return
    msg = f"Unsupported error_mode '{error_mode}'. Expected 'skip' or 'raise'."
    raise ValueError(msg)


def _validate_duplicate_policy(policy: str) -> None:
    if policy in {"overwrite", "skip", "error", "union"}:
        return
    msg = (
        f"Unsupported on_duplicate policy '{policy}'. "
        "Expected 'overwrite', 'skip', 'error', or 'union'."
    )
    raise ValueError(msg)
