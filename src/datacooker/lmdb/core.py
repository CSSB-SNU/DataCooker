"""Core LMDB workflow helpers for recipe-driven database pipelines."""

from __future__ import annotations

import logging
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
    LoadFunc,
    SerializeFunc,
    TransformFunc,
)

logger = logging.getLogger(__name__)

_ENV_CACHE: dict[tuple[str, bool, bool], Any] = {}
ErrorMode = Literal["skip", "raise"]
DuplicatePolicy = Literal["overwrite", "skip", "error"]


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
    deserialize: DeserializeFunc,
) -> dict[str, Any]:
    """Read a decoded entry from an LMDB database."""
    value = read_lmdb_raw(env_path, key)
    if value is None:
        msg = f"Key '{key}' not found in LMDB database at '{env_path}'."
        raise KeyError(msg)
    return deserialize(value)


def build_lmdb(
    *data_list: Path,
    env_path: Path,
    recipe: Path,
    serialize: SerializeFunc,
    key_func: Callable[[Path], str] = default_lmdb_key,
    inputs: Mapping[str, Any] | None = None,
    metadata_recipe: Path | None = None,
    metadata_input: Mapping[str, Any] | None = None,
    load_func: LoadFunc | None = None,
    transform_func: TransformFunc | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    map_size: int = int(1e12),
    skip_existing: bool = False,
    test_run: bool = True,
    error_mode: ErrorMode = "skip",
    **extra_kwargs: Any,
) -> LmdbWriteReport:
    """Build an LMDB database from raw files and a DataCooker recipe."""
    _validate_error_mode(error_mode)
    lmdb_module = _require_lmdb()
    env = lmdb_module.open(str(env_path), map_size=int(map_size))
    metadata_dict = load_metadata(
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        extra_kwargs=extra_kwargs,
    )
    base_inputs = dict(inputs or {})
    source_paths = list(data_list)
    pending_paths = (
        filter_pending_lmdb_paths(source_paths, env_path, key_func=key_func)
        if skip_existing
        else source_paths
    )
    skipped_existing = len(source_paths) - len(pending_paths)

    def _process_file(data_file: Path) -> tuple[bytes, bytes, Exception | None]:
        key = key_func(data_file)
        try:
            data_dict = dict(base_inputs)
            data_dict.update(metadata_dict)
            output = parse_file(
                recipe_path=recipe,
                file_path=data_file,
                load_func=load_func,
                inputs=data_dict,
                transform_func=transform_func,
                **extra_kwargs,
            )
            return key.encode(), serialize(output), None
        except Exception as error:  # noqa: BLE001
            return key.encode(), b"", error

    run_serial_test(_process_file, pending_paths, enabled=test_run)
    stats = _parallel_write(
        env=env,
        items=pending_paths,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        process_item=_process_file,
        error_mode=error_mode,
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
    serialize: SerializeFunc,
    deserialize: DeserializeFunc,
    parameters: Mapping[str, Any] | None = None,
    metadata_recipe: Path | None = None,
    metadata_input: Mapping[str, Any] | None = None,
    convert_func: ConvertFunc | None = None,
    transform_func: TransformFunc | None = None,
    split_entries: bool = False,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    map_size: int = int(1e12),
    skip_existing: bool = True,
    test_run: bool = True,
    error_mode: ErrorMode = "skip",
    **extra_kwargs: Any,
) -> LmdbWriteReport:
    """Rebuild an LMDB database from an existing one through a DataCooker recipe."""
    _validate_error_mode(error_mode)
    lmdb_module = _require_lmdb()
    new_env = lmdb_module.open(str(new_env_path), map_size=int(map_size))
    metadata_dict = load_metadata(
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        extra_kwargs=extra_kwargs,
    )
    parameter_dict = dict(parameters or {})

    old_keys = extract_lmdb_keys(old_env_path)
    parsed_keys = set(extract_lmdb_keys(new_env_path))
    pending_keys = [key for key in old_keys if not skip_existing or key not in parsed_keys]
    skipped_existing = len(old_keys) - len(pending_keys)
    logger.info("Already parsed %d entries into %s.", len(parsed_keys), new_env_path)
    logger.info("To be parsed %d entries.", len(pending_keys))

    old_env = lmdb_module.open(str(old_env_path), readonly=True, lock=False)

    def _process_key(key: str) -> tuple[bytes, bytes, Exception | None] | None:
        try:
            with old_env.begin() as txn:
                raw_value = txn.get(key.encode())
            if raw_value is None:
                return None

            data = deserialize(bytes(raw_value))
            data = convert_func(data) if convert_func is not None else data
            rebuilt = rebuild_entry(
                recipe=recipe,
                data=data,
                metadata_dict=metadata_dict,
                parameter_dict=parameter_dict,
                transform_func=transform_func,
                split_entries=split_entries,
                extra_kwargs=extra_kwargs,
            )
            if rebuilt is None:
                return None
            return key.encode(), serialize(rebuilt), None
        except Exception as error:  # noqa: BLE001
            return key.encode(), b"", error

    run_serial_test(_process_key, pending_keys, enabled=test_run)
    stats = _parallel_write(
        env=new_env,
        items=pending_keys,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        process_item=_process_key,
        prefer="threads",
        error_mode=error_mode,
    )
    old_env.close()
    new_env.close()
    return LmdbWriteReport(
        attempted=len(old_keys),
        written=stats.written,
        skipped_existing=skipped_existing,
        skipped_empty=stats.skipped_empty,
        failed=len(stats.failed_keys),
        failed_keys=stats.failed_keys,
    )


def merge_lmdb_shards(
    shard_paths: Sequence[Path],
    merged_env_path: Path,
    *,
    map_size: int = int(1e12),
    overwrite: bool = False,
    on_duplicate: DuplicatePolicy = "overwrite",
) -> LmdbWriteReport:
    """Merge multiple LMDB shard databases into one database."""
    _validate_duplicate_policy(on_duplicate)
    lmdb_module = _require_lmdb()
    if merged_env_path.exists() and not overwrite:
        msg = f"{merged_env_path} already exists. Use overwrite=True to replace it."
        raise FileExistsError(msg)

    merged_env = lmdb_module.open(str(merged_env_path), map_size=map_size)
    total_keys = 0
    duplicate_keys = 0
    seen_keys: set[bytes] = set()
    try:
        for shard_path in shard_paths:
            logger.info("Merging shard: %s", shard_path)
            shard_env = lmdb_module.open(str(shard_path), readonly=True, lock=False)
            try:
                with shard_env.begin() as shard_txn, merged_env.begin(
                    write=True
                ) as merged_txn:
                    for key, value in shard_txn.cursor():
                        normalized_key = bytes(key)
                        if normalized_key in seen_keys:
                            duplicate_keys += 1
                            if on_duplicate == "skip":
                                continue
                            if on_duplicate == "error":
                                key_text = normalized_key.decode()
                                msg = (
                                    f"Duplicate LMDB key '{key_text}' encountered while "
                                    f"merging shards into '{merged_env_path}'."
                                )
                                raise KeyError(msg)
                        seen_keys.add(normalized_key)
                        merged_txn.put(key, value)
                        total_keys += 1
            finally:
                shard_env.close()
    finally:
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
    key = (str(path), readonly, lock)
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
    prefer: str = "processes",
    error_mode: ErrorMode = "skip",
) -> _WriteStats:
    _validate_error_mode(error_mode)
    written = 0
    skipped_empty = 0
    failed_keys: list[str] = []
    for chunk_results in iter_parallel_chunks(
        items=items,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        process_item=process_item,
        prefer=prefer,
        log_label="Processing files",
    ):
        pending_writes: list[tuple[bytes, bytes]] = []
        with env.begin(write=True) as txn:
            for result in chunk_results:
                if result is None:
                    skipped_empty += 1
                    continue
                key, payload, error = result
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


def _validate_error_mode(error_mode: str) -> None:
    if error_mode in {"skip", "raise"}:
        return
    msg = f"Unsupported error_mode '{error_mode}'. Expected 'skip' or 'raise'."
    raise ValueError(msg)


def _validate_duplicate_policy(policy: str) -> None:
    if policy in {"overwrite", "skip", "error"}:
        return
    msg = (
        f"Unsupported on_duplicate policy '{policy}'. "
        "Expected 'overwrite', 'skip', or 'error'."
    )
    raise ValueError(msg)
