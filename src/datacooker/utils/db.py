"""Generic LMDB helpers for recipe-driven database workflows."""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from functools import lru_cache
from importlib import import_module
from pathlib import Path
from typing import Any

from datacooker.api import parse_dict, parse_file
from datacooker.protocols import (
    ConvertFunc,
    DeserializeFunc,
    LoadFunc,
    SerializeFunc,
    TransformFunc,
)

logger = logging.getLogger(__name__)

_ENV_CACHE: dict[tuple[str, bool, bool], Any] = {}


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
    inputs: Mapping[str, Any] | None = None,
    metadata_recipe: Path | None = None,
    metadata_input: Mapping[str, Any] | None = None,
    load_func: LoadFunc | None = None,
    transform_func: TransformFunc | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    map_size: int = int(1e12),
    test_run: bool = True,
    **extra_kwargs: Any,
) -> None:
    """Build an LMDB database from raw files and a DataCooker recipe."""
    lmdb_module = _require_lmdb()
    _require_joblib()
    env = lmdb_module.open(str(env_path), map_size=int(map_size))
    metadata_dict = _load_metadata(
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        extra_kwargs=extra_kwargs,
    )
    base_inputs = dict(inputs or {})

    def _process_file(data_file: Path) -> tuple[bytes, bytes, Exception | None]:
        key = data_file.name.split(".")[0]
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

    _run_serial_test(_process_file, data_list, enabled=test_run)
    _parallel_write(
        env=env,
        items=list(data_list),
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        process_item=_process_file,
    )
    env.close()


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
    test_run: bool = True,
    **extra_kwargs: Any,
) -> None:
    """Rebuild an LMDB database from an existing one through a DataCooker recipe."""
    lmdb_module = _require_lmdb()
    _require_joblib()
    new_env = lmdb_module.open(str(new_env_path), map_size=int(map_size))
    metadata_dict = _load_metadata(
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        extra_kwargs=extra_kwargs,
    )
    parameter_dict = dict(parameters or {})

    old_keys = extract_lmdb_keys(old_env_path)
    parsed_keys = set(extract_lmdb_keys(new_env_path))
    pending_keys = [key for key in old_keys if key not in parsed_keys]
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
            rebuilt = _rebuild_entry(
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

    _run_serial_test(_process_key, pending_keys, enabled=test_run)
    _parallel_write(
        env=new_env,
        items=pending_keys,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        process_item=_process_key,
        prefer="threads",
    )
    old_env.close()
    new_env.close()


def merge_lmdb_shards(
    shard_paths: Sequence[Path],
    merged_env_path: Path,
    *,
    map_size: int = int(1e12),
    overwrite: bool = False,
) -> None:
    """Merge multiple LMDB shard databases into one database."""
    lmdb_module = _require_lmdb()
    if merged_env_path.exists() and not overwrite:
        msg = f"{merged_env_path} already exists. Use overwrite=True to replace it."
        raise FileExistsError(msg)

    merged_env = lmdb_module.open(str(merged_env_path), map_size=map_size)
    total_keys = 0
    try:
        for shard_path in shard_paths:
            logger.info("Merging shard: %s", shard_path)
            shard_env = lmdb_module.open(str(shard_path), readonly=True, lock=False)
            try:
                with shard_env.begin() as shard_txn, merged_env.begin(
                    write=True
                ) as merged_txn:
                    for key, value in shard_txn.cursor():
                        merged_txn.put(key, value)
                        total_keys += 1
            finally:
                shard_env.close()
    finally:
        merged_env.sync()
        merged_env.close()

    logger.info("[Done] Merged %d shards into %s", len(shard_paths), merged_env_path)
    logger.info("Total keys merged: %d", total_keys)


def _load_metadata(
    *,
    metadata_recipe: Path | None,
    metadata_input: Mapping[str, Any] | None,
    extra_kwargs: Mapping[str, Any],
) -> dict[str, Any]:
    if metadata_recipe is None:
        return {}
    if metadata_input is None:
        msg = "metadata_input must be provided if metadata_recipe is specified."
        raise ValueError(msg)
    return parse_dict(
        recipe_path=metadata_recipe,
        datadict=metadata_input,
        **extra_kwargs,
    )


def _run_serial_test(
    process_item: Any,
    items: Sequence[Any],
    *,
    enabled: bool,
) -> None:
    if not enabled or not items:
        return

    for item in items[:10]:
        result = process_item(item)
        if result is not None and result[2] is not None:
            raise result[2]


def _parallel_write(
    *,
    env: Any,
    items: Sequence[Any],
    chunk_size: int,
    n_jobs: int,
    process_item: Any,
    prefer: str = "processes",
) -> None:
    parallel, delayed = _require_joblib()
    for start in range(0, len(items), chunk_size):
        logger.info(
            "Processing files %d to %d / %d",
            start,
            min(start + chunk_size, len(items)),
            len(items),
        )
        chunk = items[start : start + chunk_size]
        results = parallel(n_jobs=n_jobs, verbose=10, prefer=prefer)(
            delayed(process_item)(item) for item in chunk
        )

        with env.begin(write=True) as txn:
            for result in results:
                if result is None:
                    continue
                key, payload, error = result
                if error is not None:
                    logger.error("Error processing %s: %s", key.decode(), error)
                    continue
                txn.put(key, payload)


def _require_lmdb() -> Any:
    try:
        return _load_lmdb_module()
    except ModuleNotFoundError as exc:
        msg = (
            "LMDB support requires the optional 'lmdb' dependency. "
            "Install DataCooker with LMDB support before using datacooker.utils.db."
        )
        raise ModuleNotFoundError(msg) from exc


def _require_joblib() -> tuple[Any, Any]:
    _require_lmdb()
    try:
        joblib_module = _load_joblib_module()
    except ModuleNotFoundError as exc:
        msg = (
            "LMDB workflow helpers require the optional 'joblib' dependency. "
            "Install DataCooker with joblib support before using datacooker.utils.db."
        )
        raise ModuleNotFoundError(msg) from exc
    return joblib_module.Parallel, joblib_module.delayed


@lru_cache(maxsize=1)
def _load_lmdb_module() -> Any:
    return import_module("lmdb")


@lru_cache(maxsize=1)
def _load_joblib_module() -> Any:
    return import_module("joblib")


def _rebuild_entry(
    *,
    recipe: Path,
    data: dict[str, Any],
    metadata_dict: Mapping[str, Any],
    parameter_dict: Mapping[str, Any],
    transform_func: TransformFunc | None,
    split_entries: bool,
    extra_kwargs: Mapping[str, Any],
) -> dict[str, Any] | None:
    if not split_entries:
        datadict = dict(data)
        datadict.update(metadata_dict)
        datadict.update(parameter_dict)
        rebuilt = parse_dict(
            recipe_path=recipe,
            datadict=datadict,
            transform_func=transform_func,
            **extra_kwargs,
        )
        if all(value is None for value in rebuilt.values()):
            return None
        return rebuilt

    output_dict: dict[str, Any] = {}
    for data_key, inner_dict in data.items():
        if not isinstance(inner_dict, Mapping):
            msg = (
                "split_entries=True requires each LMDB entry to deserialize into a "
                "mapping of sub-entry mappings."
            )
            raise TypeError(msg)
        datadict = dict(inner_dict)
        datadict.update(metadata_dict)
        datadict.update(parameter_dict)
        rebuilt = parse_dict(
            recipe_path=recipe,
            datadict=datadict,
            transform_func=transform_func,
            **extra_kwargs,
        )
        if all(value is None for value in rebuilt.values()):
            continue
        output_dict[data_key] = rebuilt

    return output_dict or None
