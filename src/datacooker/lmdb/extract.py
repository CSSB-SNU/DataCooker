"""LMDB-wide extraction helpers built on top of the core LMDB APIs."""

from __future__ import annotations

import logging
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from datacooker._parallel import iter_parallel_chunks, require_joblib, run_serial_test
from datacooker.api import parse_dict
from datacooker.lmdb._shared import load_metadata
from datacooker.lmdb.core import extract_lmdb_keys, read_lmdb_raw
from datacooker.protocols import ConvertFunc, DeserializeFunc, TransformFunc
from datacooker.readers import ReaderHooks, decode_payload, resolve_key_transform

logger = logging.getLogger(__name__)


def _require_deserializer(reader: ReaderHooks) -> None:
    if reader.deserializer is None:
        msg = "extract_lmdb_records requires a deserializer."
        raise ValueError(msg)


def _read_existing_payload(env_path: Path, key: str) -> bytes:
    raw_data = read_lmdb_raw(env_path, key)
    if raw_data is None:
        msg = f"Key '{key}' not found in LMDB database at '{env_path}'."
        raise KeyError(msg)
    return raw_data


def extract_lmdb_records(
    env_path: Path,
    recipe: Path,
    *,
    deserialize: DeserializeFunc | None = None,
    inputs: Mapping[str, Any] | None = None,
    metadata_recipe: Path | None = None,
    metadata_input: Mapping[str, Any] | None = None,
    reader: ReaderHooks | None = None,
    adapter: ConvertFunc | None = None,
    key_transform: TransformFunc | None = None,
    deserializer: DeserializeFunc | None = None,
    convert_func: ConvertFunc | None = None,
    transform_func: TransformFunc | None = None,
    merge_recipe: Path | None = None,
    merge_inputs: Mapping[str, Any] | None = None,
    merge_input_name: str = "data_dict",
    chunk_size: int = 100,
    n_jobs: int = -1,
    test_run: bool = True,
    prefer: str = "processes",
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Extract recipe outputs from every LMDB record and optionally merge them."""
    resolved_reader = ReaderHooks.from_mapping(
        None,
        adapter=adapter,
        deserializer=deserializer,
        key_transform=key_transform,
        convert_func=convert_func,
        deserialize=deserialize,
        transform_func=transform_func,
    )
    if reader is not None:
        resolved_reader = ReaderHooks(
            loader=reader.loader,
            adapter=reader.adapter or resolved_reader.adapter,
            deserializer=reader.deserializer or resolved_reader.deserializer,
            key_transform=reader.key_transform or resolved_reader.key_transform,
        )
    if resolved_reader.deserializer is None:
        msg = "extract_lmdb_records requires a deserializer or reader.deserializer."
        raise ValueError(msg)
    require_joblib(feature_name="datacooker.lmdb")
    if chunk_size < 1:
        msg = f"chunk_size must be >= 1, got {chunk_size}."
        raise ValueError(msg)

    key_list = extract_lmdb_keys(env_path)
    metadata_dict = load_metadata(
        metadata_recipe=metadata_recipe,
        metadata_input=metadata_input,
        extra_kwargs=extra_kwargs,
    )
    base_inputs = dict(inputs or {})
    if not key_list:
        return _merge_lmdb_record_outputs(
            {},
            merge_recipe=merge_recipe,
            merge_inputs=merge_inputs,
            merge_input_name=merge_input_name,
            extra_kwargs=extra_kwargs,
        )

    def _process_key(key: str) -> tuple[str, dict[str, Any], Exception | None]:
        try:
            _require_deserializer(resolved_reader)
            raw_data = _read_existing_payload(env_path, key)
            data = decode_payload(raw_data, reader=resolved_reader)
            datadict = {"key": key, "db_data": data}
            datadict.update(base_inputs)
            datadict.update(metadata_dict)
            result = parse_dict(
                datadict=datadict,
                recipe_path=recipe,
                key_transform=resolve_key_transform(reader=resolved_reader),
                **extra_kwargs,
            )
        except Exception as error:  # noqa: BLE001
            return key, {}, error
        return key, result, None

    run_serial_test(_process_key, key_list, enabled=test_run)

    merged_output: dict[str, Any] = {}
    for chunk_results in iter_parallel_chunks(
        items=key_list,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        process_item=_process_key,
        prefer=prefer,
        log_label="Extracting LMDB keys",
    ):
        for key, result, error in chunk_results:
            if error is not None:
                logger.error("Error extracting %s: %s", key, error)
                continue
            merged_output[key] = result

    return _merge_lmdb_record_outputs(
        merged_output,
        merge_recipe=merge_recipe,
        merge_inputs=merge_inputs,
        merge_input_name=merge_input_name,
        extra_kwargs=extra_kwargs,
    )


def _merge_lmdb_record_outputs(
    output: dict[str, Any],
    *,
    merge_recipe: Path | None,
    merge_inputs: Mapping[str, Any] | None,
    merge_input_name: str,
    extra_kwargs: Mapping[str, Any],
) -> dict[str, Any]:
    if merge_recipe is None:
        return output

    merge_datadict = {merge_input_name: output}
    if merge_inputs is not None:
        merge_datadict.update(merge_inputs)
    return parse_dict(
        recipe_path=merge_recipe,
        datadict=merge_datadict,
        **extra_kwargs,
    )
