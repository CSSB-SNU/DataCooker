"""Input-boundary helpers for loading and adapting workflow inputs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from datacooker.protocols import ConvertFunc, DeserializeFunc, LoadFunc, TransformFunc


@dataclass(frozen=True, slots=True)
class ReaderHooks:
    """Collect the read-side hooks used at workflow boundaries."""

    loader: LoadFunc | None = None
    adapter: ConvertFunc | None = None
    deserializer: DeserializeFunc | None = None
    key_transform: TransformFunc | None = None

    @classmethod
    def from_mapping(
        cls,
        mapping: Mapping[str, Any] | None = None,
        *,
        loader: LoadFunc | None = None,
        adapter: ConvertFunc | None = None,
        deserializer: DeserializeFunc | None = None,
        key_transform: TransformFunc | None = None,
        load_func: LoadFunc | None = None,
        convert_func: ConvertFunc | None = None,
        deserialize: DeserializeFunc | None = None,
        transform_func: TransformFunc | None = None,
    ) -> ReaderHooks:
        """Build reader hooks from a config mapping plus compatibility aliases."""
        reader_dict = dict(mapping or {})
        return cls(
            loader=loader or load_func or reader_dict.get("loader") or reader_dict.get("load_func"),
            adapter=adapter or convert_func or reader_dict.get("adapter") or reader_dict.get("convert_func"),
            deserializer=(
                deserializer
                or deserialize
                or reader_dict.get("deserializer")
                or reader_dict.get("deserialize")
            ),
            key_transform=(
                key_transform
                or transform_func
                or reader_dict.get("key_transform")
                or reader_dict.get("transform_func")
            ),
        )


def dot_path(key: str) -> tuple[str, ...]:
    """Split a dotted key into an execution-context key path."""
    return tuple(key.split("."))


def load_inputs(
    file_path: Path,
    *,
    reader: ReaderHooks | None = None,
    loader: LoadFunc | None = None,
    base_inputs: Mapping[str, Any] | None = None,
    **extra_kwargs: Any,
) -> dict[str, Any]:
    """Load a file and merge the result with provided base inputs."""
    resolved_reader = ReaderHooks.from_mapping(None, loader=loader)
    if reader is not None:
        resolved_reader = ReaderHooks(
            loader=reader.loader or resolved_reader.loader,
            adapter=reader.adapter,
            deserializer=reader.deserializer,
            key_transform=reader.key_transform,
        )

    payload = dict(base_inputs or {})
    if resolved_reader.loader is not None:
        payload.update(resolved_reader.loader(file_path))
    else:
        payload["file_path"] = file_path
    payload.update(extra_kwargs)
    return payload


def decode_payload(
    payload: bytes,
    *,
    reader: ReaderHooks | None = None,
    deserializer: DeserializeFunc | None = None,
    adapter: ConvertFunc | None = None,
) -> dict[str, Any]:
    """Deserialize raw bytes and optionally adapt the resulting dictionary."""
    resolved_reader = ReaderHooks.from_mapping(
        None,
        deserializer=deserializer,
        adapter=adapter,
    )
    if reader is not None:
        resolved_reader = ReaderHooks(
            loader=reader.loader,
            adapter=reader.adapter or resolved_reader.adapter,
            deserializer=reader.deserializer or resolved_reader.deserializer,
            key_transform=reader.key_transform,
        )
    if resolved_reader.deserializer is None:
        msg = "decode_payload requires a deserializer."
        raise ValueError(msg)

    data = resolved_reader.deserializer(payload)
    if resolved_reader.adapter is not None:
        data = resolved_reader.adapter(data)
    return data


def resolve_key_transform(
    *,
    reader: ReaderHooks | None = None,
    key_transform: TransformFunc | None = None,
    transform_func: TransformFunc | None = None,
) -> TransformFunc | None:
    """Resolve the key transform from v2 or compatibility arguments."""
    if key_transform is not None:
        return key_transform
    if reader is not None and reader.key_transform is not None:
        return reader.key_transform
    return transform_func


def normalize_input_sequence(data_list: Any, *, output_name: str) -> list[dict[str, Any]]:
    """Validate that a split-recipe output is a sequence of mappings."""
    if data_list is None:
        msg = f"Expected '{output_name}' in split workflow output."
        raise KeyError(msg)
    if isinstance(data_list, (str, bytes, bytearray)) or not isinstance(data_list, Sequence):
        msg = (
            f"Expected '{output_name}' to be a sequence of mappings, "
            f"got {type(data_list).__name__}."
        )
        raise TypeError(msg)

    normalized: list[dict[str, Any]] = []
    for index, item in enumerate(data_list):
        if not isinstance(item, Mapping):
            msg = (
                f"Expected '{output_name}[{index}]' to be a mapping, "
                f"got {type(item).__name__}."
            )
            raise TypeError(msg)
        normalized.append(dict(item))
    return normalized
