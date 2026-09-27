"""Input-boundary helpers for loading and adapting workflow inputs."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from datacooker.protocols import ConvertFunc, DeserializeFunc, LoadFunc, TransformFunc


@dataclass(frozen=True, slots=True)
class InputHooks:
    """The file-input boundary: turn a path into a workflow input dict."""

    loader: LoadFunc | None = None


@dataclass(frozen=True, slots=True)
class PayloadHooks:
    """The payload boundary: turn stored bytes into a workflow data dict.

    ``deserializer`` performs ``bytes -> dict``; the optional ``adapter`` then
    reshapes the resulting ``dict -> dict``.
    """

    deserializer: DeserializeFunc | None = None
    adapter: ConvertFunc | None = None


@dataclass(frozen=True, slots=True)
class ReaderHooks:
    """Aggregate of the read-side hooks used across workflow boundaries.

    The individual concerns are also available as cohesive views: ``.input``
    (file loading) and ``.payload`` (byte decoding + adapting). Simple helpers
    accept those focused types directly; multi-stage pipelines (LMDB build /
    rebuild / extract) keep using the aggregate for convenience.
    """

    loader: LoadFunc | None = None
    adapter: ConvertFunc | None = None
    deserializer: DeserializeFunc | None = None
    key_transform: TransformFunc | None = None

    @property
    def input(self) -> InputHooks:
        """Return the file-input boundary view."""
        return InputHooks(loader=self.loader)

    @property
    def payload(self) -> PayloadHooks:
        """Return the payload (deserialize + adapt) boundary view."""
        return PayloadHooks(deserializer=self.deserializer, adapter=self.adapter)

    def merge(self, override: ReaderHooks | None) -> ReaderHooks:
        """Return a copy where every non-None field of ``override`` wins.

        This replaces the field-by-field merge that used to be hand-written at
        each boundary, which was easy to get subtly wrong.
        """
        if override is None:
            return self
        return ReaderHooks(
            loader=override.loader if override.loader is not None else self.loader,
            adapter=override.adapter if override.adapter is not None else self.adapter,
            deserializer=(
                override.deserializer
                if override.deserializer is not None
                else self.deserializer
            ),
            key_transform=(
                override.key_transform
                if override.key_transform is not None
                else self.key_transform
            ),
        )

    @classmethod
    def from_mapping(
        cls,
        mapping: Mapping[str, Any] | None = None,
        *,
        loader: LoadFunc | None = None,
        adapter: ConvertFunc | None = None,
        deserializer: DeserializeFunc | None = None,
        key_transform: TransformFunc | None = None,
    ) -> ReaderHooks:
        """Build reader hooks from a config mapping and/or explicit callables."""
        reader_dict = dict(mapping or {})
        return cls(
            loader=loader or reader_dict.get("loader"),
            adapter=adapter or reader_dict.get("adapter"),
            deserializer=deserializer or reader_dict.get("deserializer"),
            key_transform=key_transform or reader_dict.get("key_transform"),
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
    resolved_reader = ReaderHooks(loader=loader).merge(reader)

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
    resolved_reader = ReaderHooks(
        deserializer=deserializer,
        adapter=adapter,
    ).merge(reader)
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
) -> TransformFunc | None:
    """Resolve the key transform from an explicit value or reader hooks."""
    if key_transform is not None:
        return key_transform
    if reader is not None and reader.key_transform is not None:
        return reader.key_transform
    return None


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
