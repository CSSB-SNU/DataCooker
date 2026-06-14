"""Output-boundary helpers for materializing workflow results."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from datacooker.protocols import ProjectFunc, SerializeFunc


@dataclass(frozen=True, slots=True)
class WriterHooks:
    """Collect the write-side hooks used at workflow boundaries."""

    serializer: SerializeFunc | None = None
    materializer: ProjectFunc | None = None

    def merge(self, override: WriterHooks | None) -> WriterHooks:
        """Return a copy where every non-None field of ``override`` wins."""
        if override is None:
            return self
        return WriterHooks(
            serializer=(
                override.serializer if override.serializer is not None else self.serializer
            ),
            materializer=(
                override.materializer
                if override.materializer is not None
                else self.materializer
            ),
        )

    @classmethod
    def from_mapping(
        cls,
        mapping: Mapping[str, Any] | None = None,
        *,
        serializer: SerializeFunc | None = None,
        materializer: ProjectFunc | None = None,
    ) -> WriterHooks:
        """Build writer hooks from a config mapping and/or explicit callables."""
        writer_dict = dict(mapping or {})
        return cls(
            serializer=serializer or writer_dict.get("serializer"),
            materializer=materializer or writer_dict.get("materializer"),
        )


def encode_output(
    data: dict[str, Any],
    *,
    writer: WriterHooks | None = None,
    serializer: SerializeFunc | None = None,
) -> bytes:
    """Serialize a workflow output dictionary to bytes."""
    resolved_writer = WriterHooks(serializer=serializer).merge(writer)
    if resolved_writer.serializer is None:
        msg = "encode_output requires a serializer."
        raise ValueError(msg)
    return resolved_writer.serializer(data)


def write_output(
    data: dict[str, Any],
    *,
    writer: WriterHooks | None = None,
    materializer: ProjectFunc | None = None,
    output_path: Path | None = None,
) -> None:
    """Materialize workflow results through an optional output writer."""
    resolved_writer = WriterHooks(materializer=materializer).merge(writer)
    if resolved_writer.materializer is None:
        return
    if output_path is None:
        msg = "output_path is required when a materializer is configured."
        raise ValueError(msg)
    resolved_writer.materializer(data=data, output_path=output_path)
