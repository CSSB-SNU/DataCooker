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

    @classmethod
    def from_mapping(
        cls,
        mapping: Mapping[str, Any] | None = None,
        *,
        serializer: SerializeFunc | None = None,
        materializer: ProjectFunc | None = None,
        serialize: SerializeFunc | None = None,
        project_func: ProjectFunc | None = None,
    ) -> WriterHooks:
        """Build writer hooks from a config mapping plus compatibility aliases."""
        writer_dict = dict(mapping or {})
        return cls(
            serializer=(
                serializer
                or serialize
                or writer_dict.get("serializer")
                or writer_dict.get("serialize")
            ),
            materializer=(
                materializer
                or project_func
                or writer_dict.get("materializer")
                or writer_dict.get("project_func")
            ),
        )


def encode_output(
    data: dict[str, Any],
    *,
    writer: WriterHooks | None = None,
    serializer: SerializeFunc | None = None,
) -> bytes:
    """Serialize a workflow output dictionary to bytes."""
    resolved_writer = WriterHooks.from_mapping(None, serializer=serializer)
    if writer is not None:
        resolved_writer = WriterHooks(
            serializer=writer.serializer or resolved_writer.serializer,
            materializer=writer.materializer,
        )
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
    resolved_writer = WriterHooks.from_mapping(None, materializer=materializer)
    if writer is not None:
        resolved_writer = WriterHooks(
            serializer=writer.serializer,
            materializer=writer.materializer or resolved_writer.materializer,
        )
    if resolved_writer.materializer is None:
        return
    if output_path is None:
        msg = "output_path is required when a materializer is configured."
        raise ValueError(msg)
    resolved_writer.materializer(data=data, output_path=output_path)
