"""Public protocol types used by DataCooker."""

from __future__ import annotations

from collections.abc import Sequence
from os import PathLike
from pathlib import Path
from typing import Any, Protocol


class LoadFunc(Protocol):
    """Protocol for file loaders used by workflow convenience wrappers."""

    def __call__(self, file_path: Path) -> dict[str, Any]:
        """Load a file into a workflow input dictionary."""
        ...


class TransformFunc(Protocol):
    """Protocol for key transforms used by the execution context."""

    def __call__(self, key: str) -> Sequence[str]:
        """Transform a flat key into a namespaced key path."""
        ...


class ConvertFunc(Protocol):
    """Protocol for data conversion functions."""

    def __call__(self, data: dict[str, Any]) -> dict[str, Any]:
        """Convert an input data dictionary into another data dictionary."""
        ...


class SerializeFunc(Protocol):
    """Protocol for converting workflow data into database bytes."""

    def __call__(self, data: dict[str, Any]) -> bytes:
        """Serialize a workflow output dictionary."""
        ...


class DeserializeFunc(Protocol):
    """Protocol for converting database bytes into workflow data."""

    def __call__(self, payload: bytes) -> dict[str, Any]:
        """Deserialize a workflow output dictionary."""
        ...


class ProjectFunc(Protocol):
    """Protocol for projection functions that materialize workflow results."""

    def __call__(
        self,
        *,
        data: dict[str, Any],
        output_path: PathLike[str] | str | Path,
    ) -> None:
        """Materialize workflow results to an external destination."""
        ...


class KeyFunc(Protocol):
    """Protocol for deriving database keys from filesystem paths."""

    def __call__(self, path: Path) -> str:
        """Return the database key for a given file path."""
        ...


FileReader = LoadFunc
DataAdapter = ConvertFunc
KeyTransform = TransformFunc
PayloadWriter = SerializeFunc
PayloadReader = DeserializeFunc
ResultWriter = ProjectFunc
