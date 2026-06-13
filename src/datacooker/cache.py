from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator, Sequence


class ParsingCache:
    """
    Store execution inputs and intermediate outputs for a workflow run.

    This cache uses a key_transform function to interpret string keys.
    By default it treats keys as flat strings, but custom transforms
    can allow nested structures (e.g. dot notation).

    """

    def __init__(
        self,
        key_transform: Callable[[str], Sequence[str]] | None = None,
    ) -> None:
        self._storage: dict[str, Any] = {}
        if key_transform is None:
            self._key_transform = lambda k: (k,)
        else:
            self._key_transform = lambda k: tuple(key_transform(k))

    def add_data(self, name: str, data: object) -> None:
        """Store Data with a given name."""
        parts = self._key_transform(name)
        cur = self._storage
        for part in parts[:-1]:
            cur = cur.setdefault(part, {})
            if not isinstance(cur, dict):
                msg = f"Cannot create nested key {name}, {part} is not a dict."
                raise TypeError(msg)
        if parts[-1] in cur:
            msg = f"Data with name '{name}' already exists in context."
            raise KeyError(msg)
        cur[parts[-1]] = data

    def __contains__(self, name: str) -> bool:
        """Return True if name exists in context."""
        parts = self._key_transform(name)
        cur: Any = self._storage
        for part in parts:
            if not isinstance(cur, dict) or part not in cur:
                return False
            cur = cur[part]
        return True

    def __getitem__(self, name: str) -> object:
        """Get data by name."""
        parts = self._key_transform(name)
        cur: Any = self._storage
        for part in parts:
            if not isinstance(cur, dict) or part not in cur:
                msg = f"Data with name '{name}' not found in context."
                raise KeyError(msg)
            cur = cur[part]
        return cur

    def get(self, name: str, default: object = None) -> object:
        """Get data by name, returning a default when the key is missing."""
        if name in self:
            return self[name]
        return default

    def keys(self) -> list[str]:
        """Return a list of all keys (flattened back to strings)."""
        result: list[str] = []

        def _collect(d: dict[str, Any], prefix: str = "") -> None:
            for k, v in d.items():
                new_key = f"{prefix}.{k}" if prefix else k
                if isinstance(v, dict):
                    _collect(v, new_key)
                else:
                    result.append(new_key)

        _collect(self._storage)
        return result

    def __iter__(self) -> Iterator[str]:
        """Iterate over flattened keys stored in the cache."""
        return iter(self.keys())


ExecutionContext = ParsingCache
