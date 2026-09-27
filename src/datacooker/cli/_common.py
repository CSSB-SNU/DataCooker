"""Shared helpers for DataCooker's optional CLI entrypoints."""

from __future__ import annotations

from pathlib import Path
from typing import Any


class _Missing:
    pass


MISSING = _Missing()


def pop_config_value(
    config: dict[str, Any],
    *names: str,
    default: Any = MISSING,
) -> Any:
    """Pop the first present config value among several aliases."""
    for name in names:
        if name in config:
            return config.pop(name)
    if default is not MISSING:
        return default
    joined_names = ", ".join(names)
    msg = f"Missing required config value. Expected one of: {joined_names}."
    raise KeyError(msg)


def coerce_path(config: dict[str, Any], *names: str) -> Path | None:
    """Return a config value as a ``Path`` if present."""
    for name in names:
        if name not in config:
            continue
        value = config[name]
        if value is None:
            return None
        path = Path(value) if isinstance(value, str) else value
        config[name] = path
        return path
    return None


def resolve_optional_int(
    cli_value: int | None,
    config: dict[str, Any],
    candidate_keys: tuple[str, ...],
) -> int | None:
    """Resolve an optional integer from CLI first, then config aliases."""
    if cli_value is not None:
        return cli_value
    for key in candidate_keys:
        value = config.get(key)
        if value is not None:
            return int(value)
    return None


def with_shard_suffix(path: Path, shard_idx: int) -> Path:
    """Append a shard suffix to an LMDB path."""
    return path.with_name(f"{path.stem}_shard{shard_idx}{path.suffix}")
