"""OmegaConf-backed config loading helpers."""

from __future__ import annotations

from functools import lru_cache
from importlib import import_module
from pathlib import Path
from typing import Any

from datacooker.utils.importing import resolve_object

DEFAULT_CALLABLE_KEYS = (
    "loader",
    "load_func",
    "adapter",
    "transform_func",
    "key_transform",
    "convert_func",
    "materializer",
    "project_func",
    "serializer",
    "serialize",
    "deserializer",
    "deserialize",
    "key_func",
)


def load_config(
    config_path: Path,
    *,
    callable_keys: tuple[str, ...] = DEFAULT_CALLABLE_KEYS,
) -> dict[str, Any]:
    """Load an OmegaConf YAML file and resolve known callable fields."""
    register_default_resolvers()
    omegaconf = _load_omegaconf()
    cfg = omegaconf.load(config_path)
    cfg_dict = omegaconf.to_container(cfg, resolve=True)
    if not isinstance(cfg_dict, dict):
        msg = "Configuration file must contain a dictionary at the top level."
        raise TypeError(msg)

    config: dict[str, Any] = {str(key): value for key, value in cfg_dict.items()}
    return _resolve_config_callables(config, callable_keys=callable_keys)


def register_default_resolvers() -> None:
    """Register the default path resolver used by config-driven CLIs."""
    omegaconf = _load_omegaconf()
    try:
        omegaconf.register_new_resolver("p", lambda value: Path(value))
    except ValueError:
        return


@lru_cache(maxsize=1)
def _load_omegaconf() -> Any:
    try:
        return import_module("omegaconf").OmegaConf
    except ModuleNotFoundError as exc:
        msg = (
            "datacooker.config requires the optional 'omegaconf' dependency. "
            "Install DataCooker with the 'config' or 'cli' extras to use this feature."
        )
        raise ModuleNotFoundError(msg) from exc


def _resolve_config_callables(
    value: Any,
    *,
    callable_keys: tuple[str, ...],
    current_key: str | None = None,
) -> Any:
    if isinstance(value, dict):
        return {
            str(key): _resolve_config_callables(
                item,
                callable_keys=callable_keys,
                current_key=str(key),
            )
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [
            _resolve_config_callables(
                item,
                callable_keys=callable_keys,
                current_key=current_key,
            )
            for item in value
        ]
    if (
        isinstance(value, str)
        and current_key is not None
        and (current_key in callable_keys or current_key.endswith("_func"))
    ):
        return resolve_object(value)
    return value
