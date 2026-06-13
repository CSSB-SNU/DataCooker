"""Internal shared helpers for recipe-driven LMDB workflows."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Any

from datacooker.api import parse_dict
from datacooker.protocols import TransformFunc


def load_metadata(
    *,
    metadata_recipe: Path | None,
    metadata_input: Mapping[str, Any] | None,
    extra_kwargs: Mapping[str, Any],
) -> dict[str, Any]:
    """Load optional metadata once and return it as a plain input mapping."""
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


def rebuild_entry(
    *,
    recipe: Path,
    data: dict[str, Any],
    metadata_dict: Mapping[str, Any],
    parameter_dict: Mapping[str, Any],
    transform_func: TransformFunc | None,
    split_entries: bool,
    extra_kwargs: Mapping[str, Any],
) -> dict[str, Any] | None:
    """Apply a recipe to one deserialized LMDB entry."""
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
