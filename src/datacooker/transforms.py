"""Reusable, dependency-free Transform instructions for DataCooker recipes.

Small generic building blocks usable as recipe ``instruction`` callables in any
workflow — the "Transform" counterpart to :mod:`datacooker.readers` and
:mod:`datacooker.writers`. Everything here is stdlib-only; numpy/domain-heavy
transforms belong in the downstream package that owns those dependencies.

``NDArray`` only appears in annotations and is imported under ``TYPE_CHECKING``
(with ``from __future__ import annotations``) so importing this module never
requires numpy at runtime.
"""

from __future__ import annotations

import re
from collections.abc import Callable
from typing import TYPE_CHECKING, TypeVar

if TYPE_CHECKING:
    from numpy.typing import NDArray

InputType = TypeVar("InputType", str, int, float)


def single_value_instruction(
    *,
    dtype: type[InputType],
) -> Callable[..., InputType]:
    """Return an instruction that coerces a one-element input to a single value.

    The returned function 'remembers' the dtype via closure.
    """

    def _single_value(
        data: list[InputType] | NDArray,
    ) -> InputType:
        formatted_data = [dtype(datum) for datum in data]
        if len(formatted_data) != 1:
            msg = f"Expected single value, got {len(formatted_data)}"
            raise ValueError(msg)

        return formatted_data[0]

    return _single_value


def extract_float_single(*args: str | None) -> float | None:
    """Extract a single valid float from multiple string-like inputs."""
    pattern = re.compile(r"^-?\d+(?:\.\d+)?$")

    # Build mask: True if value is a valid float string
    _list = [
        float(a[0])
        for a in args
        if isinstance(a, list) and len(a) == 1 and pattern.match(a[0]) is not None
    ]
    if len(_list) == 0:
        return None
    # get largest float value
    _list.sort(reverse=True)
    return _list[0]


def key_stack(**kwargs: list[InputType] | NDArray) -> dict[str, list[InputType] | NDArray]:
    """Stack the provided keyword inputs into a single dict keyed by argument name."""
    return kwargs


def merge_dict(
    *args: list[InputType] | NDArray,
) -> dict:
    """Merge several grouped dicts (key -> column -> values), erroring on column collisions."""
    _merged_dict = {}
    for _dict in args:
        if _dict is None:
            continue
        for row1, inner_dict in _dict.items():
            if row1 not in _merged_dict:
                _merged_dict[row1] = {}
            for col, values in inner_dict.items():
                if col in _merged_dict[row1]:
                    msg = f"Duplicate column '{col}' for row '{row1}'"
                    raise ValueError(msg)
                _merged_dict[row1][col] = values
    return _merged_dict
