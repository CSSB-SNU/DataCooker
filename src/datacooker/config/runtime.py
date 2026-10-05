"""Explicit process-local reader inputs and deterministic raw-key selection."""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any


def activate_environment(config: dict[str, Any]) -> None:
    """Apply declared reader settings in a single build/planner/verification process."""
    for key, value in config.get("runtime_environment", {}).items():
        if not isinstance(key, str) or not isinstance(value, (str, Path)):
            msg = "runtime_environment must map variable names to strings or paths"
            raise TypeError(msg)
        os.environ[key] = str(value)


def canonical_inputs(paths: list[Path], key_func: Any, policy: str = "error") -> tuple[list[Path], list[dict[str, str]]]:
    """Resolve equal keys before asynchronous processing; never depend on finish order."""
    if policy not in {"error", "first_path"}:
        msg = f"Unknown duplicate_input_policy: {policy}"
        raise ValueError(msg)
    selected: dict[str, Path] = {}
    duplicates: list[dict[str, str]] = []
    for path in sorted(paths, key=str):
        key = str(key_func(path))
        if key in selected:
            if policy == "error":
                msg = f"Duplicate input key {key}: {selected[key]} and {path}"
                raise ValueError(msg)
            duplicates.append({"key": key, "selected": str(selected[key]), "excluded": str(path)})
        else:
            selected[key] = path
    return list(selected.values()), duplicates
