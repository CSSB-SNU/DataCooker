"""Build-need predicates ``fn(config) -> bool``.

datacooker has no ``if`` primitive: the only conditional is an instruction that returns
``None`` (the engine then propagates the skip downstream). These are the reusable
predicates a caller's instruction consults to make that decision — e.g. StructCooker's
``build-all`` cooks a meta-recipe whose per-DB instruction returns ``None`` (skip) when
:func:`output_absent` is False, so a repeated run does only the missing work.

``config`` is the loaded engine config (paths already coerced), so a predicate can look
at the target, source, parameters, or anything else it needs.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any


def output_absent(config: dict[str, Any]) -> bool:
    """Report absent outputs; refuse to infer verified completion from existing artifacts."""
    projection = config.get("output_data_path") or config.get("output_path")
    target = projection or config.get("new_env_path") or config.get("env_path")
    if target and (Path(target).exists() or Path(f"{target}.index.tsv").exists()):
        msg = "Output existence cannot establish completion; use the verified release receipt path"
        raise RuntimeError(msg)
    return True


def always(config: dict[str, Any]) -> bool:  # noqa: ARG001
    """Return True always (a build is always needed) — the do-nothing default."""
    return True
