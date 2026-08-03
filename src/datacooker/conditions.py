"""Built-in conditions for the pipeline's declarative ``condition:`` gate.

A condition is a predicate ``fn(config) -> bool``: the pipeline runs only if it
returns True, else it no-ops. It is declared as a dotted path -- in the engine config
(``condition: datacooker.conditions.output_absent``) or via ``--condition`` -- exactly
like a recipe / reader / writer reference, so an ``if`` is first-class and fully
general: any importable predicate works, not a hard-coded flag.

``config`` is the loaded engine config (paths already coerced), so a predicate can
look at the target, source, parameters, or anything else it needs.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any


def output_absent(config: dict[str, Any]) -> bool:
    """Run only if this DB's output is not built yet.

    The pipeline's final stage writes ``<target>.index.tsv``, so its presence means the
    DB is fully built; returning False then skips it. This is the condition that makes a
    repeated ``build-all`` do only the missing work (idempotent re-run).
    """
    target = config.get("new_env_path") or config.get("env_path")
    return not (target and Path(f"{target}.index.tsv").exists())


def always(config: dict[str, Any]) -> bool:  # noqa: ARG001
    """Always run -- the effective default when no ``condition`` is set."""
    return True
