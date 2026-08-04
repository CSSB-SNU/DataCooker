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
    """Return True if this DB's output is not built yet (so it still needs building).

    For an LMDB build/rebuild the final stage writes ``<target>.index.tsv``, so its
    presence means the DB is fully built. For a projection op (materialize/extract) the
    output is a plain file at ``output_data_path`` (a.k.a. ``output_path`` -- the same
    aliases the workflow CLI accepts), so its own presence is the signal; checking it
    **before** the LMDB target is deliberate, so an extract op is not fooled by its
    always-present *source* index. Returning False lets the caller skip it -- this is
    what makes a repeated ``build-all`` do only the missing work (idempotent re-run).

    A **parallel** op (``split_recipe``) writes a whole directory of files with no single
    sentinel, so it has no ``output_data_path``/target here and this returns True: it is
    always rebuilt. That is safe (never wrongly skips) and cheap in practice because those
    recipes (hmmsearch, a3m cleanup) skip work whose output already exists internally.
    """
    projection = config.get("output_data_path") or config.get("output_path")
    if projection:
        return not Path(projection).exists()
    target = config.get("new_env_path") or config.get("env_path")
    return not (target and Path(f"{target}.index.tsv").exists())


def always(config: dict[str, Any]) -> bool:  # noqa: ARG001
    """Return True always (a build is always needed) — the do-nothing default."""
    return True
