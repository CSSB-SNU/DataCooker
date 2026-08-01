"""Executor protocol shared by the local and SLURM backends."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable


@dataclass(frozen=True, slots=True)
class TaskHandle:
    """A submitted (or would-be) unit of work. ``job_id`` is the backend's handle
    used to chain dependants (SLURM job id); ``None`` for a synchronous local run
    or a dry run."""

    name: str
    job_id: str | None = None
    depends_on: tuple[str, ...] = field(default_factory=tuple)


@runtime_checkable
class Executor(Protocol):
    """Runs the pipeline's stages. All methods return a TaskHandle whose
    ``job_id`` later stages depend on (``afterok`` semantics)."""

    def run_tier_shards(
        self,
        *,
        name: str,
        argv_for_shard,           # (shard_idx) -> list[str] engine command
        n_shards: int,
        mem_gb: int,
        cores: int,
        depends_on: tuple[str, ...] = (),
    ) -> TaskHandle:
        """Run all shards of one size tier (one array on SLURM)."""
        ...

    def run_once(
        self,
        *,
        name: str,
        argv: list[str],
        mem_gb: int,
        cores: int,
        depends_on: tuple[str, ...] = (),
    ) -> TaskHandle:
        """Run a single command (merge, index) after its dependencies."""
        ...
