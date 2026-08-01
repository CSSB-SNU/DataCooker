"""Execution backends for the planning-first pipeline.

An executor takes the planner's tiers and actually runs them -- locally (tiers in
sequence, for small DBs / dev) or on SLURM (one array per tier, merge + index
chained via ``afterok``). The pipeline driver is backend-agnostic; it hands the
executor a plan and the per-tier engine command.
"""
from datacooker.executors.base import Executor, TaskHandle

__all__ = ["Executor", "TaskHandle"]
