"""Structured DataCooker exception hierarchy."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any


@dataclass(eq=False)
class DataCookerError(Exception):
    """Base exception for all DataCooker errors."""

    message: str

    def __post_init__(self) -> None:
        """Initialize the base ``Exception`` payload."""
        super().__init__(self.message)

    def __str__(self) -> str:
        """Return the rendered exception message."""
        return self.message


@dataclass(eq=False)
class InvalidRecipeError(DataCookerError):
    """Raised when a recipe definition is structurally invalid."""


@dataclass(eq=False)
class DuplicateTargetError(InvalidRecipeError):
    """Raised when multiple steps define the same target."""

    target_name: str = ""
    existing_step: str | None = None


@dataclass(eq=False)
class MissingDependencyError(InvalidRecipeError):
    """Raised when a required dependency cannot be resolved."""

    dependency_name: str = ""
    requesting_target: str | None = None
    dependency_chain: tuple[str, ...] = field(default_factory=tuple)
    available_inputs: tuple[str, ...] = field(default_factory=tuple)


@dataclass(eq=False)
class CycleError(InvalidRecipeError):
    """Raised when the workflow graph contains a cycle."""

    cycle: tuple[str, ...] = field(default_factory=tuple)


@dataclass(eq=False)
class UnknownTargetError(DataCookerError, KeyError):
    """Raised when a requested target does not exist in the recipe."""

    target_name: str = ""
    available_targets: tuple[str, ...] = field(default_factory=tuple)


@dataclass(eq=False)
class MissingTargetError(DataCookerError, KeyError):
    """Raised when a requested target was not produced during execution."""

    target_name: str = ""
    requested_targets: tuple[str, ...] = field(default_factory=tuple)
    available_outputs: tuple[str, ...] = field(default_factory=tuple)


@dataclass(eq=False)
class InstructionOutputError(DataCookerError):
    """Raised when a step returns an output shape that does not match its targets."""

    produced_targets: tuple[str, ...] = field(default_factory=tuple)
    actual_value: Any = None
    expected_output_count: int = 0
    actual_output_count: int | None = None


@dataclass(eq=False)
class StepExecutionError(DataCookerError):
    """Raised when a workflow instruction raises an exception."""

    target_name: str = ""
    step_description: str | None = None
    dependency_chain: tuple[str, ...] = field(default_factory=tuple)
    instruction_name: str | None = None
    original_exception: Exception | None = None
