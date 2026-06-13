"""DataCooker exception hierarchy."""

from __future__ import annotations


class DataCookerError(Exception):
    """Base exception for all DataCooker errors."""


class InvalidRecipeError(DataCookerError):
    """Raised when a recipe definition is structurally invalid."""


class DuplicateTargetError(InvalidRecipeError):
    """Raised when multiple steps define the same target."""


class MissingDependencyError(InvalidRecipeError):
    """Raised when a required dependency cannot be resolved."""


class CycleError(InvalidRecipeError):
    """Raised when the workflow graph contains a cycle."""


class UnknownTargetError(DataCookerError, KeyError):
    """Raised when a requested target does not exist in the recipe."""


class MissingTargetError(DataCookerError, KeyError):
    """Raised when a requested target was not produced during execution."""
