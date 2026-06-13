from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from typing import Any, TypeAlias, get_args

from .errors import (
    CycleError,
    DuplicateTargetError,
    MissingDependencyError,
    UnknownTargetError,
)


@dataclass(frozen=True)
class Variable:
    """Represents a target variable in a recipe step."""

    name: str
    type: type[Any]


VariableSet: TypeAlias = tuple[Variable, ...]
VariableMap: TypeAlias = Mapping[str, Variable]

RawVariableSet: TypeAlias = tuple[tuple[str, type[Any]], ...]
RawVariableMap: TypeAlias = Mapping[str, tuple[str, type[Any]]]


@dataclass(frozen=True, slots=True)
class Inputs:
    """Represents input variables for a recipe step."""

    args: VariableSet = field(default_factory=tuple)
    kwargs: VariableMap = field(default_factory=dict)
    params: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Recipe:
    """A single step in a data processing recipe."""

    targets: VariableSet
    instruction: Callable
    inputs: Inputs


class RecipeError(KeyError):
    """Custom error for missing recipe targets."""


class RecipeBook:
    """Recipe Builder for defining data processing workflows."""

    def __init__(self) -> None:
        self.steps: list[Recipe] = []

    def _check_duplicate_targets(self, targets: VariableSet) -> None:
        for target in targets:
            if target.name in self:
                msg = f"Target '{target.name}' is already defined in the recipe."
                raise DuplicateTargetError(msg)

    def _coerce_to_variable_set(self, varset: RawVariableSet | None) -> VariableSet:
        if not varset:
            return ()
        return tuple(Variable(*t) for t in varset)

    def _coerce_to_variable_map(self, varset: RawVariableMap | None) -> VariableMap:
        if not varset:
            return {}
        return {key: Variable(*t) for key, t in varset.items()}

    def _single_add(
        self,
        targetset: RawVariableSet,
        instruction: Callable,
        inputs: dict[str, Any],
    ) -> None:
        arg_vars = self._coerce_to_variable_set(inputs.get("args"))
        kwarg_vars = self._coerce_to_variable_map(inputs.get("kwargs"))
        params_dict = inputs.get("params", {})
        final_inputs = Inputs(args=arg_vars, kwargs=kwarg_vars, params=params_dict)
        final_targetset = self._coerce_to_variable_set(targetset)

        self._check_duplicate_targets(final_targetset)
        step = Recipe(
            targets=final_targetset,
            instruction=instruction,
            inputs=final_inputs,
        )
        self.steps.append(step)

    def add(
        self,
        targets: RawVariableSet | list[RawVariableSet],
        instruction: Callable,
        inputs: dict[str, Any] | list[dict[str, Any]],
    ) -> "RecipeBook":
        """Add a new step to the recipe."""
        if isinstance(targets, list):
            if not isinstance(inputs, list):
                msg = "Targets and inputs must both be scalars or both be lists."
                raise TypeError(msg)
            if len(targets) != len(inputs):
                msg = (
                    "When providing lists of targets and inputs, "
                    "they must be of equal length."
                )
                raise ValueError(msg)
            for targetset, input_bundle in zip(targets, inputs, strict=True):
                self._single_add(targetset, instruction, input_bundle)
            return self

        if isinstance(inputs, list):
            msg = "Targets and inputs must both be scalars or both be lists."
            raise TypeError(msg)

        self._single_add(targets, instruction, inputs)

        return self

    def __contains__(self, target_name: str) -> bool:
        """Check if any step contains a target with this name."""
        return any(t.name == target_name for step in self.steps for t in step.targets)

    def __getitem__(self, target_name: str) -> Recipe:
        """Retrieve a recipe step by target name."""
        for step in self.steps:
            for t in step.targets:
                if t.name == target_name:
                    return step
        msg = f"Recipe for target '{target_name}' not found."
        raise RecipeError(msg)

    def targets(self) -> list[Variable]:
        """Return a list of all target names defined in the recipe book."""
        all_targets = []
        for step in self.steps:
            all_targets.extend(list(step.targets))
        return all_targets

    def target_names(self) -> list[str]:
        """Return the declared output target names in declaration order."""
        return [target.name for target in self.targets()]

    def dependencies(self, target_name: str) -> tuple[Variable, ...]:
        """Return the direct dependencies for a given target."""
        recipe = self[target_name]
        return (*recipe.inputs.args, *recipe.inputs.kwargs.values())

    def required_inputs(self, targets: str | list[str] | None = None) -> set[str]:
        """Return unresolved external inputs required for the requested targets."""
        requested_targets = self._normalize_requested_targets(targets)
        required: set[str] = set()
        visited: set[str] = set()

        def walk(target_name: str) -> None:
            if target_name in visited:
                return
            visited.add(target_name)
            for dependency in self.dependencies(target_name):
                if _is_wildcard_pattern(dependency.name):
                    continue
                if dependency.name in self:
                    walk(dependency.name)
                elif not _allows_none(dependency.type):
                    required.add(dependency.name)

        for target_name in requested_targets:
            walk(target_name)
        return required

    def validate(
        self,
        *,
        available_inputs: set[str] | None = None,
        targets: str | list[str] | None = None,
    ) -> None:
        """Validate the static workflow graph and its reachable dependencies."""
        requested_targets = self._normalize_requested_targets(targets)
        adjacency: dict[str, set[str]] = {name: set() for name in requested_targets}
        visited: set[str] = set()

        def build_graph(target_name: str) -> None:
            if target_name in visited:
                return
            visited.add(target_name)
            adjacency.setdefault(target_name, set())
            for dependency in self.dependencies(target_name):
                if _is_wildcard_pattern(dependency.name):
                    continue
                if dependency.name in self:
                    adjacency[target_name].add(dependency.name)
                    build_graph(dependency.name)

        for target_name in requested_targets:
            build_graph(target_name)

        if available_inputs is not None:
            missing_inputs = self.required_inputs(requested_targets) - set(available_inputs)
            if missing_inputs:
                missing = ", ".join(sorted(missing_inputs))
                msg = f"Missing required workflow inputs: {missing}."
                raise MissingDependencyError(msg)

        not_visited = 0
        visiting = 1
        visited_state = 2
        state = dict.fromkeys(adjacency, not_visited)
        trail: list[str] = []

        def visit(target_name: str) -> None:
            current_state = state[target_name]
            if current_state == visiting:
                cycle = " -> ".join([*trail, target_name])
                msg = f"Cycle detected in recipe graph: {cycle}."
                raise CycleError(msg)
            if current_state == visited_state:
                return
            state[target_name] = visiting
            trail.append(target_name)
            for dependency_name in adjacency[target_name]:
                visit(dependency_name)
            trail.pop()
            state[target_name] = visited_state

        for target_name in requested_targets:
            visit(target_name)

    def _normalize_requested_targets(
        self,
        targets: str | list[str] | None = None,
    ) -> list[str]:
        if targets is None:
            requested_targets = self.target_names()
        elif isinstance(targets, str):
            requested_targets = [targets]
        else:
            requested_targets = list(targets)

        for target_name in requested_targets:
            if target_name not in self:
                msg = f"Target '{target_name}' is not declared in the recipe."
                raise UnknownTargetError(msg)
        return requested_targets


def _allows_none(annotation: object) -> bool:
    return type(None) in get_args(annotation)


def _is_wildcard_pattern(name: str) -> bool:
    return any(char in name for char in "*?[")
