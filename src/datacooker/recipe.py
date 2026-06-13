"""Recipe declarations and graph introspection utilities."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, TypeAlias, TypeGuard, cast, get_args

from .errors import (
    CycleError,
    DuplicateTargetError,
    InvalidRecipeError,
    MissingDependencyError,
    UnknownTargetError,
)


@dataclass(frozen=True, slots=True)
class Variable:
    """Represents a named artifact or input used by a workflow step."""

    name: str
    type: type[Any]


@dataclass(frozen=True, slots=True)
class Inputs:
    """Represents positional inputs, named inputs, and static parameters."""

    args: tuple[Variable, ...] = field(default_factory=tuple)
    kwargs: Mapping[str, Variable] = field(default_factory=dict)
    params: Mapping[str, Any] = field(default_factory=dict)


@dataclass(frozen=True, slots=True)
class Recipe:
    """Represents a single workflow step."""

    targets: tuple[Variable, ...]
    instruction: Callable[..., Any]
    inputs: Inputs

    @property
    def target_names(self) -> tuple[str, ...]:
        """Return the output target names declared by this step."""
        return tuple(variable.name for variable in self.targets)

    def describe(self) -> str:
        """Return a compact human-readable summary for the step."""
        outputs = ", ".join(self.target_names)
        positional = [variable.name for variable in self.inputs.args]
        keyword = [
            f"{name}={variable.name}" for name, variable in self.inputs.kwargs.items()
        ]
        static_params = [f"{name}={value!r}" for name, value in self.inputs.params.items()]
        inputs = positional + keyword + static_params
        dependency_text = ", ".join(inputs) if inputs else "<none>"
        return f"{outputs} <- {dependency_text}"


class RecipeError(KeyError):
    """Raised when a declared recipe target cannot be found."""


VariableSet: TypeAlias = tuple[Variable, ...]
VariableMap: TypeAlias = Mapping[str, Variable]
RawVariable: TypeAlias = tuple[str, type[Any]]
RawVariableSet: TypeAlias = tuple[RawVariable, ...]
RawVariableMap: TypeAlias = Mapping[str, RawVariable]
VariableLike: TypeAlias = Variable | RawVariable
RAW_VARIABLE_PARTS = 2


def variable(name: str, data_type: type[Any] = object) -> Variable:
    """Create a typed variable declaration."""
    return Variable(name=name, type=data_type)


class RecipeBook:
    """Declarative builder and validator for static workflow graphs."""

    def __init__(self) -> None:
        self._steps: list[Recipe] = []
        self._recipes_by_target: dict[str, Recipe] = {}
        self._default_targets: list[str] | None = None

    @property
    def steps(self) -> tuple[Recipe, ...]:
        """Return the declared steps in insertion order."""
        return tuple(self._steps)

    @property
    def default_targets(self) -> list[str] | None:
        """Return the default target list used when no explicit targets are given."""
        if self._default_targets is None:
            return None
        return list(self._default_targets)

    def copy(self) -> RecipeBook:
        """Return a shallow copy of the recipe book."""
        copied = RecipeBook()
        copied._steps = list(self._steps)
        copied._recipes_by_target = dict(self._recipes_by_target)
        copied._default_targets = self.default_targets
        return copied

    def set_default_targets(
        self,
        targets: str | Sequence[str] | None,
    ) -> RecipeBook:
        """Set or clear the default target selection for this recipe book."""
        if targets is None:
            self._default_targets = None
            return self

        normalized_targets = self._normalize_requested_targets(targets, allow_empty=True)
        self._default_targets = normalized_targets
        return self

    def step(
        self,
        outputs: VariableLike | Sequence[VariableLike],
        instruction: Callable[..., Any],
        *,
        args: VariableLike | Sequence[VariableLike] | None = None,
        kwargs: Mapping[str, VariableLike] | None = None,
        params: Mapping[str, Any] | None = None,
    ) -> RecipeBook:
        """Add a step using a clearer named-argument declaration style."""
        final_targets = _coerce_variable_sequence(outputs)
        if not final_targets:
            msg = "Recipe steps must declare at least one target."
            raise InvalidRecipeError(msg)

        step_inputs = Inputs(
            args=_coerce_optional_variable_sequence(args),
            kwargs=_coerce_variable_map(kwargs),
            params=dict(params or {}),
        )
        self._register_step(
            Recipe(
                targets=final_targets,
                instruction=instruction,
                inputs=step_inputs,
            )
        )
        return self

    def add(
        self,
        targets: RawVariableSet | list[RawVariableSet],
        instruction: Callable[..., Any],
        inputs: dict[str, Any] | list[dict[str, Any]],
    ) -> RecipeBook:
        """Backward-compatible step declaration API."""
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
                self._add_legacy_step(targetset, instruction, input_bundle)
            return self

        if isinstance(inputs, list):
            msg = "Targets and inputs must both be scalars or both be lists."
            raise TypeError(msg)

        self._add_legacy_step(targets, instruction, inputs)
        return self

    def _add_legacy_step(
        self,
        targetset: RawVariableSet,
        instruction: Callable[..., Any],
        inputs: Mapping[str, Any],
    ) -> None:
        self.step(
            outputs=targetset,
            instruction=instruction,
            args=inputs.get("args"),
            kwargs=inputs.get("kwargs"),
            params=inputs.get("params"),
        )

    def _register_step(self, recipe: Recipe) -> None:
        for target in recipe.targets:
            if target.name in self._recipes_by_target:
                msg = f"Target '{target.name}' is already defined in the recipe."
                raise DuplicateTargetError(msg)

        self._steps.append(recipe)
        for target in recipe.targets:
            self._recipes_by_target[target.name] = recipe

    def __contains__(self, target_name: str) -> bool:
        """Return whether the recipe book declares a target."""
        return target_name in self._recipes_by_target

    def __getitem__(self, target_name: str) -> Recipe:
        """Return the step that produces the requested target."""
        try:
            return self._recipes_by_target[target_name]
        except KeyError as exc:
            msg = f"Recipe for target '{target_name}' not found."
            raise RecipeError(msg) from exc

    def __len__(self) -> int:
        """Return the number of declared steps."""
        return len(self._steps)

    def targets(self) -> list[Variable]:
        """Return all declared targets in insertion order."""
        return [target for step in self._steps for target in step.targets]

    def target_names(self) -> list[str]:
        """Return all declared target names in insertion order."""
        return [target.name for target in self.targets()]

    def dependencies(self, target_name: str) -> tuple[Variable, ...]:
        """Return the direct variable dependencies for the target."""
        recipe = self[target_name]
        return (*recipe.inputs.args, *recipe.inputs.kwargs.values())

    def required_inputs(self, targets: str | Sequence[str] | None = None) -> set[str]:
        """Return unresolved external inputs needed for the requested targets."""
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

    def missing_inputs(
        self,
        available_inputs: Sequence[str] | set[str],
        targets: str | Sequence[str] | None = None,
    ) -> set[str]:
        """Return missing external inputs for the requested targets."""
        available = set(available_inputs)
        return self.required_inputs(targets) - available

    def execution_order(self, targets: str | Sequence[str] | None = None) -> list[Recipe]:
        """Return the topologically sorted steps needed for the requested targets."""
        requested_targets = self._normalize_requested_targets(targets)
        order: list[Recipe] = []
        emitted_steps: set[int] = set()
        not_visited = 0
        visiting = 1
        visited_state = 2
        state = dict.fromkeys(self.target_names(), not_visited)
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
            recipe = self[target_name]
            for dependency in self.dependencies(target_name):
                if _is_wildcard_pattern(dependency.name):
                    continue
                if dependency.name in self:
                    visit(dependency.name)
            trail.pop()
            state[target_name] = visited_state

            recipe_id = id(recipe)
            if recipe_id not in emitted_steps:
                emitted_steps.add(recipe_id)
                order.append(recipe)

        for target_name in requested_targets:
            visit(target_name)

        return order

    def describe(
        self,
        *,
        targets: str | Sequence[str] | None = None,
        available_inputs: Sequence[str] | set[str] | None = None,
    ) -> str:
        """Return a human-readable workflow summary for the requested targets."""
        requested_targets = self._normalize_requested_targets(targets)
        required_inputs = sorted(self.required_inputs(requested_targets))
        missing_inputs = (
            sorted(self.missing_inputs(available_inputs, requested_targets))
            if available_inputs is not None
            else []
        )
        lines = [
            "Targets: " + ", ".join(requested_targets),
            "Required inputs: " + (", ".join(required_inputs) if required_inputs else "<none>"),
        ]
        if available_inputs is not None:
            lines.append(
                "Missing inputs: "
                + (", ".join(missing_inputs) if missing_inputs else "<none>")
            )
        lines.append("Execution order:")
        for index, recipe in enumerate(self.execution_order(requested_targets), start=1):
            lines.append(f"{index}. {recipe.describe()}")
        return "\n".join(lines)

    def validate(
        self,
        *,
        available_inputs: set[str] | None = None,
        targets: str | Sequence[str] | None = None,
    ) -> None:
        """Validate missing inputs, unknown targets, and dependency cycles."""
        requested_targets = self._normalize_requested_targets(targets)
        if available_inputs is not None:
            missing_inputs = self.missing_inputs(available_inputs, requested_targets)
            if missing_inputs:
                missing = ", ".join(sorted(missing_inputs))
                msg = f"Missing required workflow inputs: {missing}."
                raise MissingDependencyError(msg)
        self.execution_order(requested_targets)

    def _normalize_requested_targets(
        self,
        targets: str | Sequence[str] | None,
        *,
        allow_empty: bool = False,
    ) -> list[str]:
        if targets is None:
            normalized = (
                list(self._default_targets)
                if self._default_targets is not None
                else self.target_names()
            )
        elif isinstance(targets, str):
            normalized = [targets]
        else:
            normalized = list(targets)

        if not allow_empty and not normalized:
            msg = "No targets were requested and the recipe has no default targets."
            raise InvalidRecipeError(msg)

        for target_name in normalized:
            if target_name not in self:
                msg = f"Target '{target_name}' is not declared in the recipe."
                raise UnknownTargetError(msg)
        return normalized


def _is_raw_variable(value: object) -> TypeGuard[RawVariable]:
    return (
        isinstance(value, tuple)
        and len(value) == RAW_VARIABLE_PARTS
        and isinstance(value[0], str)
        and isinstance(value[1], type)
    )


def _coerce_variable(value: VariableLike) -> Variable:
    if isinstance(value, Variable):
        return value
    if _is_raw_variable(value):
        return Variable(*value)
    msg = f"Invalid variable declaration: {value!r}."
    raise InvalidRecipeError(msg)


def _coerce_variable_sequence(
    values: VariableLike | Sequence[VariableLike],
) -> tuple[Variable, ...]:
    if isinstance(values, Variable) or _is_raw_variable(values):
        return (_coerce_variable(values),)
    sequence = cast("Sequence[VariableLike]", values)
    return tuple(_coerce_variable(value) for value in sequence)


def _coerce_optional_variable_sequence(
    values: VariableLike | Sequence[VariableLike] | None,
) -> tuple[Variable, ...]:
    if values is None:
        return ()
    return _coerce_variable_sequence(values)


def _coerce_variable_map(
    values: Mapping[str, VariableLike] | None,
) -> dict[str, Variable]:
    if not values:
        return {}
    return {name: _coerce_variable(value) for name, value in values.items()}


def _allows_none(annotation: object) -> bool:
    return type(None) in get_args(annotation)


def _is_wildcard_pattern(name: str) -> bool:
    return any(char in name for char in "*?[")
