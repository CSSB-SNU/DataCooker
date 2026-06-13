"""Recipe declarations and graph introspection utilities."""

from __future__ import annotations

import fnmatch
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
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
    name: str | None = None
    namespace: str | None = None
    tags: frozenset[str] = field(default_factory=frozenset)
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        """Normalize optional metadata fields into immutable step metadata."""
        object.__setattr__(self, "name", _normalize_optional_string(self.name, field_name="name"))
        object.__setattr__(
            self,
            "namespace",
            _normalize_optional_string(self.namespace, field_name="namespace"),
        )
        object.__setattr__(
            self,
            "tags",
            _normalize_string_set(self.tags, field_name="tags"),
        )
        object.__setattr__(
            self,
            "metadata",
            cast("Mapping[str, Any]", MappingProxyType(dict(self.metadata))),
        )

    @property
    def target_names(self) -> tuple[str, ...]:
        """Return the output target names declared by this step."""
        return tuple(variable.name for variable in self.targets)

    @property
    def qualified_name(self) -> str | None:
        """Return the fully-qualified logical step name when available."""
        if self.name is None:
            return self.namespace
        if self.namespace is None:
            return self.name
        return f"{self.namespace}.{self.name}"

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
        summary = f"{outputs} <- {dependency_text}"
        label = self.qualified_name
        if label is not None:
            summary = f"[{label}] {summary}"
        if self.tags:
            summary = f"{summary} [tags: {', '.join(sorted(self.tags))}]"
        return summary

    @property
    def instruction_name(self) -> str:
        """Return a readable name for the step callable."""
        return _callable_name(self.instruction)


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

    def subset(
        self,
        *,
        targets: str | Sequence[str] | None = None,
        step_names: str | Sequence[str] | None = None,
        tags: str | Sequence[str] | None = None,
        namespaces: str | Sequence[str] | None = None,
    ) -> RecipeBook:
        """Return an executable subgraph selected by targets or step metadata."""
        filters_applied = any(
            selection is not None
            for selection in (step_names, tags, namespaces)
        )
        if targets is None and not filters_applied:
            return self.copy()

        selected_targets: list[str] = []
        if targets is not None:
            selected_targets.extend(self._normalize_requested_targets(targets, allow_empty=True))

        selected_step_names = _normalize_string_selection(step_names, field_name="step_names")
        selected_tags = _normalize_string_selection(tags, field_name="tags")
        selected_namespaces = _normalize_string_selection(
            namespaces,
            field_name="namespaces",
        )

        if filters_applied:
            for recipe in self._steps:
                if _recipe_matches_filters(
                    recipe,
                    step_names=selected_step_names,
                    tags=selected_tags,
                    namespaces=selected_namespaces,
                ):
                    selected_targets.extend(recipe.target_names)

        subset_targets = _unique_names(selected_targets)
        subset_recipe = RecipeBook()
        if not subset_targets:
            return subset_recipe

        for recipe in self.execution_order(subset_targets):
            subset_recipe._register_step(recipe)
        subset_recipe._default_targets = list(subset_targets)
        return subset_recipe

    def with_namespace(self, namespace: str) -> RecipeBook:
        """Return a namespaced copy suitable for modular reuse."""
        namespace_prefix = _normalize_required_string(namespace, field_name="namespace")
        internal_targets = set(self.target_names())
        namespaced = RecipeBook()

        for recipe in self._steps:
            namespaced._register_step(
                Recipe(
                    targets=tuple(
                        variable(_namespace_target(namespace_prefix, target.name), target.type)
                        for target in recipe.targets
                    ),
                    instruction=recipe.instruction,
                    inputs=Inputs(
                        args=tuple(
                            _namespace_dependency(
                                dependency,
                                namespace_prefix,
                                internal_targets,
                            )
                            for dependency in recipe.inputs.args
                        ),
                        kwargs={
                            name: _namespace_dependency(
                                dependency,
                                namespace_prefix,
                                internal_targets,
                            )
                            for name, dependency in recipe.inputs.kwargs.items()
                        },
                        params=dict(recipe.inputs.params),
                    ),
                    name=recipe.name,
                    namespace=_join_namespace(namespace_prefix, recipe.namespace),
                    tags=recipe.tags,
                    metadata=recipe.metadata,
                )
            )

        if self._default_targets is not None:
            namespaced._default_targets = [
                _namespace_target(namespace_prefix, target_name)
                for target_name in self._default_targets
            ]
        return namespaced

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
        name: str | None = None,
        tags: str | Sequence[str] | None = None,
        namespace: str | None = None,
        metadata: Mapping[str, Any] | None = None,
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
                name=name,
                namespace=namespace,
                tags=_normalize_string_set(tags, field_name="tags"),
                metadata=dict(metadata or {}),
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
                existing_step = self._recipes_by_target[target.name].describe()
                msg = (
                    f"Target '{target.name}' is already defined in the recipe. "
                    f"Existing step: {existing_step}"
                )
                raise DuplicateTargetError(
                    msg,
                    target_name=target.name,
                    existing_step=existing_step,
                )

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
                cycle_members = (*trail, target_name)
                cycle = " -> ".join(cycle_members)
                msg = f"Cycle detected in recipe graph: {cycle}."
                raise CycleError(msg, cycle=cycle_members)
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

    def to_mermaid(
        self,
        *,
        targets: str | Sequence[str] | None = None,
        available_inputs: Sequence[str] | set[str] | None = None,
    ) -> str:
        """Render the requested workflow graph as Mermaid flowchart text."""
        graph = self._build_graph_view(targets=targets, available_inputs=available_inputs)
        lines = [
            "flowchart LR",
            "    classDef step fill:#f4f1e8,stroke:#6b5b3e,color:#2d2418;",
            "    classDef target fill:#e6f0ff,stroke:#315ea8,color:#15233f;",
            "    classDef input fill:#eef7ea,stroke:#4f7a39,color:#1d3112;",
            "    classDef missing fill:#fdebec,stroke:#ba3d4c,color:#4a1218;",
            "    classDef requested fill:#fff0c2,stroke:#b18100,color:#5e4300;",
        ]
        lines.extend(
            f'    {_node_id("input", input_name)}["{_escape_mermaid(input_name)}"]'
            for input_name in graph.external_inputs
        )

        for index, recipe in enumerate(graph.order, start=1):
            step_id = _step_node_id(index)
            label = f"step {index}\\n{recipe.describe()}"
            lines.append(f'    {step_id}{{{{"{_escape_mermaid(label)}"}}}}')
        lines.extend(
            f'    {_node_id("target", target_name)}["{_escape_mermaid(target_name)}"]'
            for target_name in graph.produced_targets
        )

        for index, recipe in enumerate(graph.order, start=1):
            step_id = _step_node_id(index)
            dependencies = (*recipe.inputs.args, *recipe.inputs.kwargs.values())
            lines.extend(
                f"    {_dependency_node_id(graph, dependency.name)} --> {step_id}"
                for dependency in dependencies
            )
            lines.extend(
                f"    {step_id} --> {_node_id('target', target.name)}"
                for target in recipe.targets
            )
        lines.extend(f"    class {_step_node_id(index)} step;" for index in range(1, len(graph.order) + 1))

        for input_name in graph.external_inputs:
            node_class = "missing" if input_name in graph.missing_inputs else "input"
            lines.append(f"    class {_node_id('input', input_name)} {node_class};")

        for target_name in graph.produced_targets:
            node_class = "requested" if target_name in graph.requested_targets else "target"
            lines.append(f"    class {_node_id('target', target_name)} {node_class};")

        return "\n".join(lines)

    def to_dot(
        self,
        *,
        targets: str | Sequence[str] | None = None,
        available_inputs: Sequence[str] | set[str] | None = None,
    ) -> str:
        """Render the requested workflow graph as Graphviz DOT text."""
        graph = self._build_graph_view(targets=targets, available_inputs=available_inputs)
        lines = [
            "digraph DataCooker {",
            '    rankdir="LR";',
            '    graph [fontname="Helvetica"];',
            '    node [fontname="Helvetica"];',
            '    edge [fontname="Helvetica"];',
        ]

        for input_name in graph.external_inputs:
            color = "#ba3d4c" if input_name in graph.missing_inputs else "#4f7a39"
            lines.append(
                f'    {_node_id("input", input_name)} '
                f'[label="{_escape_dot(input_name)}", shape=ellipse, style=filled, '
                f'fillcolor="#eef7ea", color="{color}"];'
            )

        for index, recipe in enumerate(graph.order, start=1):
            step_id = _step_node_id(index)
            label = f"step {index}\\n{recipe.describe()}"
            lines.append(
                f'    {step_id} [label="{_escape_dot(label)}", shape=box, '
                'style="rounded,filled", fillcolor="#f4f1e8", color="#6b5b3e"];'
            )

        for target_name in graph.produced_targets:
            shape = "doubleoctagon" if target_name in graph.requested_targets else "box"
            fill = "#fff0c2" if target_name in graph.requested_targets else "#e6f0ff"
            lines.append(
                f'    {_node_id("target", target_name)} '
                f'[label="{_escape_dot(target_name)}", shape={shape}, style=filled, '
                f'fillcolor="{fill}", color="#315ea8"];'
            )

        for index, recipe in enumerate(graph.order, start=1):
            step_id = _step_node_id(index)
            dependencies = (*recipe.inputs.args, *recipe.inputs.kwargs.values())
            lines.extend(
                f"    {_dependency_node_id(graph, dependency.name)} -> {step_id};"
                for dependency in dependencies
            )
            lines.extend(
                f"    {step_id} -> {_node_id('target', target.name)};"
                for target in recipe.targets
            )

        lines.append("}")
        return "\n".join(lines)

    def visualize(
        self,
        *,
        output_format: str = "mermaid",
        targets: str | Sequence[str] | None = None,
        available_inputs: Sequence[str] | set[str] | None = None,
    ) -> str:
        """Render the workflow graph in a supported visualization format."""
        if output_format == "mermaid":
            return self.to_mermaid(targets=targets, available_inputs=available_inputs)
        if output_format == "dot":
            return self.to_dot(targets=targets, available_inputs=available_inputs)
        msg = (
            f"Unsupported visualization format '{output_format}'. "
            "Expected 'mermaid' or 'dot'."
        )
        raise InvalidRecipeError(msg)

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
                raise MissingDependencyError(
                    msg,
                    dependency_name=sorted(missing_inputs)[0],
                    dependency_chain=tuple(requested_targets),
                    available_inputs=tuple(sorted(available_inputs)),
                )
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
                raise UnknownTargetError(
                    msg,
                    target_name=target_name,
                    available_targets=tuple(self.target_names()),
                )
        return normalized

    def _build_graph_view(
        self,
        *,
        targets: str | Sequence[str] | None = None,
        available_inputs: Sequence[str] | set[str] | None = None,
    ) -> _GraphView:
        requested_targets = tuple(self._normalize_requested_targets(targets))
        order = tuple(self.execution_order(requested_targets))
        produced_targets = tuple(
            target.name for recipe in order for target in recipe.targets
        )
        produced_target_set = set(produced_targets)
        external_inputs: set[str] = set()

        for recipe in order:
            for dependency in (*recipe.inputs.args, *recipe.inputs.kwargs.values()):
                if dependency.name not in produced_target_set:
                    external_inputs.add(dependency.name)

        missing_inputs = (
            frozenset(self.missing_inputs(available_inputs, requested_targets))
            if available_inputs is not None
            else frozenset()
        )
        return _GraphView(
            requested_targets=requested_targets,
            produced_targets=produced_targets,
            external_inputs=tuple(sorted(external_inputs)),
            missing_inputs=missing_inputs,
            order=order,
        )


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


def _normalize_optional_string(
    value: str | None,
    *,
    field_name: str,
) -> str | None:
    if value is None:
        return None
    return _normalize_required_string(value, field_name=field_name)


def _normalize_required_string(value: str, *, field_name: str) -> str:
    normalized = value.strip()
    if not normalized:
        msg = f"{field_name} cannot be empty."
        raise InvalidRecipeError(msg)
    return normalized


def _normalize_string_set(
    values: str | Iterable[str] | None,
    *,
    field_name: str,
) -> frozenset[str]:
    if values is None:
        return frozenset()
    normalized = _normalize_string_selection(values, field_name=field_name)
    return frozenset(normalized)


def _normalize_string_selection(
    values: str | Iterable[str] | None,
    *,
    field_name: str,
) -> tuple[str, ...]:
    if values is None:
        return ()
    raw_values = [values] if isinstance(values, str) else list(values)
    normalized = (
        _normalize_required_string(value, field_name=field_name)
        for value in raw_values
    )
    return tuple(_unique_names(normalized))


def _unique_names(values: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    unique: list[str] = []
    for value in values:
        if value not in seen:
            seen.add(value)
            unique.append(value)
    return unique


def _recipe_matches_filters(
    recipe: Recipe,
    *,
    step_names: Sequence[str],
    tags: Sequence[str],
    namespaces: Sequence[str],
) -> bool:
    if step_names and recipe.name not in step_names:
        return False
    if tags and not set(recipe.tags).intersection(tags):
        return False
    return not namespaces or any(
        _namespace_matches(recipe.namespace, namespace)
        for namespace in namespaces
    )


def _namespace_matches(recipe_namespace: str | None, selected_namespace: str) -> bool:
    if recipe_namespace is None:
        return False
    return recipe_namespace == selected_namespace or recipe_namespace.startswith(
        f"{selected_namespace}."
    )


def _join_namespace(prefix: str, suffix: str | None) -> str:
    if suffix is None:
        return prefix
    return f"{prefix}.{suffix}"


def _namespace_target(namespace: str, target_name: str) -> str:
    return f"{namespace}.{target_name}"


def _namespace_dependency(
    dependency: Variable,
    namespace: str,
    internal_targets: set[str],
) -> Variable:
    if dependency.name in internal_targets:
        return variable(_namespace_target(namespace, dependency.name), dependency.type)
    if _is_wildcard_pattern(dependency.name) and any(
        fnmatch.fnmatch(target_name, dependency.name) for target_name in internal_targets
    ):
        return variable(_namespace_target(namespace, dependency.name), dependency.type)
    return dependency


def _allows_none(annotation: object) -> bool:
    return type(None) in get_args(annotation)


def _is_wildcard_pattern(name: str) -> bool:
    return any(char in name for char in "*?[")


@dataclass(frozen=True, slots=True)
class _GraphView:
    requested_targets: tuple[str, ...]
    produced_targets: tuple[str, ...]
    external_inputs: tuple[str, ...]
    missing_inputs: frozenset[str]
    order: tuple[Recipe, ...]


def _callable_name(instruction: Callable[..., Any]) -> str:
    return getattr(instruction, "__name__", type(instruction).__name__)


def _escape_mermaid(text: str) -> str:
    return text.replace('"', '\\"')


def _escape_dot(text: str) -> str:
    return text.replace('"', '\\"')


def _node_id(kind: str, name: str) -> str:
    safe_name = "".join(char if char.isalnum() else "_" for char in name)
    return f"{kind}_{safe_name}"


def _step_node_id(index: int) -> str:
    return f"step_{index}"


def _dependency_node_id(graph: _GraphView, dependency_name: str) -> str:
    if dependency_name in graph.produced_targets:
        return _node_id("target", dependency_name)
    return _node_id("input", dependency_name)
