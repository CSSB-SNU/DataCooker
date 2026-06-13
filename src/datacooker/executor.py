"""Workflow execution engine."""

from __future__ import annotations

import fnmatch
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, get_args, overload

from .cache import ParsingCache
from .errors import (
    CycleError,
    InstructionOutputError,
    MissingDependencyError,
    MissingTargetError,
    StepExecutionError,
    UnknownTargetError,
)
from .loading import load_recipe, normalize_targets
from .recipe import RecipeBook

if TYPE_CHECKING:
    from .recipe import Recipe


class Cooker:
    """Resolve and execute a static workflow graph."""

    @overload
    def __init__(
        self,
        parse_cache: ParsingCache,
        recipebook: RecipeBook,
        targets: Sequence[str] | str | None = None,
    ) -> None: ...

    @overload
    def __init__(
        self,
        parse_cache: ParsingCache,
        recipebook: str | Path,
        targets: Sequence[str] | str | None = None,
    ) -> None: ...

    def __init__(
        self,
        parse_cache: ParsingCache,
        recipebook: RecipeBook | str | Path,
        targets: Sequence[str] | str | None = None,
    ) -> None:
        self.parse_cache = parse_cache
        if isinstance(recipebook, RecipeBook):
            self.recipebook = recipebook
            self.targets = (
                recipebook.default_targets if targets is None else normalize_targets(targets)
            )
        else:
            self.recipebook, default_targets = load_recipe(recipebook)
            self.targets = (
                default_targets if targets is None else normalize_targets(targets)
            )

    def prep(
        self,
        data_dict: Mapping[str, Any],
        fields: Sequence[str] | None = None,
    ) -> None:
        """Populate the execution context with initial inputs."""
        if fields is None:
            fields = list(data_dict.keys())
        for field in fields:
            if field in data_dict:
                self.parse_cache.add_data(field, data_dict[field])
            else:
                msg = f"Field '{field}' not found in provided input data."
                raise ValueError(msg)

    def cook(
        self,
        targets: Sequence[str] | str | None = None,
        *,
        validate: bool = True,
    ) -> None:
        """Execute the requested subgraph in dependency order."""
        requested_targets = self._resolve_requested_targets(targets)
        if validate:
            self.recipebook.validate(
                available_inputs=set(self.parse_cache.keys()),
                targets=requested_targets,
            )

        active_targets: set[str] = set()
        resolution_stack: list[str] = []

        def resolve(target_name: str, target_type: object) -> object:
            if target_name in self.parse_cache:
                return self.parse_cache[target_name]

            if target_name not in self.recipebook:
                if type(None) in get_args(target_type):
                    return None
                raise self._missing_dependency_error(target_name, resolution_stack)

            if target_name in active_targets:
                raise self._cycle_error(target_name, resolution_stack)

            active_targets.add(target_name)
            resolution_stack.append(target_name)
            recipe = self.recipebook[target_name]

            try:
                resolved_args, resolved_kwargs = self._resolve_recipe_arguments(
                    recipe,
                    resolve,
                )
                result = self._invoke_recipe(
                    recipe,
                    target_name,
                    resolution_stack,
                    resolved_args,
                    resolved_kwargs,
                )
                return self._store_recipe_result(recipe, target_name, result)
            finally:
                resolution_stack.pop()
                active_targets.remove(target_name)

        for requested_target in requested_targets:
            resolve(requested_target, self._target_type(requested_target))

    def serve(self, targets: Sequence[str] | str | None = None) -> dict[str, Any]:
        """Return workflow outputs using a stable dictionary contract."""
        requested_targets = self._resolve_requested_targets(targets)
        results: dict[str, Any] = {}
        for target_name in requested_targets:
            if target_name not in self.parse_cache:
                msg = f"Target '{target_name}' was not produced during execution."
                raise MissingTargetError(
                    msg,
                    target_name=target_name,
                    requested_targets=tuple(requested_targets),
                    available_outputs=tuple(sorted(self.parse_cache.keys())),
                )
            results[target_name] = self.parse_cache[target_name]
        return results

    def _expand_wildcard_args(self, pattern: str) -> list[Any]:
        """Resolve wildcard arguments against currently available cache keys."""
        return [
            self.parse_cache[key]
            for key in sorted(self.parse_cache.keys())
            if fnmatch.fnmatch(key, pattern)
        ]

    def _resolve_requested_targets(
        self,
        targets: Sequence[str] | str | None = None,
    ) -> list[str]:
        requested_targets = normalize_targets(targets)
        if requested_targets is None:
            if self.targets is not None:
                requested_targets = list(self.targets)
            else:
                requested_targets = self.recipebook.target_names()

        for target_name in requested_targets:
            if target_name not in self.recipebook:
                msg = f"Target '{target_name}' is not declared in the recipe."
                raise UnknownTargetError(
                    msg,
                    target_name=target_name,
                    available_targets=tuple(self.recipebook.target_names()),
                )
        return requested_targets

    def _target_type(self, target_name: str) -> object:
        return next(
            variable.type
            for variable in self.recipebook.targets()
            if variable.name == target_name
        )

    def _missing_dependency_error(
        self,
        dependency_name: str,
        resolution_stack: Sequence[str],
    ) -> MissingDependencyError:
        message = f"Required dependency '{dependency_name}' is not available."
        return MissingDependencyError(
            message,
            dependency_name=dependency_name,
            requesting_target=resolution_stack[-1] if resolution_stack else None,
            dependency_chain=(*resolution_stack, dependency_name),
            available_inputs=tuple(sorted(self.parse_cache.keys())),
        )

    def _cycle_error(
        self,
        target_name: str,
        resolution_stack: Sequence[str],
    ) -> CycleError:
        message = f"Cycle detected while resolving '{target_name}'."
        return CycleError(message, cycle=(*resolution_stack, target_name))

    def _resolve_recipe_arguments(
        self,
        recipe: Recipe,
        resolve: Callable[[str, object], object],
    ) -> tuple[list[Any], dict[str, Any]]:
        resolved_args: list[Any] = []
        for variable in recipe.inputs.args:
            if _is_wildcard_pattern(variable.name):
                resolved_args.extend(self._expand_wildcard_args(variable.name))
            else:
                resolved_args.append(resolve(variable.name, variable.type))

        resolved_kwargs = {
            key: resolve(variable.name, variable.type)
            for key, variable in recipe.inputs.kwargs.items()
        }
        final_kwargs = {**recipe.inputs.params, **resolved_kwargs}
        return resolved_args, final_kwargs

    def _invoke_recipe(
        self,
        recipe: Recipe,
        target_name: str,
        resolution_stack: Sequence[str],
        resolved_args: Sequence[Any],
        resolved_kwargs: Mapping[str, Any],
    ) -> object:
        try:
            return recipe.instruction(*resolved_args, **resolved_kwargs)
        except Exception as exc:
            msg = (
                f"Instruction '{recipe.instruction_name}' failed while producing "
                f"target '{target_name}'."
            )
            raise StepExecutionError(
                msg,
                target_name=target_name,
                step_description=recipe.describe(),
                dependency_chain=tuple(resolution_stack),
                instruction_name=recipe.instruction_name,
                original_exception=exc,
            ) from exc

    def _store_recipe_result(
        self,
        recipe: Recipe,
        target_name: str,
        result: object,
    ) -> object:
        produced_targets = [variable.name for variable in recipe.targets]
        if len(produced_targets) == 1:
            self.parse_cache.add_data(produced_targets[0], result)
            return result

        if not isinstance(result, tuple) or len(result) != len(produced_targets):
            msg = (
                f"Instruction for targets {produced_targets} returned {type(result)}, "
                f"but a tuple of length {len(produced_targets)} was expected."
            )
            raise InstructionOutputError(
                msg,
                produced_targets=tuple(produced_targets),
                actual_value=result,
                expected_output_count=len(produced_targets),
                actual_output_count=len(result) if isinstance(result, tuple) else None,
            )

        resolved_output: object | None = None
        for produced_name, produced_value in zip(produced_targets, result, strict=True):
            self.parse_cache.add_data(produced_name, produced_value)
            if produced_name == target_name:
                resolved_output = produced_value
        return resolved_output


def _is_wildcard_pattern(name: str) -> bool:
    return any(char in name for char in "*?[")
