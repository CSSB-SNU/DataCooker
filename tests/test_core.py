from __future__ import annotations

import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from datacooker import (
    CycleError,
    Inputs,
    InstructionOutputError,
    MissingDependencyError,
    Recipe,
    RecipeBook,
    StepExecutionError,
    UnknownTargetError,
    Variable,
    describe,
    execute,
    load_recipe,
    parse_dict,
    parse_file,
    variable,
    visualize,
)


class DataCookerCoreTests(unittest.TestCase):
    def test_step_api_supports_clearer_declarations(self) -> None:
        recipe = RecipeBook().step(
            outputs=variable("sum", int),
            instruction=lambda a, b: a + b,
            args=[variable("a", int), variable("b", int)],
        )

        result = execute(recipe, {"a": 2, "b": 5})

        self.assertEqual(result, {"sum": 7})

    def test_execute_returns_dict_for_single_target(self) -> None:
        recipe = RecipeBook().add(
            (("sum", int),),
            lambda a, b: a + b,
            {"args": (("a", int), ("b", int))},
        )

        result = execute(recipe, {"a": 2, "b": 5}, targets="sum")

        self.assertEqual(result, {"sum": 7})

    def test_partial_execution_only_runs_requested_subgraph(self) -> None:
        calls = {"unused": 0}
        recipe = RecipeBook()
        recipe.add(
            (("sum", int),),
            lambda a, b: a + b,
            {"args": (("a", int), ("b", int))},
        )
        recipe.add(
            (("unused", int),),
            lambda value: calls.__setitem__("unused", calls["unused"] + 1) or value,
            {"args": (("sum", int),)},
        )

        result = execute(recipe, {"a": 1, "b": 2}, targets="sum")

        self.assertEqual(result, {"sum": 3})
        self.assertEqual(calls["unused"], 0)

    def test_parse_dict_returns_requested_targets(self) -> None:
        recipe = RecipeBook().add(
            (("product", int),),
            lambda a, b: a * b,
            {"args": (("a", int), ("b", int))},
        )

        result = parse_dict(recipe, {"a": 3, "b": 4})

        self.assertEqual(result, {"product": 12})

    def test_parse_file_merges_loader_output_and_explicit_inputs(self) -> None:
        recipe = RecipeBook().add(
            (("sum", int),),
            lambda a, b, offset: a + b + offset,
            {"args": (("a", int), ("b", int)), "kwargs": {"offset": ("offset", int)}},
        )

        def load_func(file_path: Path) -> dict[str, int]:
            del file_path
            return {"a": 2, "b": 3}

        result = parse_file(
            recipe,
            Path("dummy.txt"),
            load_func,
            inputs={"offset": 4},
            targets="sum",
        )

        self.assertEqual(result, {"sum": 9})

    def test_recipe_file_uses_default_targets(self) -> None:
        recipe_module = textwrap.dedent(
            """
            from datacooker import RecipeBook, variable

            RECIPE = RecipeBook().step(
                outputs=variable("sum", int),
                instruction=lambda a, b: a + b,
                args=[variable("a", int), variable("b", int)],
            ).set_default_targets(["sum"])
            """
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            recipe_path = Path(tmpdir) / "recipe_mod.py"
            recipe_path.write_text(recipe_module)
            result = parse_dict(recipe_path, {"a": 4, "b": 6})

        self.assertEqual(result, {"sum": 10})

    def test_load_recipe_returns_default_targets(self) -> None:
        recipe_module = textwrap.dedent(
            """
            from datacooker import RecipeBook

            RECIPE = RecipeBook().add(
                (("sum", int),),
                lambda a, b: a + b,
                {"args": (("a", int), ("b", int))},
            )
            TARGETS = ["sum"]
            """
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            recipe_path = Path(tmpdir) / "recipe_mod.py"
            recipe_path.write_text(recipe_module)
            loaded_recipe, default_targets = load_recipe(recipe_path)

        self.assertIsInstance(loaded_recipe, RecipeBook)
        self.assertEqual(default_targets, ["sum"])

    def test_missing_input_raises_validation_error(self) -> None:
        recipe = RecipeBook().add(
            (("sum", int),),
            lambda a, b: a + b,
            {"args": (("a", int), ("b", int))},
        )

        with self.assertRaises(MissingDependencyError):
            execute(recipe, {"a": 1})

    def test_optional_raw_variable_allows_missing_external_input(self) -> None:
        recipe = RecipeBook().add(
            (("summary", str),),
            lambda required, optional: f"{required}:{optional}",
            {"args": (("required", str), ("optional", str | None))},
        )

        result = execute(recipe, {"required": "value"})

        self.assertEqual(result, {"summary": "value:None"})

    def test_generic_raw_variable_annotation_is_accepted(self) -> None:
        recipe = RecipeBook().add(
            (("size", int),),
            lambda payload: len(payload),
            {"args": (("payload", dict[str, int]),)},
        )

        result = execute(recipe, {"payload": {"a": 1, "b": 2}})

        self.assertEqual(result, {"size": 2})

    def test_unknown_target_raises_error(self) -> None:
        recipe = RecipeBook().add(
            (("sum", int),),
            lambda a, b: a + b,
            {"args": (("a", int), ("b", int))},
        )

        with self.assertRaises(UnknownTargetError):
            execute(recipe, {"a": 1, "b": 2}, targets="missing")

    def test_cycle_detection_raises_error(self) -> None:
        recipe = RecipeBook()
        recipe.add(
            (("a", int),),
            lambda b: b,
            {"args": (("b", int),)},
        )
        recipe.add(
            (("b", int),),
            lambda a: a,
            {"args": (("a", int),)},
        )

        with self.assertRaises(CycleError):
            execute(recipe, {})

    def test_wildcard_arguments_collect_matching_values(self) -> None:
        recipe = RecipeBook().add(
            (("items", list),),
            lambda *values: list(values),
            {"args": (("item_*", int),)},
        )

        result = execute(recipe, {"item_b": 2, "item_a": 1})

        self.assertEqual(result, {"items": [1, 2]})

    def test_recipebook_default_targets_apply_to_execute(self) -> None:
        recipe = (
            RecipeBook()
            .step(
                outputs=variable("sum", int),
                instruction=lambda a, b: a + b,
                args=[variable("a", int), variable("b", int)],
            )
            .step(
                outputs=variable("label", str),
                instruction=lambda total: f"sum={total}",
                args=[variable("sum", int)],
            )
            .set_default_targets("label")
        )

        result = execute(recipe, {"a": 1, "b": 4})

        self.assertEqual(result, {"label": "sum=5"})

    def test_describe_includes_required_and_missing_inputs(self) -> None:
        recipe = (
            RecipeBook()
            .step(
                outputs=variable("sum", int),
                instruction=lambda a, b: a + b,
                args=[variable("a", int), variable("b", int)],
            )
            .step(
                outputs=variable("label", str),
                instruction=lambda total, prefix: f"{prefix}:{total}",
                args=[variable("sum", int)],
                kwargs={"prefix": variable("prefix", str)},
            )
        )

        summary = describe(
            recipe,
            targets="label",
            available_inputs={"a", "b"},
        )

        self.assertIn("Targets: label", summary)
        self.assertIn("Required inputs: a, b, prefix", summary)
        self.assertIn("Missing inputs: prefix", summary)
        self.assertIn("1. sum <- a, b", summary)
        self.assertIn("2. label <- sum, prefix=prefix", summary)

    def test_visualize_mermaid_renders_step_graph(self) -> None:
        recipe = (
            RecipeBook()
            .step(
                outputs=variable("sum", int),
                instruction=lambda a, b: a + b,
                args=[variable("a", int), variable("b", int)],
            )
            .step(
                outputs=variable("label", str),
                instruction=lambda total: f"sum={total}",
                args=[variable("sum", int)],
            )
            .set_default_targets("label")
        )

        rendered = visualize(recipe, available_inputs={"a", "b"})

        self.assertIn("flowchart LR", rendered)
        self.assertIn('step_1{{"step 1<br/>sum <- a, b"}}', rendered)
        self.assertIn('step_2{{"step 2<br/>label <- sum"}}', rendered)
        self.assertIn("target_sum --> step_2", rendered)
        self.assertIn("class target_label requested;", rendered)

    def test_visualize_dot_renders_requested_targets(self) -> None:
        recipe = RecipeBook().step(
            outputs=variable("sum", int),
            instruction=lambda a, b: a + b,
            args=[variable("a", int), variable("b", int)],
        )

        rendered = recipe.to_dot(targets="sum", available_inputs={"a"})

        self.assertIn("digraph DataCooker {", rendered)
        self.assertIn('target_sum [label="sum", shape=doubleoctagon', rendered)
        self.assertIn('input_b [label="b", shape=ellipse, style=filled', rendered)
        self.assertIn('color="#ba3d4c"', rendered)

    def test_execution_order_is_dependency_sorted(self) -> None:
        recipe = RecipeBook()
        recipe.step(
            outputs=variable("sum", int),
            instruction=lambda a, b: a + b,
            args=[variable("a", int), variable("b", int)],
        )
        recipe.step(
            outputs=(variable("double", int), variable("triple", int)),
            instruction=lambda total: (total * 2, total * 3),
            args=[variable("sum", int)],
        )

        order = recipe.execution_order(["double", "triple"])

        self.assertEqual([step.target_names for step in order], [("sum",), ("double", "triple")])

    def test_recipe_validate_reports_required_inputs(self) -> None:
        recipe = RecipeBook()
        recipe.add(
            (("sum", int),),
            lambda a, b: a + b,
            {"args": (("a", int), ("b", int))},
        )
        recipe.add(
            (("label", str),),
            lambda total: f"sum={total}",
            {"args": (("sum", int),)},
        )

        self.assertEqual(recipe.required_inputs("label"), {"a", "b"})

    def test_public_declaration_types_are_exported(self) -> None:
        self.assertIsInstance(variable("value", int), Variable)
        self.assertIsInstance(Inputs(), Inputs)
        recipe = Recipe(
            targets=(variable("out", int),),
            instruction=lambda value: value,
            inputs=Inputs(args=(variable("value", int),)),
        )
        self.assertEqual(recipe.describe(), "out <- value")

    def test_missing_dependency_error_exposes_context(self) -> None:
        recipe = RecipeBook().step(
            outputs=variable("sum", int),
            instruction=lambda a, b: a + b,
            args=[variable("a", int), variable("b", int)],
        )

        with self.assertRaises(MissingDependencyError) as caught:
            execute(recipe, {"a": 1})

        error = caught.exception
        self.assertEqual(error.dependency_name, "b")
        self.assertEqual(error.available_inputs, ("a",))
        self.assertEqual(error.dependency_chain, ("sum",))

    def test_step_execution_error_wraps_original_exception(self) -> None:
        recipe = RecipeBook().step(
            outputs=variable("sum", int),
            instruction=lambda _a: 1 / 0,
            args=[variable("a", int)],
        )

        with self.assertRaises(StepExecutionError) as caught:
            execute(recipe, {"a": 1})

        error = caught.exception
        self.assertEqual(error.target_name, "sum")
        self.assertEqual(error.step_description, "sum <- a")
        self.assertEqual(error.dependency_chain, ("sum",))
        self.assertEqual(type(error.original_exception), ZeroDivisionError)

    def test_instruction_output_error_exposes_expected_shape(self) -> None:
        recipe = RecipeBook().step(
            outputs=(variable("left", int), variable("right", int)),
            instruction=lambda value: value,
            args=[variable("value", int)],
        )

        with self.assertRaises(InstructionOutputError) as caught:
            execute(recipe, {"value": 3})

        error = caught.exception
        self.assertEqual(error.produced_targets, ("left", "right"))
        self.assertEqual(error.expected_output_count, 2)
        self.assertIsNone(error.actual_output_count)

    def test_step_metadata_is_stored_on_recipe(self) -> None:
        recipebook = RecipeBook().step(
            outputs=variable("normalized", int),
            instruction=lambda value: value,
            args=[variable("raw", int)],
            name="normalize",
            namespace="prep.clean",
            tags=["core", "prep"],
            metadata={"owner": "tests"},
        )

        step = recipebook.steps[0]

        self.assertEqual(step.name, "normalize")
        self.assertEqual(step.namespace, "prep.clean")
        self.assertEqual(step.qualified_name, "prep.clean.normalize")
        self.assertEqual(step.tags, frozenset({"core", "prep"}))
        self.assertEqual(dict(step.metadata), {"owner": "tests"})
        self.assertIn("[prep.clean.normalize]", step.describe())

    def test_subset_by_tags_keeps_dependency_closure(self) -> None:
        recipe = (
            RecipeBook()
            .step(
                outputs=variable("sum", int),
                instruction=lambda a, b: a + b,
                args=[variable("a", int), variable("b", int)],
                name="sum",
                tags="core",
            )
            .step(
                outputs=variable("label", str),
                instruction=lambda total: f"sum={total}",
                args=[variable("sum", int)],
                name="label",
                tags="report",
            )
        )

        subset = recipe.subset(tags="report")
        result = execute(subset, {"a": 2, "b": 3})

        self.assertEqual(subset.default_targets, ["label"])
        self.assertEqual([step.target_names for step in subset.steps], [("sum",), ("label",)])
        self.assertEqual(result, {"label": "sum=5"})

    def test_with_namespace_prefixes_internal_targets_only(self) -> None:
        recipe = (
            RecipeBook()
            .step(
                outputs=variable("sum", int),
                instruction=lambda a, b: a + b,
                args=[variable("a", int), variable("b", int)],
                name="sum",
            )
            .step(
                outputs=variable("label", str),
                instruction=lambda total, prefix: f"{prefix}:{total}",
                args=[variable("sum", int)],
                kwargs={"prefix": variable("prefix", str)},
                name="label",
                namespace="report",
            )
            .set_default_targets("label")
        )

        namespaced = recipe.with_namespace("demo")
        result = execute(
            namespaced,
            {"a": 2, "b": 3, "prefix": "sum"},
        )

        self.assertEqual(namespaced.default_targets, ["demo.label"])
        self.assertEqual(namespaced.target_names(), ["demo.sum", "demo.label"])
        self.assertEqual(namespaced.steps[1].namespace, "demo.report")
        self.assertEqual(result, {"demo.label": "sum:5"})

    def test_execute_supports_tag_filtered_subgraphs(self) -> None:
        recipe = (
            RecipeBook()
            .step(
                outputs=variable("embedding", str),
                instruction=lambda text: text.upper(),
                args=[variable("text", str)],
                name="embed",
                tags="model",
            )
            .step(
                outputs=variable("length", int),
                instruction=lambda text: len(text),
                args=[variable("text", str)],
                name="length",
                tags="stats",
            )
        )

        result = execute(recipe, {"text": "abc"}, tags="model")

        self.assertEqual(result, {"embedding": "ABC"})


if __name__ == "__main__":
    unittest.main()
