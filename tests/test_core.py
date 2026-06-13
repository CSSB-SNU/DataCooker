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
    MissingDependencyError,
    Recipe,
    RecipeBook,
    UnknownTargetError,
    Variable,
    describe,
    execute,
    load_recipe,
    parse_dict,
    parse_file,
    variable,
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


if __name__ == "__main__":
    unittest.main()
