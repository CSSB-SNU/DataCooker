from __future__ import annotations

import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from datacooker import (
    CycleError,
    MissingDependencyError,
    RecipeBook,
    UnknownTargetError,
    execute,
    load_recipe,
    parse_dict,
    parse_file,
)


class DataCookerCoreTests(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main()
