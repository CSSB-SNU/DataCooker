from __future__ import annotations

import sys
import tempfile
import textwrap
import unittest
from importlib.util import find_spec
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from datacooker import run_parallel_workflow, run_workflow

CONFIG_AVAILABLE = find_spec("omegaconf") is not None
JOBLIB_AVAILABLE = find_spec("joblib") is not None

if CONFIG_AVAILABLE:
    from datacooker.config import load_config


class DataCookerConfigTests(unittest.TestCase):
    @unittest.skipUnless(CONFIG_AVAILABLE, "omegaconf is not installed")
    def test_load_config_resolves_path_resolver_and_callables(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recipe_path = root / "example.py"
            config_path = root / "config.yaml"
            config_path.write_text(
                textwrap.dedent(
                    f"""
                    recipe: ${{p:{recipe_path}}}
                    project_func: pathlib.Path
                    transform_func: pathlib.Path
                    """
                ),
                encoding="utf-8",
            )

            config = load_config(config_path)

            self.assertEqual(config["recipe"], recipe_path)
            self.assertIs(config["project_func"], Path)
            self.assertIs(config["transform_func"], Path)

    @unittest.skipUnless(CONFIG_AVAILABLE, "omegaconf is not installed")
    def test_load_config_resolves_nested_callables(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            config_path = root / "config.yaml"
            config_path.write_text(
                textwrap.dedent(
                    """
                    steps:
                      - transform_func: pathlib.Path
                    hooks:
                      project_func: pathlib.Path
                    """
                ),
                encoding="utf-8",
            )

            config = load_config(config_path)

        self.assertIs(config["steps"][0]["transform_func"], Path)
        self.assertIs(config["hooks"]["project_func"], Path)

    def test_run_workflow_can_project_results(self) -> None:
        projected: dict[str, object] = {}

        def project_func(*, data: dict[str, object], output_path: Path) -> None:
            projected["data"] = data
            projected["output_path"] = output_path

        recipe = (
            """
            from datacooker import RecipeBook, variable

            RECIPE = RecipeBook().step(
                outputs=variable("sum", int),
                instruction=lambda a, b: a + b,
                args=[variable("a", int), variable("b", int)],
            )
            """
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recipe_path = root / "recipe_mod.py"
            recipe_path.write_text(textwrap.dedent(recipe), encoding="utf-8")
            output_path = root / "out.json"

            results = run_workflow(
                recipe_path,
                inputs={"a": 2, "b": 3},
                project_func=project_func,
                output_path=output_path,
            )

        self.assertEqual(results, {"sum": 5})
        self.assertEqual(projected["data"], {"sum": 5})
        self.assertEqual(projected["output_path"], output_path)

    @unittest.skipUnless(JOBLIB_AVAILABLE, "joblib is not installed")
    def test_run_parallel_workflow_executes_split_work_items(self) -> None:
        split_recipe = """
            from datacooker import RecipeBook, variable

            RECIPE = RecipeBook().step(
                outputs=variable("data_list", list),
                instruction=lambda base: [{"value": base}, {"value": base + 1}],
                args=[variable("base", int)],
            )
        """
        worker_recipe = """
            from datacooker import RecipeBook, variable

            RECIPE = RecipeBook().step(
                outputs=variable("double", int),
                instruction=lambda value, offset: value * 2 + offset,
                args=[variable("value", int), variable("offset", int)],
            )
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            split_recipe_path = root / "split_recipe.py"
            worker_recipe_path = root / "worker_recipe.py"
            split_recipe_path.write_text(textwrap.dedent(split_recipe), encoding="utf-8")
            worker_recipe_path.write_text(textwrap.dedent(worker_recipe), encoding="utf-8")

            report = run_parallel_workflow(
                split_recipe_path,
                worker_recipe_path,
                inputs={"base": 4, "offset": 1},
                chunk_size=1,
                n_jobs=1,
                test_run=False,
            )

        self.assertEqual(report.attempted, 2)
        self.assertEqual(report.failed, 0)
        self.assertEqual(report.outputs(), [{"double": 9}, {"double": 11}])
