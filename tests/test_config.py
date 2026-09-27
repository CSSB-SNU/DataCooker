from __future__ import annotations

import sys
import tempfile
import textwrap
import unittest
from importlib.util import find_spec
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from datacooker import (
    ReaderHooks,
    WriterHooks,
    run_lmdb_extract,
    run_recipe,
    run_recipe_batch,
)
from datacooker.lmdb import build_lmdb

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
                    writer:
                      materializer: pathlib.Path
                    reader:
                      key_transform: pathlib.Path
                    """
                ),
                encoding="utf-8",
            )

            config = load_config(config_path)

            self.assertEqual(config["recipe"], recipe_path)
            self.assertIs(config["writer"]["materializer"], Path)
            self.assertIs(config["reader"]["key_transform"], Path)

    @unittest.skipUnless(CONFIG_AVAILABLE, "omegaconf is not installed")
    def test_load_config_resolves_nested_callables(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            config_path = root / "config.yaml"
            config_path.write_text(
                textwrap.dedent(
                    """
                    steps:
                      - key_transform: pathlib.Path
                    hooks:
                      materializer: pathlib.Path
                    """
                ),
                encoding="utf-8",
            )

            config = load_config(config_path)

        self.assertIs(config["steps"][0]["key_transform"], Path)
        self.assertIs(config["hooks"]["materializer"], Path)

    def test_run_recipe_can_materialize_results(self) -> None:
        projected: dict[str, object] = {}

        def materializer(*, data: dict[str, object], output_path: Path) -> None:
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

            results = run_recipe(
                recipe_path,
                inputs={"a": 2, "b": 3},
                writer=WriterHooks(materializer=materializer),
                output_path=output_path,
            )

        self.assertEqual(results, {"sum": 5})
        self.assertEqual(projected["data"], {"sum": 5})
        self.assertEqual(projected["output_path"], output_path)

    @unittest.skipUnless(JOBLIB_AVAILABLE, "joblib is not installed")
    def test_run_recipe_batch_executes_split_work_items(self) -> None:
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

            report = run_recipe_batch(
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

    @unittest.skipUnless(JOBLIB_AVAILABLE, "joblib is not installed")
    def test_run_lmdb_extract_accepts_reader_writer_hooks(self) -> None:
        if find_spec("lmdb") is None:
            self.skipTest("lmdb is not installed")

        written: dict[str, object] = {}

        def materializer(*, data: dict[str, object], output_path: Path) -> None:
            written["data"] = data
            written["output_path"] = output_path

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source_recipe = root / "source_recipe.py"
            source_recipe.write_text(
                textwrap.dedent(
                    """
                    from datacooker import RecipeBook, variable

                    RECIPE = RecipeBook().step(
                        outputs=variable("value", int),
                        instruction=lambda value: value,
                        args=[variable("value", int)],
                    )
                    """
                ),
                encoding="utf-8",
            )
            extract_recipe = root / "extract_recipe.py"
            extract_recipe.write_text(
                textwrap.dedent(
                    """
                    from datacooker import RecipeBook, variable

                    RECIPE = RecipeBook().step(
                        outputs=variable("double", int),
                        instruction=lambda db_data: db_data["value"] * 2,
                        args=[variable("db_data", dict)],
                    )
                    """
                ),
                encoding="utf-8",
            )
            input_file = root / "value.txt"
            input_file.write_text("7", encoding="utf-8")
            env_path = root / "values.lmdb"

            build_lmdb(
                input_file,
                env_path=env_path,
                recipe=source_recipe,
                serializer=lambda data: str(data["value"]).encode("utf-8"),
                loader=lambda file_path: {"value": int(file_path.read_text())},
                n_jobs=1,
                test_run=False,
            )

            output_path = root / "stats.txt"
            result = run_lmdb_extract(
                env_path,
                extract_recipe,
                reader=ReaderHooks(deserializer=lambda payload: {"value": int(payload.decode("utf-8"))}),
                writer=WriterHooks(materializer=materializer),
                output_path=output_path,
                n_jobs=1,
                test_run=False,
            )

        self.assertEqual(result, {"value": {"double": 14}})
        self.assertEqual(written["data"], {"value": {"double": 14}})
        self.assertEqual(written["output_path"], output_path)
