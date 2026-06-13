from __future__ import annotations

import json
import sys
import tempfile
import textwrap
import unittest
from importlib.util import find_spec
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from datacooker.lmdb import (
    LmdbWriteReport,
    build_lmdb,
    count_lmdb_entries,
    default_lmdb_key,
    extract_lmdb_keys,
    extract_lmdb_records,
    filter_pending_lmdb_paths,
    read_all_lmdb_raw,
    read_lmdb,
    read_lmdb_raw,
    rebuild_lmdb,
)
from datacooker.processing import parallel_process, parallel_process_report
from datacooker.utils import (
    resolve_node_config,
    resolve_object,
    scan_paths,
    shard_items,
)

LMDB_AVAILABLE = find_spec("lmdb") is not None
JOBLIB_AVAILABLE = find_spec("joblib") is not None
DB_UTILS_AVAILABLE = LMDB_AVAILABLE and JOBLIB_AVAILABLE


def _serialize(data: dict[str, object]) -> bytes:
    return json.dumps(data, sort_keys=True).encode("utf-8")


def _deserialize(payload: bytes) -> dict[str, object]:
    return json.loads(payload.decode("utf-8"))


class DataCookerUtilsTests(unittest.TestCase):
    @unittest.skipUnless(
        DB_UTILS_AVAILABLE,
        "lmdb/joblib are not installed in this environment",
    )
    def test_build_and_read_lmdb(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recipe_path = root / "recipe_mod.py"
            recipe_path.write_text(
                textwrap.dedent(
                    """
                    from datacooker import RecipeBook, variable

                    RECIPE = RecipeBook().step(
                        outputs=variable("double", int),
                        instruction=lambda value: value * 2,
                        args=[variable("value", int)],
                    ).set_default_targets("double")
                    """
                )
            )
            first = root / "a.txt"
            second = root / "b.txt"
            first.write_text("2")
            second.write_text("5")
            env_path = root / "numbers.lmdb"

            report = build_lmdb(
                first,
                second,
                env_path=env_path,
                recipe=recipe_path,
                serialize=_serialize,
                load_func=lambda file_path: {"value": int(file_path.read_text())},
                chunk_size=1,
                n_jobs=1,
                test_run=False,
            )

            self.assertIsInstance(report, LmdbWriteReport)
            self.assertEqual(report.attempted, 2)
            self.assertEqual(report.written, 2)
            self.assertEqual(report.failed, 0)
            self.assertEqual(sorted(extract_lmdb_keys(env_path)), ["a", "b"])
            self.assertEqual(default_lmdb_key(first), "a")
            self.assertEqual(count_lmdb_entries(env_path), 2)
            self.assertEqual(
                [path.name for path in filter_pending_lmdb_paths([first, second], env_path)],
                [],
            )
            self.assertEqual(read_lmdb_raw(env_path, "a"), _serialize({"double": 4}))
            self.assertEqual(
                read_all_lmdb_raw(env_path),
                {
                    "a": _serialize({"double": 4}),
                    "b": _serialize({"double": 10}),
                },
            )
            self.assertEqual(
                read_lmdb(env_path, "a", deserialize=_deserialize),
                {"double": 4},
            )
            self.assertEqual(
                extract_lmdb_records(
                    env_path,
                    recipe_path := recipe_path,
                    deserialize=_deserialize,
                    chunk_size=1,
                    n_jobs=1,
                    test_run=False,
                ),
                {
                    "a": {"double": 4},
                    "b": {"double": 10},
                },
            )

    @unittest.skipUnless(
        DB_UTILS_AVAILABLE,
        "lmdb/joblib are not installed in this environment",
    )
    def test_rebuild_lmdb_supports_flat_entries(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source_recipe = root / "source_recipe.py"
            source_recipe.write_text(
                textwrap.dedent(
                    """
                    from datacooker import RecipeBook, variable

                    RECIPE = RecipeBook().step(
                        outputs=variable("double", int),
                        instruction=lambda value: value * 2,
                        args=[variable("value", int)],
                    )
                    """
                )
            )
            rebuild_recipe = root / "rebuild_recipe.py"
            rebuild_recipe.write_text(
                textwrap.dedent(
                    """
                    from datacooker import RecipeBook, variable

                    RECIPE = RecipeBook().step(
                        outputs=variable("triple", int),
                        instruction=lambda double: double + 1,
                        args=[variable("double", int)],
                    )
                    """
                )
            )
            data_path = root / "value.txt"
            data_path.write_text("4")
            old_env = root / "old.lmdb"
            new_env = root / "new.lmdb"

            build_lmdb(
                data_path,
                env_path=old_env,
                recipe=source_recipe,
                serialize=_serialize,
                load_func=lambda file_path: {"value": int(file_path.read_text())},
                n_jobs=1,
                test_run=False,
            )

            report = rebuild_lmdb(
                old_env_path=old_env,
                new_env_path=new_env,
                recipe=rebuild_recipe,
                serialize=_serialize,
                deserialize=_deserialize,
                n_jobs=1,
                test_run=False,
            )

            self.assertEqual(report.attempted, 1)
            self.assertEqual(report.written, 1)
            self.assertEqual(report.failed, 0)
            self.assertEqual(
                read_lmdb(new_env, "value", deserialize=_deserialize),
                {"triple": 9},
            )

    @unittest.skipUnless(
        DB_UTILS_AVAILABLE,
        "lmdb/joblib are not installed in this environment",
    )
    def test_build_lmdb_can_skip_existing_entries(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recipe_path = root / "recipe_mod.py"
            recipe_path.write_text(
                textwrap.dedent(
                    """
                    from datacooker import RecipeBook, variable

                    RECIPE = RecipeBook().step(
                        outputs=variable("value", int),
                        instruction=lambda value: value,
                        args=[variable("value", int)],
                    )
                    """
                )
            )
            first = root / "a.txt"
            second = root / "b.txt"
            first.write_text("1")
            second.write_text("2")
            env_path = root / "numbers.lmdb"

            build_lmdb(
                first,
                env_path=env_path,
                recipe=recipe_path,
                serialize=_serialize,
                load_func=lambda file_path: {"value": int(file_path.read_text())},
                n_jobs=1,
                test_run=False,
            )
            report = build_lmdb(
                first,
                second,
                env_path=env_path,
                recipe=recipe_path,
                serialize=_serialize,
                load_func=lambda file_path: {"value": int(file_path.read_text())},
                skip_existing=True,
                n_jobs=1,
                test_run=False,
            )

            self.assertEqual(report.attempted, 2)
            self.assertEqual(report.written, 1)
            self.assertEqual(report.skipped_existing, 1)
            self.assertEqual(count_lmdb_entries(env_path), 2)

    def test_sharding_helpers(self) -> None:
        self.assertEqual(resolve_node_config(node_rank=1, node_count=3), (1, 3))
        self.assertEqual(shard_items(list(range(8)), node_rank=1, node_count=3), [1, 4, 7])

    def test_resolve_object(self) -> None:
        self.assertIs(resolve_object("pathlib.Path"), Path)

    def test_scan_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "a.txt").write_text("a")
            nested = root / "nested"
            nested.mkdir()
            (nested / "b.txt").write_text("b")
            (nested / "c.csv").write_text("c")

            self.assertEqual(
                sorted(path.relative_to(root).as_posix() for path in scan_paths(root, pattern="*.txt")),
                ["a.txt", "nested/b.txt"],
            )

    @unittest.skipUnless(JOBLIB_AVAILABLE, "joblib is not installed in this environment")
    def test_parallel_process(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recipe_path = root / "recipe_mod.py"
            recipe_path.write_text(
                textwrap.dedent(
                    """
                    from datacooker import RecipeBook, variable

                    RECIPE = RecipeBook().step(
                        outputs=variable("double", int),
                        instruction=lambda value, bias: value * 2 + bias,
                        args=[variable("value", int), variable("bias", int)],
                    )
                    """
                )
            )

            results = parallel_process(
                [{"value": 2}, {"value": 5}],
                inputs={"bias": 1},
                recipe=recipe_path,
                n_jobs=1,
                chunk_size=1,
                test_run=False,
            )

            self.assertEqual(results, [({"double": 5}, None), ({"double": 11}, None)])

    @unittest.skipUnless(JOBLIB_AVAILABLE, "joblib is not installed in this environment")
    def test_parallel_process_report(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            recipe_path = root / "recipe_mod.py"
            recipe_path.write_text(
                textwrap.dedent(
                    """
                    from datacooker import RecipeBook, variable

                    RECIPE = RecipeBook().step(
                        outputs=variable("double", int),
                        instruction=lambda value: value * 2,
                        args=[variable("value", int)],
                    )
                    """
                )
            )

            report = parallel_process_report(
                [{"value": 2}, {"value": 5}],
                inputs=None,
                recipe=recipe_path,
                n_jobs=1,
                chunk_size=1,
                test_run=False,
            )

            self.assertEqual(report.total, 2)
            self.assertEqual(report.assigned, 2)
            self.assertEqual(report.attempted, 2)
            self.assertEqual(report.succeeded, 2)
            self.assertEqual(report.failed, 0)
            self.assertEqual(report.outputs(), [{"double": 4}, {"double": 10}])
