from __future__ import annotations

import json
import tempfile
import textwrap
import unittest
from importlib.util import find_spec
from pathlib import Path

from datacooker import (
    build_lmdb,
    count_lmdb_entries,
    extract_lmdb_keys,
    read_all_lmdb_raw,
    read_lmdb,
    read_lmdb_raw,
    rebuild_lmdb,
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

            build_lmdb(
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

            self.assertEqual(sorted(extract_lmdb_keys(env_path)), ["a", "b"])
            self.assertEqual(count_lmdb_entries(env_path), 2)
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

            rebuild_lmdb(
                old_env_path=old_env,
                new_env_path=new_env,
                recipe=rebuild_recipe,
                serialize=_serialize,
                deserialize=_deserialize,
                n_jobs=1,
                test_run=False,
            )

            self.assertEqual(
                read_lmdb(new_env, "value", deserialize=_deserialize),
                {"triple": 9},
            )

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
