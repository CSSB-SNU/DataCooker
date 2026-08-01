from __future__ import annotations

import json
import sys
import tempfile
import unittest
from importlib.util import find_spec
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from datacooker.lmdb import rebuild_lmdb

LMDB_AVAILABLE = find_spec("lmdb") is not None

# One source entry holding two sub-entries, the shape split_entries expects.
SOURCE = {
    b"1abc": {"1_1_.": {"payload": "a"}, "2_1_.": {"payload": "b"}},
    b"2xyz": {"1_1_.": {"payload": "c"}},
}

RECIPE = """
from datacooker import RecipeBook

recipe = RecipeBook()


def _wrap(payload: str) -> dict:
    return {"kept": payload.upper()}


recipe.step(
    outputs=(("wrapped", dict),),
    instruction=_wrap,
    kwargs={"payload": ("payload", str)},
)

RECIPE = recipe
TARGETS = ["wrapped"]
"""


def _write_lmdb(path: Path, items: dict[bytes, dict]) -> None:
    import lmdb

    env = lmdb.open(str(path), map_size=10_000_000)
    with env.begin(write=True) as txn:
        for key, value in items.items():
            txn.put(key, json.dumps(value).encode())
    env.sync()
    env.close()


def _read_lmdb(path: Path) -> dict[str, dict]:
    import lmdb

    env = lmdb.open(str(path), readonly=True, lock=False)
    try:
        with env.begin() as txn:
            return {k.decode(): json.loads(bytes(v)) for k, v in txn.cursor()}
    finally:
        env.close()


def _deserialize(raw: bytes) -> dict:
    return json.loads(raw)


def _serialize(data: dict) -> bytes:
    return json.dumps(data).encode()


@unittest.skipUnless(LMDB_AVAILABLE, "lmdb is not installed")
class ExplodeRebuildTests(unittest.TestCase):
    def _rebuild(self, root: Path, dst_name: str, **kwargs: object) -> object:
        recipe = root / "recipe.py"
        recipe.write_text(RECIPE, encoding="utf-8")
        src = root / "src.lmdb"
        if not src.exists():
            _write_lmdb(src, SOURCE)
        return rebuild_lmdb(
            src,
            root / dst_name,
            recipe,
            deserializer=_deserialize,
            serializer=_serialize,
            split_entries=True,
            map_size=10_000_000,
            test_run=True,
            **kwargs,  # type: ignore[arg-type]
        )

    def test_without_explode_sub_entries_stay_in_one_record(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._rebuild(root, "dst.lmdb")

            self.assertEqual(
                _read_lmdb(root / "dst.lmdb"),
                {
                    "1abc": {
                        "1_1_.": {"wrapped": {"kept": "A"}},
                        "2_1_.": {"wrapped": {"kept": "B"}},
                    },
                    "2xyz": {"1_1_.": {"wrapped": {"kept": "C"}}},
                },
            )

    def test_explode_writes_one_record_per_sub_entry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            report = self._rebuild(root, "dst.lmdb", explode_entries=True)

            self.assertEqual(report.attempted, 2)  # source keys
            self.assertEqual(report.written, 3)  # exploded records
            self.assertEqual(
                _read_lmdb(root / "dst.lmdb"),
                {
                    "1abc_1_1_.": {"wrapped": {"kept": "A"}},
                    "1abc_2_1_.": {"wrapped": {"kept": "B"}},
                    "2xyz_1_1_.": {"wrapped": {"kept": "C"}},
                },
            )

    def test_explode_value_key_unwraps_the_named_target(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._rebuild(
                root,
                "dst.lmdb",
                explode_entries=True,
                explode_value_key="wrapped",
            )

            self.assertEqual(
                _read_lmdb(root / "dst.lmdb"),
                {
                    "1abc_1_1_.": {"kept": "A"},
                    "1abc_2_1_.": {"kept": "B"},
                    "2xyz_1_1_.": {"kept": "C"},
                },
            )

    def test_explode_resume_skips_already_exploded_source_keys(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self._rebuild(root, "dst.lmdb", explode_entries=True)
            report = self._rebuild(root, "dst.lmdb", explode_entries=True)

            # Both source keys are already present as exploded records, so a
            # resume must recognise them from their destination keys' prefix.
            self.assertEqual(report.skipped_existing, 2)
            self.assertEqual(report.written, 0)

    def test_explode_requires_split_entries(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            recipe = root / "recipe.py"
            recipe.write_text(RECIPE, encoding="utf-8")
            src = root / "src.lmdb"
            _write_lmdb(src, SOURCE)

            with self.assertRaises(ValueError):
                rebuild_lmdb(
                    src,
                    root / "dst.lmdb",
                    recipe,
                    deserializer=_deserialize,
                    serializer=_serialize,
                    split_entries=False,
                    explode_entries=True,
                    map_size=10_000_000,
                )


if __name__ == "__main__":
    unittest.main()
