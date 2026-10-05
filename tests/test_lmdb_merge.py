import tempfile
import unittest
from pathlib import Path

from datacooker.lmdb import merge_lmdb_shards
from tests.test_lmdb_rekey import _read_lmdb, _write_lmdb


class MergeLmdbTests(unittest.TestCase):
    def test_overwrite_removes_stale_keys_and_sidecars(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, dest = root / "source", root / "dest"
            _write_lmdb(source, {b"a": b"fresh"})
            _write_lmdb(dest, {b"z": b"stale"})
            sidecar = Path(str(dest) + ".index.tsv")
            sidecar.write_text("old")
            merge_lmdb_shards([source], dest, overwrite=True, map_size=10000000)
            self.assertEqual(_read_lmdb(dest), {b"a": b"fresh"})
            self.assertFalse(sidecar.exists())

    def test_failed_merge_preserves_previous_database(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            a, b, dest = root / "a", root / "b", root / "dest"
            _write_lmdb(a, {b"a": b"first"})
            _write_lmdb(b, {b"a": b"duplicate"})
            _write_lmdb(dest, {b"z": b"original"})
            with self.assertRaises(KeyError):
                merge_lmdb_shards([a, b], dest, overwrite=True,
                                  on_duplicate="error", map_size=10000000)
            self.assertEqual(_read_lmdb(dest), {b"z": b"original"})

    def test_missing_shard_preserves_previous_database(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source, dest = root / "source", root / "dest"
            _write_lmdb(source, {b"a": b"first"})
            _write_lmdb(dest, {b"z": b"original"})
            import lmdb
            with self.assertRaises(lmdb.Error):
                merge_lmdb_shards([source, root / "missing"], dest,
                                  overwrite=True, map_size=10000000)
            self.assertEqual(_read_lmdb(dest), {b"z": b"original"})
            self.assertEqual(_read_lmdb(source), {b"a": b"first"})
