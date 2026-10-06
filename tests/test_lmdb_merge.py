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


def test_merge_map_size_never_below_the_shards(tmp_path):
    from datacooker.cli.lmdb import required_map_size

    for i, n in enumerate((1000, 3000)):
        (tmp_path / f"s{i}").mkdir()
        (tmp_path / f"s{i}" / "data.mdb").write_bytes(b"x" * n)
    assert required_map_size([tmp_path / "s0", tmp_path / "s1"]) == 5000


def test_patch_replaces_and_adds_in_single_and_sharded_dbs(tmp_path):
    import lmdb
    from datacooker.lmdb.sharded import patch_lmdb, publish_shards, open_env

    def make(path, items):
        with lmdb.open(str(path), map_size=1 << 24) as env, env.begin(write=True) as t:
            for k, v in items.items():
                t.put(k.encode(), v.encode())

    make(tmp_path / "single", {"a": "old", "b": "keep"})
    make(tmp_path / "p1", {"a": "new", "c": "added"})
    assert patch_lmdb(tmp_path / "single", tmp_path / "p1") == {"patched": 2, "replaced": 1, "added": 1, "entries": 3}
    with lmdb.open(str(tmp_path / "single"), readonly=True, lock=False) as env, env.begin() as t:
        assert t.get(b"a") == b"new" and t.get(b"b") == b"keep" and t.get(b"c") == b"added"

    make(tmp_path / "s0", {"a": "old", "b": "keep"})
    make(tmp_path / "s1", {"d": "keep"})
    publish_shards([tmp_path / "s0", tmp_path / "s1"], tmp_path / "coll")
    make(tmp_path / "p2", {"a": "new", "c": "added"})
    assert patch_lmdb(tmp_path / "coll", tmp_path / "p2")["entries"] == 4
    with open_env(tmp_path / "coll", readonly=True, lock=False) as env, env.begin() as t:
        assert t.get(b"a") == b"new" and t.get(b"c") == b"added" and t.get(b"d") == b"keep"


def test_patch_takes_several_patch_lmdbs(tmp_path):
    import lmdb
    from datacooker.lmdb.sharded import patch_lmdb

    def make(path, items):
        with lmdb.open(str(path), map_size=1 << 24) as env, env.begin(write=True) as t:
            for k, v in items.items():
                t.put(k.encode(), v.encode())

    make(tmp_path / "db", {"a": "old", "b": "keep"})
    make(tmp_path / "p0", {"a": "new"})
    make(tmp_path / "p1", {"c": "added"})
    assert patch_lmdb(tmp_path / "db", [tmp_path / "p0", tmp_path / "p1"])["entries"] == 3
