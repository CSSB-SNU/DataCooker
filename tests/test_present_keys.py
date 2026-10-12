"""A key-list rebuild looks its keys up instead of listing the whole source."""
import lmdb

from datacooker.lmdb import core
from datacooker.lmdb.sharded import publish_shards


def _write(path, keys):
    with lmdb.open(str(path), map_size=1 << 20) as env, env.begin(write=True) as txn:
        for k in keys:
            txn.put(k.encode(), b"v")


def test_present_keys_in_plain_and_sharded_lmdb(tmp_path, monkeypatch):
    _write(tmp_path / "a", ["k1", "k2"])
    _write(tmp_path / "b", ["k3"])
    publish_shards([tmp_path / "a", tmp_path / "b"], tmp_path / "logical")
    monkeypatch.setattr(core, "extract_lmdb_keys", lambda *_: (_ for _ in ()).throw(AssertionError("listed")))
    assert core.present_lmdb_keys(tmp_path / "a", ["k2", "k9", "k1"]) == ["k2", "k1"]
    assert core.present_lmdb_keys(tmp_path / "logical", ["k3", "k1", "zz"]) == ["k3", "k1"]
    assert core.present_lmdb_keys(tmp_path / "missing", ["k1"]) == []
