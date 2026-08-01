from __future__ import annotations

import sys
import tempfile
import unittest
from importlib.util import find_spec
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from datacooker.lmdb import rekey_lmdb, suffix_key_map

LMDB_AVAILABLE = find_spec("lmdb") is not None


def _write_lmdb(path: Path, items: dict[bytes, bytes]) -> None:
    import lmdb

    env = lmdb.open(str(path), map_size=10_000_000)
    with env.begin(write=True) as txn:
        for key, value in items.items():
            txn.put(key, value)
    env.sync()
    env.close()


def _read_lmdb(path: Path) -> dict[bytes, bytes]:
    import lmdb

    env = lmdb.open(str(path), readonly=True, lock=False)
    try:
        with env.begin() as txn:
            return {bytes(k): bytes(v) for k, v in txn.cursor()}
    finally:
        env.close()


@unittest.skipUnless(LMDB_AVAILABLE, "lmdb is not installed")
class RekeyLmdbTests(unittest.TestCase):
    def test_suffix_rekey_copies_values_verbatim(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src, dst = root / "src.lmdb", root / "dst.lmdb"
            _write_lmdb(src, {b"MGYP01": b"\x00value-1", b"MGYP02": b"\xffvalue-2"})

            report = rekey_lmdb(
                src,
                dst,
                key_map=suffix_key_map("_1"),
                map_size=10_000_000,
            )

            self.assertEqual(report.attempted, 2)
            self.assertEqual(report.written, 2)
            self.assertEqual(report.duplicate_keys, 0)
            self.assertEqual(
                _read_lmdb(dst),
                {b"MGYP01_1": b"\x00value-1", b"MGYP02_1": b"\xffvalue-2"},
            )

    def test_identity_copy_when_key_map_is_none(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src, dst = root / "src.lmdb", root / "dst.lmdb"
            payload = {b"a": b"1", b"b": b"2", b"c": b"3"}
            _write_lmdb(src, payload)

            rekey_lmdb(src, dst, map_size=10_000_000)

            self.assertEqual(_read_lmdb(dst), payload)

    def test_key_map_returning_none_drops_entry(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src, dst = root / "src.lmdb", root / "dst.lmdb"
            _write_lmdb(src, {b"keep": b"1", b"drop": b"2"})

            report = rekey_lmdb(
                src,
                dst,
                key_map=lambda k: None if k == "drop" else k,
                map_size=10_000_000,
            )

            self.assertEqual(report.skipped_empty, 1)
            self.assertEqual(_read_lmdb(dst), {b"keep": b"1"})

    def test_duplicate_error_policy_raises(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src, dst = root / "src.lmdb", root / "dst.lmdb"
            _write_lmdb(src, {b"a": b"1", b"b": b"2"})

            with self.assertRaises(KeyError):
                rekey_lmdb(
                    src,
                    dst,
                    key_map=lambda _k: "collide",
                    on_duplicate="error",
                    map_size=10_000_000,
                )

    def test_duplicate_skip_keeps_first(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src, dst = root / "src.lmdb", root / "dst.lmdb"
            _write_lmdb(src, {b"a": b"first", b"b": b"second"})

            report = rekey_lmdb(
                src,
                dst,
                key_map=lambda _k: "collide",
                on_duplicate="skip",
                map_size=10_000_000,
            )

            self.assertEqual(report.duplicate_keys, 1)
            self.assertEqual(_read_lmdb(dst), {b"collide": b"first"})

    def test_duplicate_overwrite_keeps_last(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src, dst = root / "src.lmdb", root / "dst.lmdb"
            _write_lmdb(src, {b"a": b"first", b"b": b"second"})

            report = rekey_lmdb(
                src,
                dst,
                key_map=lambda _k: "collide",
                on_duplicate="overwrite",
                map_size=10_000_000,
            )

            self.assertEqual(report.duplicate_keys, 1)
            self.assertEqual(_read_lmdb(dst), {b"collide": b"second"})

    def test_skip_existing_resumes_without_duplicating(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src, dst = root / "src.lmdb", root / "dst.lmdb"
            _write_lmdb(src, {b"k1": b"v1", b"k2": b"v2"})
            # Pretend a previous interrupted run already wrote k1's new key.
            _write_lmdb(dst, {b"k1_1": b"v1"})

            report = rekey_lmdb(
                src,
                dst,
                key_map=suffix_key_map("_1"),
                skip_existing=True,
                map_size=10_000_000,
            )

            self.assertEqual(report.skipped_existing, 1)
            self.assertEqual(report.written, 1)
            self.assertEqual(_read_lmdb(dst), {b"k1_1": b"v1", b"k2_1": b"v2"})


if __name__ == "__main__":
    unittest.main()
