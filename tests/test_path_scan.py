from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

from datacooker.lmdb.index import measure_files
from datacooker.utils.paths import scan_paths


class PathScanTests(TestCase):
    """Verify structured discovery without changing recursive filename matching."""

    def test_structured_glob_avoids_intermediate_hits(self) -> None:
        """A fixed-depth glob selects final alignments only."""
        with TemporaryDirectory() as directory:
            root = Path(directory)
            final = root / "P" / "001" / "P001" / "P001.a3m"
            scratch = final.parent / "hhblits" / "P001.extra.a3m"
            scratch.parent.mkdir(parents=True)
            final.write_text("query")
            scratch.write_text("intermediate")
            self.assertEqual(scan_paths(root, pattern="P/*/P*/P*.a3m"), [final])
            self.assertEqual(set(scan_paths(root, pattern="P*.a3m")), {final, scratch})

    def test_parallel_sizes_preserve_order_and_surface_missing_files(self) -> None:
        """Metadata concurrency retains planner identity and fails on vanished input."""
        with TemporaryDirectory() as directory:
            root = Path(directory)
            paths = [root / "second", root / "first"]
            paths[0].write_bytes(b"abc")
            paths[1].write_bytes(b"")
            sizes = measure_files(paths)
            self.assertEqual([(s.key, s.value_bytes) for s in sizes], [(str(paths[0]), 3), (str(paths[1]), 0)])
            self.assertEqual(measure_files([]), [])
            with self.assertRaises(FileNotFoundError):
                measure_files([root / "missing"])
