from __future__ import annotations

import unittest
from pathlib import Path
from unittest.mock import patch

from datacooker.lmdb.extract import extract_lmdb_records


class ExtractErrorTests(unittest.TestCase):
    def test_failed_record_aborts_instead_of_publishing_partial_output(self) -> None:
        with (
            patch("datacooker.lmdb.extract.extract_lmdb_keys", return_value=["bad"]),
            patch("datacooker.lmdb.extract.iter_parallel_chunks", return_value=[[("bad", {}, ValueError("bad record"))]]),
            self.assertRaisesRegex(ValueError, "bad record"),
        ):
            extract_lmdb_records(Path("source"), Path("recipe.py"),
                                 deserializer=lambda _raw: {}, test_run=False)
