"""Filesystem helpers for recipe-driven workflows."""

from __future__ import annotations

import fnmatch
import os
from concurrent.futures import ThreadPoolExecutor
from itertools import islice
from pathlib import Path


def scan_paths(
    root: Path,
    *,
    pattern: str = "*",
) -> list[Path]:
    """Scan filenames recursively, or resolve a root-relative glob containing '/'.

    A structured pattern can avoid traversing unrelated intermediate directories
    in large datasets. Plain filename patterns retain recursive matching.
    """
    if "/" in pattern:
        parts = pattern.split("/")
        # Split fixed-depth globs at their first wildcard directory. Independent
        # shards can overlap remote directory/stat latency without descending into
        # irrelevant search intermediates. Preserve glob order in the result.
        split = next((i for i, part in enumerate(parts[:-1]) if any(c in part for c in "*?[")), None)
        if split is not None and "**" not in parts:
            directories = iter(root.glob("/".join(parts[:split + 1])))
            remainder = "/".join(parts[split + 1:])

            def matching_files(directory: Path) -> list[Path]:
                return [path for path in directory.glob(remainder) if path.is_file()]

            with ThreadPoolExecutor() as pool:
                paths: list[Path] = []
                while batch := list(islice(directories, 128)):
                    for matches in pool.map(matching_files, batch):
                        paths.extend(matches)
                return paths
        return [path for path in root.glob(pattern) if path.is_file()]
    matches: list[Path] = []

    def _scan(directory: Path) -> None:
        with os.scandir(directory) as entries:
            for entry in entries:
                entry_path = Path(entry.path)
                if entry.is_dir(follow_symlinks=False):
                    _scan(entry_path)
                elif fnmatch.fnmatch(entry.name, pattern):
                    matches.append(entry_path)

    _scan(root)
    return matches
