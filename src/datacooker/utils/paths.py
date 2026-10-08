"""Filesystem helpers for recipe-driven workflows."""

from __future__ import annotations

import fnmatch
import os
from concurrent.futures import ThreadPoolExecutor
from itertools import islice
from pathlib import Path

# Concurrent directory listings for a recursive scan (DATACOOKER_SCAN_THREADS).
SCAN_THREADS = int(os.environ.get("DATACOOKER_SCAN_THREADS", "64"))


def scan_paths(
    root: Path,
    *,
    pattern: str = "*",
) -> list[Path]:
    """Scan filenames recursively, or resolve a root-relative glob containing '/'.

    A plain filename pattern may list alternatives separated by ``|``
    (``name.pdb|name.pdb.zst``): a file matches if it matches any of them. That is how a
    pattern takes two exact names without a trailing ``*``, which would also match the
    partial files a downloader leaves behind (``name.pdb1056791784``).

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
    alternatives = pattern.split("|")

    def _scan(directory: str) -> tuple[list[str], list[str]]:
        files: list[str] = []
        subdirs: list[str] = []
        with os.scandir(directory) as entries:
            for entry in entries:
                if entry.is_dir(follow_symlinks=False):
                    subdirs.append(entry.path)
                elif any(fnmatch.fnmatch(entry.name, alt) for alt in alternatives):
                    files.append(entry.path)
        return files, subdirs

    # Level by level, each directory listed on a thread: on a network filesystem a
    # listing is a round trip, and a serial walk of a 16M-directory tree (OpenFold long
    # raw/) takes ~60 h against ~1 h at 64 in flight. Sorted, so the result does not
    # depend on listing order.
    matches: list[str] = []
    frontier = [str(root)]
    with ThreadPoolExecutor(max_workers=SCAN_THREADS) as pool:
        while frontier:
            following: list[str] = []
            for files, subdirs in pool.map(_scan, frontier, chunksize=64):
                matches.extend(files)
                following.extend(subdirs)
            frontier = following
    return [Path(p) for p in sorted(matches)]
