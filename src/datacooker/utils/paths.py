"""Filesystem helpers for recipe-driven workflows."""

from __future__ import annotations

import fnmatch
import os
from pathlib import Path


def scan_paths(
    root: Path,
    *,
    pattern: str = "*",
) -> list[Path]:
    """Recursively scan ``root`` for files matching a shell-style pattern."""
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
