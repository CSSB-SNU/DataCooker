from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

from .cifmol import CIFMol
from .utils import to_cif

__all__ = [
    "CIFMol",
    "to_cif",
]
