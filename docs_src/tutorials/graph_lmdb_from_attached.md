# Graph LMDB (Attached)

## Purpose
Serialize graphs from already-attached `CIFMol` entries (metadata included) into byte blobs for LMDB.

## Inputs
- `cifmol` (`list[CIFMol] | None`): list of attached CIFMol entries (with metadata).

## Outputs
- `graph_bytes` (`dict[str, bytes]`): serialized graphs keyed by identifier.

## Steps
1. `extract_graph_per_cifmol_attached` processes each `cifmol` and emits serialized graphs.

## Usage
```python
from pathlib import Path
from datacooker.core import parse

results = parse(
    recipe_path=Path("pipelines/recipe/graph_lmdb_from_attached.py"),
    file_path=Path("data/cifmol_attached.pkl"),
    load_func=lambda p: {"cifmol": load_attached_cifmol_list(p)},  # user-provided
    targets=["graph_bytes"],
)
```
