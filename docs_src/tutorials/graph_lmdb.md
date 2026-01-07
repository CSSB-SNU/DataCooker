# Graph LMDB

## Purpose
Build graph byte blobs per `CIFMol` with clustering context for LMDB storage.

## Inputs
- `cifmol` (`list[CIFMol] | None`): list of CIFMol entries.
- `seq_to_seq_hash_list` (`dict | None`): sequence → list of sequence hashes.
- `seq_hash_to_cluster` (`dict | None`): sequence hash → cluster identifier.

## Outputs
- `graph_bytes` (`dict[str, bytes]`): serialized graphs keyed by identifier.

## Steps
1. `extract_graph_per_cifmol` walks each `cifmol`, adds cluster context, and emits serialized graphs.

## Usage
```python
from pathlib import Path
from datacooker.core import parse

results = parse(
    recipe_path=Path("pipelines/recipe/graph_lmdb.py"),
    file_path=Path("data/cifmol.pkl"),  # loader should provide cifmol list + clustering dicts
    load_func=load_cifmol_with_clusters,  # user-provided
    targets=["graph_bytes"],
)
graphs = results["graph_bytes"]
```
