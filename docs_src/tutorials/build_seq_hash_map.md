# Build Sequence Hash Map

## Purpose
Create a hash map from FASTA entries derived from CIFMol to support clustering/deduplication.

## Inputs
- `fasta_dict` (`dict | None`): mapping of identifiers to FASTA strings.

## Outputs
- `seq_hash_map` (`dict`): hash → FASTA mapping.

## Steps
1. `build_seq_hash_map` consumes `fasta_dict` and returns `seq_hash_map`.

## Usage
```python
from pathlib import Path
from datacooker.core import parse

results = parse(
    recipe_path=Path("pipelines/recipe/build_seq_hash_map.py"),
    file_path=Path("data/fasta.json"),  # loader should provide fasta_dict
    load_func=lambda p: {"fasta_dict": load_fasta_dict(p)},  # user-provided helper
    targets=["seq_hash_map"],
)
seq_hash_map = results["seq_hash_map"]
```
