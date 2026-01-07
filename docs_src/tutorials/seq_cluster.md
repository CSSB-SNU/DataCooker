# Sequence Clustering

## Purpose
Split FASTA entries, cluster proteins/antibodies, and write consolidated cluster mappings.

## Inputs
- `tmp_dir` (`Path`): working directory for intermediate files.
- `seq_hash_map` (`Path`): path to the sequence hash map file.
- (later steps) `fasta_path_dict` (`dict`): produced by the first step.

## Outputs
- `fasta_path_dict` (`dict`): paths to separated FASTA files.
- `protein_cluster_dict`, `protein_d_cluster_dict`, `antibody_cluster_dict` (`dict`): clustering results.
- `cluster_dict` (`dict`): consolidated cluster mapping.

## Steps
1. `separate_sequences` writes per-category FASTA files → `fasta_path_dict`.
2. `protein_cluster` clusters protein sequences → `protein_cluster_dict`, `protein_d_cluster_dict`.
3. `antibody_cluster` clusters antibody sequences → `antibody_cluster_dict`.
4. `write_cluster` merges clusters and writes final `cluster_dict`.

## Usage
```python
from pathlib import Path
from datacooker.core import parse

results = parse(
    recipe_path=Path("pipelines/recipe/seq_cluster.py"),
    file_path=Path("."),  # loader should supply tmp_dir and seq_hash_map
    load_func=lambda _: {"tmp_dir": Path("/tmp/clusters"), "seq_hash_map": Path("data/seq_hash_map.json")},
    targets=["cluster_dict"],
)
clusters = results["cluster_dict"]
```
