# Train/Valid Graph Split

## Purpose
Build whole graphs from an edge TSV, split into train/valid components, extract edge lists, and compute statistics.

## Inputs
- `edge_tsv_path` (`Path`): TSV with graph edges.
- `ignore_nodes` (`list`): nodes to drop from the graph.
- `train_ratio` (`float`): ratio for splitting edges into train vs valid.

## Outputs
- `whole_graph`, `polymer_graph` (`Graph`): constructed graphs.
- `subgraphs` (`dict`), `edge_wo_ligand_counts` (`dict`): component splits and ligand-filtered edge counts.
- `train_edges`, `valid_edges` (`list`): split edge identifiers.
- `train_edge_list`, `valid_edge_list` (`list`): extracted edge rows.
- `train_edge_statistics`, `valid_edge_statistics` (`list`): per-split category counts.

## Steps
1. `build_whole_graph` constructs `whole_graph` and `polymer_graph` from TSV.
2. `split_graph_by_components` splits into subgraphs and counts edges without ligands.
3. `split_train_valid` divides edges by `train_ratio`.
4. `extract_edges` pulls edge rows for train/valid lists.
5. `count_category_count` computes edge statistics for each split.

## Usage
```python
from pathlib import Path
from datacooker.core import parse

results = parse(
    recipe_path=Path("pipelines/recipe/train_valid_graph_split.py"),
    file_path=Path("data/edges.tsv"),
    load_func=lambda p: {"edge_tsv_path": p, "ignore_nodes": [], "train_ratio": 0.8},
    targets=["train_edge_list", "valid_edge_list", "train_edge_statistics", "valid_edge_statistics"],
)
```
