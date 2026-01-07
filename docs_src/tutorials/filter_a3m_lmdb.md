# Filter A3M LMDB

## Purpose
Trim residue/chain feature containers from an A3M LMDB to a manageable depth for AF3 training.

## Inputs
- `_residue_container` (`FeatureContainer`): residue-level features.
- `_chain_container` (`FeatureContainer`): chain-level features.

## Outputs
- `residue_container` (`FeatureContainer`): filtered residue features.
- `chain_container` (`FeatureContainer`): filtered chain features.

## Steps
1. `filter_a3m` applies max MSA depth (16,384) to both containers.

## Usage
```python
from pathlib import Path
from datacooker.core import parse
from biomol.core import FeatureContainer

results = parse(
    recipe_path=Path("pipelines/recipe/filter_a3m_lmdb.py"),
    file_path=Path("data/a3m.lmdb"),  # loader should yield _residue_container/_chain_container
    load_func=lambda p: {
        "_residue_container": FeatureContainer(...),
        "_chain_container": FeatureContainer(...),
    },
    targets=["residue_container", "chain_container"],
)
```
