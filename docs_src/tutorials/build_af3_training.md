# AF3 Training LMDB

## Purpose
Filter CIFMol records for AF3 training, attach metadata, and drop signal peptides before packaging.

## Inputs
- `cifmol` (`CIFMol`): base CIF molecule entry.
- `seq2seqID` (`dict`): sequence → seqID mapping.
- `seqID2clusterID` (`dict`): seqID → clusterID mapping.
- `signalp_dict` (`dict`): SignalP predictions keyed by seqID.

## Outputs
- `cifmol_filtered_by_resolution_date` (`CIFMol`): entries passing resolution/date thresholds.
- `cifmol_no_water` (`CIFMol`): water-filtered version.
- `cifmol_attached` (`CIFMolAttached`): CIFMol with metadata attached.
- `cifmol_dict` (`dict`): final dict filtered by SignalP.

## Steps
1. `filter_by_resolution_and_date` (≤9.0Å, deposited by 2021-09-30) → `cifmol_filtered_by_resolution_date`.
2. `filter_water` removes water → `cifmol_no_water`.
3. `attach_metadata` adds seq/cluster info → `cifmol_attached`.
4. `filter_signalp` drops SignalP-positive entries → `cifmol_dict`.

## Usage
```python
from pathlib import Path
from datacooker.core import parse

results = parse(
    recipe_path=Path("pipelines/recipe/build_af3_training.py"),
    file_path=Path("data/sample.cif"),  # loader must yield cifmol + metadata dicts
    load_func=lambda p: {
        "cifmol": load_cifmol(p),  # user-provided loader
        "seq2seqID": {...},
        "seqID2clusterID": {...},
        "signalp_dict": {...},
    },
    targets=["cifmol_dict"],
)
```
