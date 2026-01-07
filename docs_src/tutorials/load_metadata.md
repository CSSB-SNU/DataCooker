# Load Metadata

## Purpose
Load sequence/cluster mapping TSVs and SignalP predictions for downstream filtering and attachment.

## Inputs
- `seqID2seq_path` (`Path`): TSV mapping seqID → sequence.
- `clusterID2seqID_path` (`Path`): TSV mapping clusterID → comma-separated seqIDs.
- `signalp_dir` (`Path | None`): directory containing SignalP outputs.

## Outputs
- `seqID2seq` (`dict`), `clusterID2seqID` (`dict`): raw mappings.
- `seq2seqID` (`dict`), `seqID2clusterID` (`dict`): reversed mappings.
- `signalp_dict` (`dict`): SignalP predictions keyed by seqID.

## Steps
1. `load_tsv` reads seqID→seq and clusterID→seqID TSVs.
2. `reverse_dict` flips both mappings to get `seq2seqID` and `seqID2clusterID`.
3. `load_signalp` reads SignalP outputs from `signalp_dir`.

## Usage
```python
from pathlib import Path
from datacooker.core import parse

results = parse(
    recipe_path=Path("pipelines/recipe/load_metadata.py"),
    file_path=Path("."),  # loader can ignore and return paths below
    load_func=lambda _: {
        "seqID2seq_path": Path("data/seqID2seq.tsv"),
        "clusterID2seqID_path": Path("data/clusterID2seqID.tsv"),
        "signalp_dir": Path("data/signalp"),
    },
    targets=["seq2seqID", "seqID2clusterID", "signalp_dict"],
)
```
