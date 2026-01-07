# Extract FASTA

## Purpose
Generate FASTA strings from `CIFMol` records.

## Inputs
- `cifmol` (`CIFMol | None`): CIF molecule entry to extract sequences from.

## Outputs
- `fasta` (`str`): FASTA-formatted text.

## Steps
1. `build_fasta` reads sequences from `cifmol` and emits a single FASTA string.

## Usage
```python
from pathlib import Path
from datacooker.core import parse
from pipelines.utils.convert import load_cifmol  # user-provided loader returning {"cifmol": CIFMol}

results = parse(
    recipe_path=Path("pipelines/recipe/extract_fasta.py"),
    file_path=Path("data/example.cif"),
    load_func=load_cifmol,
    targets=["fasta"],
)
fasta_text = results["fasta"]
```
