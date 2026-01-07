# A3M Recipe

## Purpose
Parse A3M headers/sequences and build an MSA container for downstream feature extraction.

## Inputs
- `raw_sequences` (`str | None`): raw A3M sequence block.
- `a3m_type` (`str | None`): optional tag describing the A3M source/type.
- `headers` (`list[str] | None`): header lines aligned to the sequences.

## Outputs
- `parsed_sequences` (`str`): cleaned/parsed sequences.
- `parsed_headers` (`dict`): parsed header fields keyed by sequence.
- `msa_container` (`dict`): combined sequences + headers container.

## Steps
1. `parse_sequence` → `parsed_sequences`.
2. `parse_headers` → `parsed_headers`.
3. `build_container` merges sequences and headers into `msa_container`.

## Usage
```python
from pathlib import Path
from datacooker.core import parse

results = parse(
    recipe_path=Path("pipelines/recipe/a3m_recipe_book.py"),
    file_path=Path("data/sample.a3m"),
    load_func=lambda p: {"raw_sequences": p.read_text(), "headers": []},
    targets=["msa_container"],
)
msa_container = results["msa_container"]
```
