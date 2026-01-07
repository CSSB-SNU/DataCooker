# Pipelines & Recipes

The repository ships several recipe books under `pipelines/recipe`. Each book bundles instructions from `pipelines/instructions` into an executable workflow.

## Data parsing
- `a3m_recipe_book.py`: parse A3M/FASTA alignments (headers + sequences) into biomol feature containers; see `pipelines/instructions/a3m_instructions.py`.
- `cif_recipe_book.py`: mmCIF structure parsing with biomol features; uses helpers in `pipelines/instructions/cif_instructions.py`.
- `ccd_recipe_book.py`: ingest CCD (chemical component dictionary) records.
- `extract_fasta.py`: extract FASTA from alignment/sequence sources.

## Packaging and metadata
- `filter_a3m_lmdb.py`: filter A3M datasets and store them in LMDB.
- `build_seq_hash_map.py`: build hash maps for sequence deduplication.
- `build_metadata.py`, `load_metadata.py`: enrich and load dataset metadata.
- `build_af3_training.py`: assemble inputs for AF3 training.

## Graph pipelines
- `graph_lmdb.py`, `graph_lmdb_from_attached.py`: convert graph-like biomolecular data into LMDB shards.
- `seq_cluster.py`: build sequence clusters.
- `train_valid_graph_split.py`: construct whole graphs, split them, extract edge lists, and compute edge statistics.

## Running a recipe book
```python
from pathlib import Path
from datacooker.core import parse
from pipelines.utils.convert import load_cif  # loader -> dict[str, Any]

results = parse(
    recipe_path=Path("pipelines/recipe/cif_recipe_book.py"),
    file_path=Path("data/example.cif"),
    load_func=load_cif,
    targets=["structure_container"],  # or None for all
)
```

## Customizing
- Swap instruction functions to change preprocessing.
- Add targets with `RecipeBook.add()` to attach more derived signals.
- Use wildcard args when a step should consume multiple cached keys (e.g., `"chain_*"`).
