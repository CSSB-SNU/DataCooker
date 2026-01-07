# DataCooker Documentation

DataCooker is a small recipe engine for biomolecular data. You describe *what* you want to build (targets) and *how* to compute them (instructions), and DataCooker resolves dependencies, executes steps, and keeps intermediate results in a cache. The repository also ships ready-made pipelines for A3M/FASTA parsing, mmCIF/CCD processing, LMDB packaging, and graph splits used by downstream training jobs.

## Why DataCooker
- Composable recipe DSL with typed targets and explicit inputs
- Deterministic dependency resolution with a pluggable `ParsingCache`
- Batteries included biomolecular pipelines (A3M, CIF/CCD, LMDB, graph splits)
- Small surface area: one `RecipeBook`, one `Cooker`, optional helper `parse`/`rebuild`

## Quick start (custom recipe)
```python
from datacooker import ParsingCache, RecipeBook, Cooker

# Define a tiny two-step recipe
book = RecipeBook()
book.add(
    targets=(("clean_text", str),),
    instruction=lambda raw: raw.strip().upper(),
    inputs={"args": (("raw_text", str),)},
)
book.add(
    targets=(("length", int),),
    instruction=lambda text: len(text),
    inputs={"args": (("clean_text", str),)},
)

cache = ParsingCache()
cooker = Cooker(cache, book)
cooker.prep({"raw_text": "  ACDEfg  "})
cooker.cook()
print(cooker.serve(["clean_text", "length"]))
# {'clean_text': 'ACDEFG', 'length': 6}, ['clean_text', 'length']
```

## Pick an existing pipeline
- `pipelines/recipe/a3m_recipe_book.py`: parse A3M/FASTA alignments into biomol feature containers
- `pipelines/recipe/cif_recipe_book.py`, `pipelines/recipe/ccd_recipe_book.py`: mmCIF/CCD extraction
- `pipelines/recipe/graph_lmdb.py`, `pipelines/recipe/graph_lmdb_from_attached.py`: build LMDBs for graph data
- `pipelines/recipe/train_valid_graph_split.py`: split graphs into train/valid with basic stats

Use the helper API to run a pipeline:
```python
from pathlib import Path
from datacooker.core import parse
from pipelines.utils.convert import load_cif  # example loader

recipe_path = Path("pipelines/recipe/train_valid_graph_split.py")
results = parse(
    recipe_path=recipe_path,
    file_path=Path("data/edges.tsv"),
    load_func=load_cif,
    transform_func=None,
    targets=["train_edge_list", "valid_edge_list"],
)
```

## Where to go next
- [Getting started](getting-started.md): install and run a first recipe
- [Concepts & architecture](concepts.md): how recipes, targets, and the cache fit together
- [Pipelines](pipelines.md): what is packaged in `pipelines/`
- [API reference](api.md): the `datacooker` Python API surface
