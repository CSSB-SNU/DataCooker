# Tutorials

This section shows how to use DataCooker in practice, from a tiny recipe to a packaged pipeline.

## 1. Minimal custom recipe
Define a two-step recipe, load it into a cache-backed cooker, and serve results:
```python
from datacooker import ParsingCache, RecipeBook, Cooker

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
```

## 2. Run a packaged pipeline
The repository ships reusable recipe books under `pipelines/recipe`. To run one, pick a recipe file and call the helper `parse` API:
```python
from pathlib import Path
from datacooker.core import parse
from pipelines.utils.convert import load_cif  # choose a loader for your data

recipe_path = Path("pipelines/recipe/train_valid_graph_split.py")
results = parse(
    recipe_path=recipe_path,
    file_path=Path("data/edges.tsv"),
    load_func=load_cif,
    transform_func=None,
    targets=["train_edge_list", "valid_edge_list"],
)
print(results)
```

Common recipe books to explore:
- `pipelines/recipe/a3m_recipe_book.py`: A3M/FASTA alignments to biomolecular feature containers
- `pipelines/recipe/cif_recipe_book.py`, `pipelines/recipe/ccd_recipe_book.py`: mmCIF/CCD extraction
- `pipelines/recipe/graph_lmdb.py`, `pipelines/recipe/graph_lmdb_from_attached.py`: LMDB builders for graph data
- `pipelines/recipe/train_valid_graph_split.py`: graph splits with basic stats

## 2.5 Recipe library (what's included)
- **Data parsing**: `a3m_recipe_book.py`, `cif_recipe_book.py`, `ccd_recipe_book.py`, `extract_fasta.py`
- **Packaging/metadata**: `filter_a3m_lmdb.py`, `build_seq_hash_map.py`, `build_metadata.py`, `load_metadata.py`, `build_af3_training.py`
- **Graph pipelines**: `graph_lmdb.py`, `graph_lmdb_from_attached.py`, `seq_cluster.py`, `train_valid_graph_split.py`

## 3. Iterate on your own recipe
Use the same pattern as the minimal example: declare targets and inputs, add the instruction, and call `cooker.prep/cook/serve`. Keep recipe files alongside your data project or inside `pipelines/recipe/` to reuse shared instructions and transforms.

Next steps:
- Dive into the concepts to understand dependency resolution (`RecipeBook`, `Cooker`, `ParsingCache`).
- Browse `pipelines/instructions/` for reusable instruction functions.
