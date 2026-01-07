# Getting Started

## Install
- Python 3.10+ is required.
- Development install (recommended while writing recipes):
  ```bash
  pip install -e .
  ```
- If you use `pixi`, the repo already ships a lockfile:
  ```bash
  pixi install
  pixi shell
  ```

## Quick smoke test
Run a minimal recipe to confirm the engine works:
```bash
python - <<'PY'
from datacooker import ParsingCache, RecipeBook, Cooker
book = RecipeBook()
book.add(targets=(("x", int),), instruction=lambda: 1, inputs={})
cache = ParsingCache()
cooker = Cooker(cache, book)
cooker.prep({})
cooker.cook()
print("x =", cooker.serve(["x"]))
PY
```

## Project layout
- `src/datacooker`: core engine (`Cooker`, `RecipeBook`, `ParsingCache`, helpers)
- `pipelines/recipe`: reusable recipe books (A3M/FASTA, mmCIF/CCD, LMDB builders, graph splits)
- `pipelines/instructions`: atomic instruction functions used by the recipe books
- `pipelines/utils` and `pipelines/transforms`: shared helpers for conversion and transforms
- `configs`: Hydra/OmegaConf configs for pipeline runs

## Next steps
- Try an existing recipe book, e.g. `pipelines/recipe/a3m_recipe_book.py`
- Customize a recipe by swapping instructions or adding new targets
- See [Concepts](concepts.md) for how dependency resolution works
