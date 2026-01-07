# Getting Started

## Installation (GitHub repo)
- Requirements: Python 3.10+.
- Clone the repository:
  ```bash
  git clone https://github.com/psk6950/DataCooker.git
  cd DataCooker
  ```
- Install (editable) for development:
  ```bash
  pip install -e .
  ```
- Or use the included `pixi` environment:
  ```bash
  pixi install
  pixi shell
  ```

## Quick smoke test
Confirm the engine runs:
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
- Follow the [Tutorials](tutorials.md) to run a custom recipe and a packaged pipeline.
- Read [Concepts](concepts.md) for dependency resolution details.
