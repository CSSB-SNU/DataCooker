# Getting Started

## Requirements

- Python 3.10+

## Installation

Clone the repository and install it (editable for development):

```bash
git clone https://github.com/CSSB-SNU/DataCooker.git
cd DataCooker
pip install -e .
```

Or use the bundled `pixi` environment:

```bash
pixi install
pixi shell
```

When DataCooker is consumed as a submodule (for example by StructCooker), it is
pinned under `libs/datacooker` and installed as an editable dependency, so a
recursive checkout is enough:

```bash
git submodule update --init --recursive
```

## Quick smoke test

Confirm the engine runs end to end:

```python
from datacooker import RecipeBook, execute

book = RecipeBook()
book.add(
    targets=(("features", dict),),
    instruction=lambda sequence: {"length": len(sequence)},
    inputs={"args": (("sequence", str),)},
)
book.add(
    targets=(("summary", str),),
    instruction=lambda features: f"length={features['length']}",
    inputs={"args": (("features", dict),)},
)

print(execute(book, {"sequence": "MKTAYIAK"}, targets=["summary"]))
# {'summary': 'length=8'}
```

`execute` seeds the input dict, resolves the dependency graph (`summary`
depends on `features`, which depends on `sequence`), runs the steps, and
returns only the requested targets.

## Next steps

- [Write a recipe](tutorials/make_recipe_basic.md)
- [Build an LMDB](tutorials/build_lmdb.md)
- [Concepts](index.md#core-concepts)
