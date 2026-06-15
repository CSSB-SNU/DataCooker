# Build an LMDB

For large datasets DataCooker drives a recipe over a directory of files and
writes one LMDB record per entry. The pipeline is config-driven through
`python -m datacooker.cli.lmdb`.

## The pieces

1. **A loader** — turns one file into seed inputs (a dict merged into the
   cache).
2. **A recipe file** — a module that exposes a module-level `RECIPE`
   (`RecipeBook`) and a `TARGETS` list naming the outputs to store.
3. **A serializer** — encodes the recipe output to bytes for LMDB.
4. **A config** — wires the above together and points at the data.

## Loader and recipe

```python
# my_pkg/readers.py
from pathlib import Path
import numpy as np

def load_entry(path: Path) -> dict:
    with np.load(path, allow_pickle=True) as handle:
        return {"arrays": {k: handle[k] for k in handle.files}}

def entry_key(path: Path) -> str:
    # every file is named the same; the id is the parent folder
    return path.parent.name
```

```python
# my_pkg/recipe.py
from datacooker import RecipeBook

recipe = RecipeBook()
recipe.add(
    targets=(("record", dict),),
    instruction=lambda arrays: dict(arrays),
    inputs={"kwargs": {"arrays": ("arrays", dict)}},
)

RECIPE = recipe
TARGETS = ["record"]
```

## Config

```yaml
# configs/build.yaml
data_dir: ${p:/data/my_dataset}
file_pattern: "entry.npz"      # matched against the file name, recursively

env_path: ${p:/data/my_dataset.lmdb}
recipe: ${p:my_pkg/recipe.py}

n_jobs: -1
chunk_size: 10000
test_run: true

key_builder: my_pkg.readers.entry_key
reader:
  loader: my_pkg.readers.load_entry
writer:
  serializer: datacooker.encode_output
```

Dotted strings under `loader`, `serializer`, `key_transform`, `key_builder`,
`adapter`, `materializer`, and `deserializer` are resolved to callables when
the config is loaded.

!!! note "Keys come from `key_builder`"
    The LMDB key defaults to the file stem (`path.name.split(".")[0]`). When
    every entry shares the same file name (`entry.npz`), set `key_builder` to a
    function that derives the id from the path — e.g. the parent folder — so
    keys do not collide.

## Run

```bash
python -m datacooker.cli.lmdb build configs/build.yaml --map-size 2000000000000

# sharded (e.g. one SLURM array task per shard)
python -m datacooker.cli.lmdb build configs/build.yaml \
    --shard-idx "${SLURM_ARRAY_TASK_ID}" --n-shards 16

# merge shards and inspect
python -m datacooker.cli.lmdb merge "/data/my_dataset.lmdb_shard*" -o /data/my_dataset.lmdb
python -m datacooker.cli.lmdb count /data/my_dataset.lmdb
```

`build` scans `data_dir` recursively (matching `file_pattern` against the file
name), derives a key per file, runs the recipe, serializes the requested
`TARGETS`, and writes the record. With `test_run: true` the first items run
serially so errors surface early. Re-running with `skip_existing` set only
processes keys that are not already present.
