# CCD Recipe

## Purpose
Parse CCD chem comp tables into a clean dictionary of components, atoms, and bonds.

## Inputs
- `_chem_comp` (`str | None`): raw `_chem_comp` table.
- `_chem_comp_atom` (`str | None`): raw `_chem_comp_atom` table.
- `_chem_comp_bond` (`str | None`): raw `_chem_comp_bond` table.

## Outputs
- `_chem_comp_dict`, `_chem_comp_atom_dict`, `_chem_comp_bond_dict` (`dict`): trimmed tables with selected columns.
- `chem_comp_dict` (`dict`): parsed chem comp records (hydrogens removed, unwrapped).

## Steps
1. `get_smaller_dict` extracts selected columns from chem comp/atom/bond tables.
2. `parse_chem_comp` joins the trimmed tables into `chem_comp_dict` (removes hydrogens, unwraps bonds).

## Usage
```python
from pathlib import Path
from datacooker.core import parse

results = parse(
    recipe_path=Path("pipelines/recipe/ccd_recipe_book.py"),
    file_path=Path("data/ccd.cif"),  # loader should yield raw CCD tables
    load_func=load_ccd_tables,       # user-provided: returns dicts for three tables
    targets=["chem_comp_dict"],
)
chem_comp = results["chem_comp_dict"]
```
