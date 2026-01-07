# CIF Recipe

## Purpose
Ingest an mmCIF file, extract metadata, normalize chem comp/entity/atom tables, and build assembly/contact graph dictionaries.

## Inputs
- mmCIF category dicts: `_entry.id`, `_pdbx_database_status.recvd_initial_deposition_date`, `_refine.ls_d_res_high`, `_em_3d_reconstruction.resolution`.
- Chem comp tables: `_chem_comp`, `_chem_comp_atom`, `_chem_comp_bond`.
- Sequence/scheme tables: `_pdbx_poly_seq_scheme`, `_pdbx_nonpoly_scheme`, `_pdbx_branch_scheme`.
- Atom and entity tables: `_atom_site`, `_entity`, `_entity_poly`, `_entity_poly_seq`, `_pdbx_entity_nonpoly`, `_pdbx_entity_branch_descriptor`, `_pdbx_entity_branch_link`, `_pdbx_entity_branch_list`.
- Assembly/connection tables: `_pdbx_struct_oper_list`, `_pdbx_struct_assembly_gen`, `_struct_conn`.
- Optional: `ccd_db_path` (`Path`) to compare chem comp entries against CCD.

## Outputs
- `metadata_dict` (`dict`): id, deposition date, resolution.
- `chem_comp_dict` / `chem_comp_full_dict` (`dict`): parsed chem comp tables, optionally compared against CCD.
- `entity_dict` (`dict`): merged polymer/nonpoly/branched entities.
- `asym_dict` (`dict`): asym dictionaries with atom sites and entity attachment.
- `assembly_dict` (`dict | None`): assemblies with contact graph edges.

## Steps
1. Extract single-value metadata (deposition date, resolution) → `metadata_dict`.
2. Trim raw mmCIF tables to smaller dicts (chem comp, schemes, atom sites, entities, assemblies, struct_conn).
3. Parse chem comp and compare against CCD (if `ccd_db_path` provided).
4. Merge entity dictionaries (polymer/nonpoly/branched) and attach to asym schemes.
5. Clean atom sites/struct_conn, rearrange atoms, and build full-length asym dict.
6. Parse assembly + struct_oper, then build assemblies and extract contact graphs → `assembly_dict`.

## Usage
```python
from pathlib import Path
from datacooker.core import parse
from pipelines.utils.convert import load_cif  # loader -> dict of mmCIF category tables

results = parse(
    recipe_path=Path("pipelines/recipe/cif_recipe_book.py"),
    file_path=Path("data/example.cif"),
    load_func=load_cif,
    targets=["assembly_dict", "metadata_dict"],
)
assembly = results["assembly_dict"]
meta = results["metadata_dict"]
```
