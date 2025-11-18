from pipelines.cifmol import CIFMol


def convert_to_cifmol_dict(
    value: dict,
) -> dict[str, dict[str, CIFMol]]:
    """Convert a dictionary containing CIFMol data into a dictionary of CIFMol objects."""
    value, metadata = value["assembly_dict"], value["metadata_dict"]

    cifmol_dict: dict[str, dict[str, CIFMol]] = {}
    for cif_key, item in value.items():
        assembly_id, model_id, alt_id = cif_key.split("_")

        md = dict(metadata)
        md["assembly_id"] = assembly_id
        md["model_id"] = model_id
        md["alt_id"] = alt_id

        item = dict(item)
        item["metadata"] = md

        cifmol_dict[cif_key] = {"cifmol": CIFMol.from_dict(item)}

    return cifmol_dict
