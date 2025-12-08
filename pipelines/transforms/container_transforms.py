from biomol.core import FeatureContainer


def convert_to_msa_container(
    value: dict,
) -> dict[str, FeatureContainer]:
    """Convert a dictionary containing CIFMol data into a dictionary of CIFMol objects."""
    value = value["msa_container"]
    residue_container, chain_container = value["residue_container"], value["chain_container"]
    residue_container = FeatureContainer.from_dict(residue_container)
    chain_container = FeatureContainer.from_dict(chain_container)
    return {
        "msa_container": {
            "_residue_container": residue_container,
            "_chain_container": chain_container,
        }
    }
