from collections.abc import Callable
from datetime import date

import numpy as np

from pipelines.cifmol import CIFMol
from pipelines.instructions.seq_instructions import extract_sequence_from_cifmol


def filter_by_resolution_and_date(
    resolution_cutoff: float = 8.0,
    date_cutoff: date = date(2099, 1, 1),
) -> Callable[[CIFMol|None], CIFMol|None]:
    """Filter instruction to select entries by resolution and date."""
    def worker(cifmol: CIFMol|None) -> CIFMol|None:
        if cifmol is None:
            return None

        resolution, deposition_date = cifmol.metadata["resolution"], cifmol.metadata["deposition_date"]
        deposition_date = date.fromisoformat(deposition_date)
        if resolution is not None and resolution <= resolution_cutoff and deposition_date < date_cutoff:
            return cifmol
        return None

    return worker

def filter_water(cifmol: CIFMol|None) -> CIFMol|None:
    """Filter instruction to remove water molecules from CIFMol."""
    if cifmol is None:
        return None
    water_mask = ~np.isin(cifmol.residues.chem_comp_id, ["HOH", "DOD"])
    if water_mask.sum() == 0:
        return None
    cifmol = cifmol.residues[water_mask].extract()

    if len(cifmol.chains) == 0:
        return None

    return cifmol

def filter_signalp(
    cifmol: CIFMol|None,
    seq2seqID_dict:dict[str,str],
    signalp_dict:dict[str,tuple[int,int]],
) -> dict|None:
    """Filter instruction to remove signal peptides from CIFMol."""
    if cifmol is None:
        return None
    seq_dict = extract_sequence_from_cifmol(cifmol)
    valid_residue_indices = []
    cursor = 0
    for chain_id, seq in seq_dict.items():
        seq_id = seq2seqID_dict[seq]
        chain_cifmol = cifmol.chains[cifmol.chains.chain_id == chain_id].extract()
        if seq_id not in signalp_dict:
            valid_residue_indices.extend(list(range(cursor, cursor + len(chain_cifmol.residues))))
            cursor += len(chain_cifmol.residues)
            continue
        _, signalp_end = signalp_dict[seq_id]
        valid_residue_indices.extend(list(range(cursor + signalp_end, cursor + len(chain_cifmol.residues))))
        cursor += len(chain_cifmol.residues)
    filtered_cifmol = cifmol.residues[valid_residue_indices].extract()
    return filtered_cifmol.to_dict()
