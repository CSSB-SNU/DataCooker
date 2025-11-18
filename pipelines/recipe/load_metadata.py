from pathlib import Path

from datacooker import RecipeBook
from pipelines.instructions.metadata_instructions import (
    load_signalp,
    load_tsv,
    reverse_dict,
)

"""Rebuild a CIF lmdb to train AF3"""

metadata_recipe = RecipeBook()


metadata_recipe.add(
    targets=[
        (("seqID2seq", dict),),
        (("clusterID2seqID", dict),),
    ],
    instruction=load_tsv,
    inputs=[
        {
            "kwargs": {
                "tsv_file_path": ("seqID2seq_path", Path),
            },
            "params": {
                "split_by_comma": False,
            },
        },
        {
            "kwargs": {
                "tsv_file_path": ("clusterID2seqID_path", Path),
            },
            "params": {
                "split_by_comma": True,
            },
        },
    ],
)

metadata_recipe.add(
    targets=[
        (("seq2seqID", dict),),
        (("seqID2clusterID", dict),),
    ],
    instruction=reverse_dict,
    inputs=[
        {
            "kwargs": {
                "input_dict": ("seqID2seq", dict),
            },
        },
        {
            "kwargs": {
                "input_dict": ("clusterID2seqID", dict),
            },
        },
    ],
)


metadata_recipe.add(
    targets=[(("signalp_dict", dict),)],
    instruction=load_signalp,
    inputs=[
        {
            "kwargs": {
                "signalp_dir": ("signalp_dir", Path | None),
            },
        },
    ],
)

RECIPE = metadata_recipe
TARGETS = ["seq2seqID", "seqID2clusterID", "signalp_dict"]

