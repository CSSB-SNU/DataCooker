import numpy as np
from datacooker import RecipeBook

from pipelines.cifmol import CIFMol
from pipelines.instructions.distance_profile import residue_distance_profile

"""Build a CIFMol->residue distance profile Cooker."""

distance_profile_recipe = RecipeBook()

distance_profile_recipe.add(
    targets=[
        ("residue_distance_profile", dict[tuple[str, str], np.ndarray]),
    ],
    instruction=residue_distance_profile(),
    inputs={
        "kwargs": {
            "cifmol": ("cifmol", CIFMol),
            "bins": ("bins", np.ndarray),
            "min_distance": ("min_distance", float),
            "max_distance": ("max_distance", float),
            "seq_sep": ("seq_sep", int),
        },
    },
)

RECIPE = distance_profile_recipe
TARGETS = ["residue_distance_profile"]
