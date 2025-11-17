from datacooker import RecipeBook
from pipelines.cifmol import CIFMol
from pipelines.instructions.seq_instructions import build_fasta

"""Build a CIFMol->fasta Cooker."""

fasta_recipe = RecipeBook()

fasta_recipe.add(
    targets=[
        ("fasta", str),
    ],
    instruction=build_fasta(),
    inputs={
        "kwargs": {
            "cifmol": ("cifmol", CIFMol | None),
        },
    },
)

RECIPE = fasta_recipe
TARGETS = ["fasta"]
