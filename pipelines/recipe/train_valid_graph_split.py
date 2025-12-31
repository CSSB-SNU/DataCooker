from pathlib import Path

from networkx import Graph

from datacooker import RecipeBook
from pipelines.instructions.graph_cluster_instructions import (
    build_whole_graph,
    count_category_count,
    extract_edges,
    split_graph_by_components,
    split_train_valid,
)

"""Build a sequence clustering Cooker."""

tv_split_recipe = RecipeBook()

tv_split_recipe.add(
    targets=[
        ("whole_graph", Graph),
        ("polymer_graph", Graph),
    ],
    instruction=build_whole_graph,
    inputs={
        "kwargs": {
            "edge_tsv_path": ("edge_tsv_path", Path),
            "ignore_nodes": ("ignore_nodes", list),
        },
    },
)

tv_split_recipe.add(
    targets=[
        ("subgraphs", dict),
        ("edge_wo_ligand_counts", dict),
    ],
    instruction=split_graph_by_components,
    inputs={
        "kwargs": {
            "whole_graph": ("polymer_graph", Path),
        },
    },
)

tv_split_recipe.add(
    targets=[
        ("train_edges", list),
        ("valid_edges", list),
    ],
    instruction=split_train_valid,
    inputs={
        "kwargs": {
            "whole_graph": ("whole_graph", Graph),
            "subgraphs": ("subgraphs", dict),
            "edge_wo_ligand_counts": ("edge_wo_ligand_counts", dict),
            "train_ratio": ("train_ratio", float),
        },
    },
)

tv_split_recipe.add(
    targets=[
        (("train_edge_list", list),),
        (("valid_edge_list", list),),
    ],
    instruction=extract_edges,
    inputs=[
        {
            "kwargs": {
                "edge_tsv_path": ("edge_tsv_path", Path),
                "to_be_extracted": ("train_edges", list),
            },
        },
        {
            "kwargs": {
                "edge_tsv_path": ("edge_tsv_path", Path),
                "to_be_extracted": ("valid_edges", list),
            },
        },
    ],
)

tv_split_recipe.add(
    targets=[
        (("train_edge_statistics", list),),
        (("valid_edge_statistics", list),),
    ],
    instruction=count_category_count,
    inputs=[
        {
            "kwargs": {
                "edge_list": ("train_edges", list),
            },
        },
        {
            "kwargs": {
                "edge_list": ("valid_edges", list),
            },
        },
    ],
)


RECIPE = tv_split_recipe
TARGETS = ["train_edge_list", "valid_edge_list", "train_edge_statistics", "valid_edge_statistics"]
