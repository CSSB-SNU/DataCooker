from pathlib import Path

def read_tsv_file(tsv_path: Path) -> dict[str, list[str]]:
    """
    Read a TSV (Tab-Separated Values) file and return its contents as a list of dictionaries.

    Args:
        tsv_path: Path to the TSV file.

    Returns
    -------
        list[dict[str, str]]
            List of rows represented as dictionaries with column headers as keys.
    """
    data: dict[str, list[str]] = {}
    with tsv_path.open("r", encoding="utf-8") as file:
        lines = file.readlines()
        for line in lines:
            key1, key2, value = line.strip().split("\t")
            key = f"{key1}_{key2}"
            data[key] = value.split(",")
    return data

def cal_degree_of_node(data: dict[str, list[str]]) -> dict[str, int]:
    """
    Calculate the degree of each node in the graph represented by the edge list.

    Args:
        data: Dictionary where keys are node identifiers and values are lists of connected nodes.

    Returns
    -------
        dict[str, int]
            Dictionary mapping each node to its degree.
    """
    degree_dict: dict[str, int] = {}
    for key, values in data.items():
        key1, key2 = key.split("_")
        if key1 not in degree_dict:
            degree_dict[key1] = 0
        if key2 not in degree_dict:
            degree_dict[key2] = 0
        degree_dict[key1] += len(values)
        degree_dict[key2] += len(values)
    return degree_dict

def remove_ligand_ligand_edges(data: dict[str, list[str]]) -> dict[str, list[str]]:
    """
    Remove edges that connect two ligand nodes.

    Args:
        data: Dictionary where keys are node identifiers and values are lists of connected nodes.

    Returns
    -------
        dict[str, list[str]]
            Filtered dictionary with ligand-ligand edges removed.
    """
    filtered_data: dict[str, list[str]] = {}
    for key, values in data.items():
        key1, key2 = key.split("_")
        if key1.startswith("cL") and key2.startswith("cL"):
            continue
        filtered_data[key] = values
    return filtered_data

def find_edges_that_include(data, key):
    all_edges = list(data.keys())
    included_edges = [edge for edge in all_edges if key in edge.split("_")]
    return included_edges

if __name__ == "__main__":
    tsv_path = Path("/home/psk6950/data/BioMolDBv2_2024Oct21/metadata/graph_edges.tsv")
    data = read_tsv_file(tsv_path)
    data = remove_ligand_ligand_edges(data)
    degree_dict = cal_degree_of_node(data)

    # sort by degree
    nodes = list(degree_dict.keys())
    degrees = [degree_dict[node] for node in nodes]
    sorted_indices = sorted(range(len(degrees)), key=lambda i: degrees[i], reverse=True)
    sorted_nodes = [nodes[i] for i in sorted_indices]
    sorted_degrees = [degrees[i] for i in sorted_indices]

    sorted_ligands_indices = [i for i, node in enumerate(sorted_nodes) if node.startswith("cL")]
    top50_ligands = [sorted_nodes[i] for i in sorted_ligands_indices[:50]]
    top50_ligand_degrees = [sorted_degrees[i] for i in sorted_ligands_indices[:50]]

    test = find_edges_that_include(data, "cL0021727")
    breakpoint()
