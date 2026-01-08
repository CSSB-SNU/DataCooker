
from collections.abc import Callable

import numpy as np
from numpy import ndarray
from pipelines.cifmol import CIFMol

def neighbor_list_grid(  # noqa: PLR0915
    xyz: np.ndarray,
    d_thr: float,
    n_max: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute neighbor list using grid-based spatial partitioning."""
    n_atom = xyz.shape[0]
    nbrs = np.full((n_atom, n_max), -1, dtype=np.int64)
    counts = np.zeros(n_atom, dtype=np.int32)

    # 1) Mask invalid atoms (any NaN)
    valid = np.all(np.isfinite(xyz), axis=1)
    if not np.any(valid):
        return nbrs, counts

    # 2) Compressed array of valid points
    valid_xyz = xyz[valid]
    n_valid = valid_xyz.shape[0]

    # 3) Discretize into cells of side length d_thr (int64 coordinates)
    cell = np.floor(valid_xyz / d_thr).astype(np.int64)  # (n_valid, 3)

    # 4) Group points by cell via lexicographic sort on (x, y, z)
    order = np.lexsort((cell[:, 2], cell[:, 1], cell[:, 0]))
    cell_sorted = cell[order]

    # 5) Unique cells and their spans [start, end) in the sorted index space
    if n_valid > 1:
        change = np.any(np.diff(cell_sorted, axis=0) != 0, axis=1)
        starts = np.concatenate(([0], np.nonzero(change)[0] + 1))
    else:
        starts = np.array([0], dtype=np.int64)
    ends = np.concatenate((starts[1:], [n_valid]))
    unique_cells = cell_sorted[starts]  # (n_unique, 3)
    n_unique = unique_cells.shape[0]

    # 6) Helper: view (n,3) int64 as a structured dtype for consistent numeric lex compare
    #    (little-endian int64 x,y,z). This matches the lexsort order above.
    def as_struct3(a_int64x3: np.ndarray) -> np.ndarray:
        a = np.ascontiguousarray(a_int64x3)
        dt = np.dtype([("x", "<i8"), ("y", "<i8"), ("z", "<i8")])
        return a.view(dt).ravel()

    unique_struct = as_struct3(unique_cells)

    # 7) Precompute inverse map from compressed->sorted if needed later
    inv_order = np.empty(n_valid, dtype=np.int64)
    inv_order[order] = np.arange(n_valid)

    # 8) 27 neighbor-cell offsets (-1,0,1)^3
    offsets = (
        np.array(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij"))
        .reshape(3, -1)
        .T
    )  # (27, 3)

    # 9) Accumulate candidate (i,j) pairs in compressed indices
    pair_i = []
    pair_j = []

    # For each offset, find matching neighbor cells and generate Cartesian pairs
    for off in offsets:
        # neighbor cells for ALL existing unique cells under this offset
        nei_cells = unique_cells + off  # (n_unique, 3)
        nei_struct = as_struct3(nei_cells)

        # Search where these neighbor cells would be, under the SAME structured ordering
        pos = np.searchsorted(unique_struct, nei_struct, side="left")
        in_bounds = pos < n_unique

        ok = np.zeros_like(in_bounds, dtype=bool)
        if np.any(in_bounds):
            ok[in_bounds] = unique_struct[pos[in_bounds]] == nei_struct[in_bounds]

        if not np.any(ok):
            continue

        # Source cell ids are those indices where a neighbor cell exists
        src_c = np.nonzero(ok)[0]  # indices in [0..n_unique)
        dst_c = pos[ok]  # matching neighbor cell ids

        # Build index ranges for points in each cell's span (sorted space indices)
        src_ranges = [np.arange(starts[c], ends[c], dtype=np.int64) for c in src_c]
        dst_ranges = [np.arange(starts[c], ends[c], dtype=np.int64) for c in dst_c]

        if len(src_ranges) == 0:
            continue

        # Cartesian product per (src_cell, dst_cell) pair (vectorized at cell level)
        src_idx_sorted = np.concatenate(
            [
                np.repeat(r, len(d))
                for r, d in zip(src_ranges, dst_ranges, strict=False)
            ],
        )
        if src_idx_sorted.size == 0:
            continue
        dst_idx_sorted = np.concatenate(
            [np.tile(d, len(r)) for r, d in zip(src_ranges, dst_ranges, strict=False)],
        )

        # Map back from sorted space → compressed (unsorted) space
        src_idx = order[src_idx_sorted]
        dst_idx = order[dst_idx_sorted]

        # Drop self-pairs
        keep = src_idx != dst_idx
        if not np.any(keep):
            continue

        pair_i.append(src_idx[keep])
        pair_j.append(dst_idx[keep])

    if not pair_i:
        # No candidate pairs; return empty neighbor lists
        return nbrs, counts

    pair_i = np.concatenate(pair_i)
    pair_j = np.concatenate(pair_j)

    # 10) Distance filtering: keep pairs with ||valid_xyz[i]-valid_xyz[j]|| <= d_thr
    dvec = valid_xyz[pair_i] - valid_xyz[pair_j]
    dist2 = np.einsum("ij,ij->i", dvec, dvec)
    keep = dist2 <= (d_thr * d_thr)
    if not np.any(keep):
        return nbrs, counts

    pair_i = pair_i[keep]
    pair_j = pair_j[keep]

    # 11) Remove duplicate pairs (same (i,j) can appear via multiple offsets)
    ij = np.stack([pair_i, pair_j], axis=1).astype(np.int64)
    ij_packed = ij.view(np.dtype((np.void, ij.dtype.itemsize * 2))).ravel()
    uniq_idx = np.unique(ij_packed, return_index=True)[1]
    ij = ij[uniq_idx]
    pair_i, pair_j = ij[:, 0], ij[:, 1]

    # 12) Map compressed indices back to original n_atom-space
    valid_true_idx = np.flatnonzero(valid)
    gi = valid_true_idx[pair_i]
    gj = valid_true_idx[pair_j]

    # 13) Fill neighbor matrix: group by gi, keep up to n_max in order of (gi, gj)
    order_fill = np.lexsort((gj, gi))
    gi = gi[order_fill]
    gj = gj[order_fill]

    uniq_i, first_pos, counts_all = np.unique(gi, return_index=True, return_counts=True)
    take_counts = np.minimum(counts_all, n_max)

    if uniq_i.size > 0:
        gather_idx = np.concatenate(
            [
                np.arange(s, s + t, dtype=np.int64)
                for s, t in zip(first_pos, take_counts, strict=True)
            ],
        )
        gi_take = gi[gather_idx]
        gj_take = gj[gather_idx]

        # Relative slot [0..taken-1] within each group
        rel = np.arange(gj.size, dtype=np.int64) - np.repeat(first_pos, counts_all)
        rel = rel[gather_idx]

        nbrs[gi_take, rel] = gj_take
        counts[uniq_i] = take_counts.astype(np.int32)

    return nbrs, counts


def cdist_clipped(
    xyz1: ndarray,
    xyz2: ndarray | None = None,
    d_thr: float = 32.0,
    n_max: int = 128,  # max neighbors per atom
) -> ndarray:
    """Compute a dense (n1, n2) distance map clipped at d_thr using neighbor_list_grid."""
    n1 = xyz1.shape[0]
    n2 = 0 if xyz2 is None else xyz2.shape[0]
    # Dense map clipped at threshold
    dist = np.full((n1, n2), d_thr, dtype=xyz1.dtype)

    # Combine xyz1 + xyz2 (neighbor_list_grid works on one array)
    xyz = np.concatenate([xyz1, xyz2], axis=0) if xyz2 is not None else xyz1

    # Build neighbor list over combined space
    nbrs, _ = neighbor_list_grid(xyz, d_thr, n_max)  # (n1+n2, n_max)

    # Only keep neighbors from xyz1 → xyz2
    valid_mask = (nbrs >= n1) & (nbrs < n1 + n2)
    if not np.any(valid_mask):
        return dist

    src_idx = np.broadcast_to(
        np.arange(n1 + n2, dtype=np.int64)[:, None],
        nbrs.shape,
    )[valid_mask]

    dst_idx = nbrs[valid_mask]

    mask_1_to_2 = (src_idx < n1) & (dst_idx >= n1)
    if not np.any(mask_1_to_2):
        return dist

    i = src_idx[mask_1_to_2]
    j = dst_idx[mask_1_to_2] - n1  # shift to [0..n2)

    dvec = xyz1[i] - xyz2[j]
    dist_ij = np.sqrt(np.einsum("ij,ij->i", dvec, dvec))

    dist[i, j] = dist_ij

    return dist


def pdist_clipped(
    xyz: ndarray,
    d_thr: float = 32.0,
    n_max: int = 128,
) -> ndarray:
    """Compute a dense (n_atom, n_atom) distance map clipped at d_thr using neighbor_list_grid."""
    n_atom = xyz.shape[0]

    # Initialize with clipped distance
    dist = np.full((n_atom, n_atom), d_thr, dtype=xyz.dtype)
    np.fill_diagonal(dist, 0.0)

    # Use existing neighbor list (fast spatial grid)
    nbrs, _ = neighbor_list_grid(xyz, d_thr, n_max)  # (n_atom, n_max), -1 padded

    # Gather all valid (i, j) pairs in one shot
    valid_mask = nbrs != -1
    if not np.any(valid_mask):
        return dist

    row_idx = np.broadcast_to(
        np.arange(n_atom, dtype=np.int64)[:, None],
        nbrs.shape,
    )[valid_mask]
    col_idx = nbrs[valid_mask]

    # Compute exact distances for neighbor pairs
    dvec = xyz[row_idx] - xyz[col_idx]
    dist_ij = np.sqrt(np.einsum("ij,ij->i", dvec, dvec))

    # Write symmetric distances (no Python loops)
    dist[row_idx, col_idx] = dist_ij
    dist[col_idx, row_idx] = dist_ij

    return dist

def get_shortest_distances(
    atom_pos: ndarray,
    atom_pos_mask: ndarray,
    atom_to_res_idx: ndarray,
    min_distance: float = 2.0,
    max_distance: float = 22.0,
) -> tuple[ndarray, ndarray]:
    """Compute residue-level shortest distances from atom positions (single example)."""
    L, _ = atom_pos.shape

    # 1) Atom-level pairwise distances (L, L)
    dist = pdist_clipped(
        atom_pos,
        d_thr=max_distance,
    )

    # Apply atom mask
    mask_i = atom_pos_mask[:, None]  # (L, 1)
    mask_j = atom_pos_mask[None, :]  # (1, L)
    valid_atom_mask = mask_i & mask_j  # (L, L)

    dist = np.where(~valid_atom_mask, max_distance, dist)
    dist = np.clip(dist, min_distance, max_distance)

    # 2) Build residue existence mask (R_max)
    R_max = int(atom_to_res_idx.max().item()) + 1

    residue_exists = np.zeros((R_max,), dtype=bool)
    valid_res_idx = atom_to_res_idx[atom_pos_mask]
    residue_exists[valid_res_idx] = True

    # Residue pair mask: both residues must exist
    residue_mask = residue_exists[:, None] & residue_exists[None, :]

    # 3) Aggregate shortest distances to residue level using scatter-reduce (min)
    ri = np.broadcast_to(atom_to_res_idx[:, None], (L, L))
    rj = np.broadcast_to(atom_to_res_idx[None, :], (L, L))
    pair_idx = ri * R_max + rj  # (L, L)

    block_size = R_max * R_max

    scatter_idx_flat = pair_idx.reshape(-1)  # (L * L,)
    src = dist.reshape(-1)  # (L * L,)

    out = np.full(
        (block_size,),
        max_distance,
        dtype=dist.dtype,
    )

    np.minimum.at(out, scatter_idx_flat, src)

    residue_dists = out.reshape(R_max, R_max)

    return residue_dists, residue_mask

def _residue_distance_profile(
    cifmol: CIFMol,
    bins: np.ndarray,
    min_distance: float,
    max_distance: float,
    seq_sep: int,
) -> dict[tuple[str, str], np.ndarray]:
    """
    Compute residue–residue shortest-distance counts per chem_comp_id pair for a CIFMol.

    Returns
    -------
    dict
        {(chem_a, chem_b): counts ndarray}, where chem_a <= chem_b lexicographically.
    """
    xyz = np.asarray(cifmol.atoms.xyz.value, dtype=np.float32)
    atom_mask = np.isfinite(xyz).all(axis=1)
    if not np.any(atom_mask):
        return {}

    atom_to_res_idx = np.asarray(cifmol.index_table.atom_to_res, dtype=np.int64)
    res_to_chain = np.asarray(cifmol.index_table.res_to_chain, dtype=np.int64)
    chem_comp_ids = np.asarray(cifmol.residues.chem_comp_id.value, dtype=str)

    def _to_numeric_seq_ids(raw_seq_ids: np.ndarray) -> np.ndarray:
        seq_ids = np.full(raw_seq_ids.shape, np.nan, dtype=np.float64)
        for idx, val in enumerate(raw_seq_ids):
            try:
                seq_ids[idx] = float(val)
            except (TypeError, ValueError):
                continue
        return seq_ids

    seq_ids = _to_numeric_seq_ids(np.asarray(cifmol.residues.cif_idx.value))

    residue_dists, residue_mask = get_shortest_distances(
        atom_pos=xyz,
        atom_pos_mask=atom_mask,
        atom_to_res_idx=atom_to_res_idx,
        min_distance=min_distance,
        max_distance=max_distance,
    )

    # exclude too close
    if seq_ids.shape[0] == residue_mask.shape[0]:
        seq_i = seq_ids[:, None]
        seq_j = seq_ids[None, :]
        chain_i = res_to_chain[:, None]
        chain_j = res_to_chain[None, :]
        close_mask = (
            (chain_i == chain_j)
            & np.isfinite(seq_i)
            & np.isfinite(seq_j)
            & (np.abs(seq_i - seq_j) < seq_sep)
        )
        residue_mask = residue_mask & ~close_mask

    if residue_dists.size == 0 or chem_comp_ids.size == 0:
        return {}

    n_res = residue_dists.shape[0]
    tri_mask = np.triu(np.ones((n_res, n_res), dtype=bool), k=1)
    valid_mask = tri_mask & residue_mask
    if not np.any(valid_mask):
        return {}

    ri_idx, rj_idx = np.nonzero(valid_mask)
    distances = residue_dists[ri_idx, rj_idx]
    if distances.size == 0:
        return {}

    chem_i = chem_comp_ids[ri_idx]
    chem_j = chem_comp_ids[rj_idx]
    pair_labels = np.sort(
        np.stack([chem_i, chem_j], axis=1),
        axis=1,
    )
    pair_keys = np.char.add(np.char.add(pair_labels[:, 0], "::"), pair_labels[:, 1])

    profile: dict[tuple[str, str], np.ndarray] = {}
    uniq_keys, inv = np.unique(pair_keys, return_inverse=True)
    for key_idx, key_str in enumerate(uniq_keys):
        mask = inv == key_idx
        if not np.any(mask):
            continue
        hist, _ = np.histogram(distances[mask], bins=bins)
        chem_a, chem_b = str(key_str).split("::", maxsplit=1)
        profile[(chem_a, chem_b)] = hist.astype(np.int64)
    return profile


def residue_distance_profile() -> Callable[..., dict[tuple[str, str], np.ndarray]]:
    """Instruction wrapper for residue–residue distance histograms per CIFMol."""

    def _worker(
        cifmol: CIFMol,
        bins: np.ndarray,
        min_distance: float,
        max_distance: float,
        seq_sep: int,
    ) -> dict[tuple[str, str], np.ndarray]:
        return _residue_distance_profile(
            cifmol=cifmol,
            bins=bins,
            min_distance=min_distance,
            max_distance=max_distance,
            seq_sep=seq_sep,
        )

    return _worker
