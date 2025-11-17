from pathlib import Path

import lmdb
import numpy as np
from pipelines.utils.convert import from_bytes

def load_keys(
    env_path: Path,
) -> list[str]:
    """
    Load all keys from LMDB database.

    Args:
        env_path: Path to the LMDB environment.
    Returns
    -------
        list[str]
            List of all keys in the LMDB database.
    """
    env = lmdb.open(str(env_path), readonly=True, lock=False)
    with env.begin() as txn:
        cursor = txn.cursor()
        keys = [key.decode() for key in cursor.iternext(keys=True, values=False)]
    env.close()
    return keys

def load_data(
    env_path: Path,
    key: str|None = None, # if None, load first key
) -> dict:
    """
    Load data from LMDB database.

    Args:
        env_path: Path to the LMDB environment.
        key: Key to load. If None, load the first key.
    Returns.
    -------
        dict
            The loaded data dictionary.
    """
    env = lmdb.open(str(env_path), readonly=True, lock=False)
    with env.begin() as txn:
        if key is None:
            cursor = txn.cursor()
            if cursor.first():              # Move to first key-value pair
                key, value = cursor.item()
            else:
                msg = "LMDB is empty."
                raise KeyError(msg)
        else:
            value = txn.get(key.encode())
    env.close()
    if value is None:
        msg = f"Key '{key}' not found in LMDB database."
        raise KeyError(msg)
    return from_bytes(value)

def deep_equal(a, b, path=""):
    mismatches = []

    def _cmp(x, y, p):
        # --- 0. Type mismatch ---
        if type(x) != type(y):
            mismatches.append((p, f"Type mismatch: {type(x)} != {type(y)}"))
            return

        # --- 1. Dict ---
        if isinstance(x, dict):
            a_keys = set(x.keys())
            b_keys = set(y.keys())

            diff = a_keys ^ b_keys
            if diff:
                mismatches.append((p, f"Key mismatch: {diff}"))

            for k in a_keys & b_keys:
                new_p = f"{p}.{k}" if p else k
                _cmp(x[k], y[k], new_p)
            return

        # --- 2. List / Tuple ---
        if isinstance(x, (list, tuple)):
            if len(x) != len(y):
                mismatches.append((p, f"Length mismatch: {len(x)} != {len(y)}"))

            for i in range(min(len(x), len(y))):
                new_p = f"{p}[{i}]"
                _cmp(x[i], y[i], new_p)
            return

        # --- 3. NumPy array ---
        if isinstance(x, np.ndarray):
            if x.shape != y.shape:
                mismatches.append((p, f"Shape mismatch: {x.shape} != {y.shape}"))
                return

            # float array → NaN-aware 비교
            if np.issubdtype(x.dtype, np.floating):
                equal_mask = (x == y) | (np.isnan(x) & np.isnan(y))
            else:
                # non-float array → 일반 비교
                equal_mask = (x == y)

            if not np.all(equal_mask):
                idxs = np.argwhere(~equal_mask)
                mismatches.append((p, f"Array mismatch at indices: {idxs.tolist()}"))
            return

        # --- 4. Primitive 값 (int/float/str/etc.) ---
        # float일 때 NaN 비교 보정
        if isinstance(x, float) and isinstance(y, float):
            if np.isnan(x) and np.isnan(y):
                return
            if x != y:
                mismatches.append((p, f"Value mismatch: {x} != {y}"))
            return

        # 기본 비교
        if x != y:
            mismatches.append((p, f"Value mismatch: {x} != {y}"))
        return

    _cmp(a, b, path)

    ok = len(mismatches) == 0
    return ok, mismatches

if __name__ == "__main__":
    old_DB_path = Path("/public_data/BioMolDBv2_2024Oct21/cif.lmdb")
    new_DB_path = Path("/public_data/BioMolDBv2_2024Oct21/cif_new_shard0.lmdb")

    old_keys = load_keys(old_DB_path)
    new_keys = load_keys(new_DB_path)
    new_data = load_data(new_DB_path, "109d")
    old_data = load_data(old_DB_path, "109d")
    is_equal, message = deep_equal(old_data, new_data)
    breakpoint()

    if set(old_keys) != set(new_keys):
        print("Key sets differ.")
        print("Only in old DB:", set(old_keys) - set(new_keys))
        print("Only in new DB:", set(new_keys) - set(old_keys))
        exit(1)

    for key in old_keys:
        old_data = load_data(old_DB_path, key)
        new_data = load_data(new_DB_path, key)
        is_equal, message = deep_equal(old_data, new_data)
        if not is_equal:
            print(f"Data mismatch for key '{key}': {message}")
            breakpoint()
            exit(1)
    breakpoint()
