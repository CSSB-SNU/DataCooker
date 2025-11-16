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
    """
    Deeply compare two nested dicts/lists/tuples/arrays.
    Returns (True, "") if equal, else (False, explanation_string)
    """
    # Type check
    if type(a) != type(b):
        return False, f"Type mismatch at {path}: {type(a)} != {type(b)}"

    # Dict 비교
    if isinstance(a, dict):
        a_keys = set(a.keys())
        b_keys = set(b.keys())
        if a_keys != b_keys:
            return False, f"Key mismatch at {path}: {a_keys ^ b_keys}"

        for k in a_keys:
            ok, msg = deep_equal(a[k], b[k], f"{path}.{k}")
            if not ok:
                return False, msg
        return True, ""

    # List/tuple 비교
    elif isinstance(a, (list, tuple)):
        if len(a) != len(b):
            return False, f"Length mismatch at {path}: {len(a)} != {len(b)}"

        for i, (x, y) in enumerate(zip(a, b)):
            ok, msg = deep_equal(x, y, f"{path}[{i}]")
            if not ok:
                return False, msg
        return True, ""

    # Numpy array 비교
    elif isinstance(a, np.ndarray):
        if not np.array_equal(a, b):
            return False, f"Array mismatch at {path}"
        return True, ""

    # Primitive 값 비교
    else:
        if a != b:
            return False, f"Value mismatch at {path}: {a} != {b}"
        return True, ""

if __name__ == "__main__":
    old_DB_path = Path("/public_data/CCD/biomol_CCD.lmdb")
    new_DB_path = Path("/public_data/CCD/biomol_CCD_test.lmdb")

    old_keys = load_keys(old_DB_path)
    new_keys = load_keys(new_DB_path)

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
