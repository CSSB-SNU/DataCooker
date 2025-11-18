import logging
import sys
from importlib import import_module
from pathlib import Path
from typing import Any

import click
import lmdb
from omegaconf import OmegaConf

from pipelines.utils.lmdb import rebuild_lmdb

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(name)s - %(message)s",
    stream=sys.stdout,
    force=True,
)
OmegaConf.register_new_resolver("p", lambda x: Path(x))

def dotted_to_obj(path: str) -> Any:
    """Convert a dotted path string to a Python object."""
    module_name, attr_name = path.rsplit(".", 1)
    module = import_module(module_name)
    return getattr(module, attr_name)


def load_config(config_path: Path) -> dict[str, Any]:
    """Load configuration from a YAML file using OmegaConf."""
    cfg = OmegaConf.load(config_path)
    cfg_dict = OmegaConf.to_container(cfg, resolve=True)
    if not isinstance(cfg_dict, dict):
        msg = "Configuration file must contain a dictionary at the top level."
        raise TypeError(msg)
    config: dict[str, Any] = dict(cfg_dict)

    if "load_func" in config and isinstance(config["load_func"], str):
        config["load_func"] = dotted_to_obj(config["load_func"])
    if "transform_func" in config and isinstance(config["transform_func"], str):
        config["transform_func"] = dotted_to_obj(config["transform_func"])

    return config

def load_data_list(data_dir: Path, pattern: str = "*.cif*") -> list[Path]:
    """Load a list of data file paths from a directory."""
    return list(data_dir.rglob(pattern))


# ==============================================================
# Command Group
# ==============================================================
@click.group()
def cli() -> None:
    """Build and merge LMDB databases from CIF files."""


# ==============================================================
# 1. Build Command
# ==============================================================
@cli.command("rebuild")
@click.argument("config", type=click.Path(exists=True, path_type=Path))
@click.option("--map-size", "-m", type=float, default=1e10, show_default=True)
def rebuild(
    config: Path,
    map_size: float,
) -> None:
    """
    Rebuild an LMDB database from an existing LMDB database with transformations.

    Example:
        python build_lmdb.py rebuild --config configs/ccd_lmdb.yaml --map-size 1e12 --shard-idx 0 --n-shards 4
    """
    map_size = int(map_size)
    config = load_config(config)

    rebuild_lmdb(
        old_env_path=config["old_env_path"],
        new_env_path=config["env_path"],
        **config,
        map_size=map_size,
    )

    # check the key count
    env = lmdb.open(str(config["env_path"]), readonly=True, lock=False)
    with env.begin() as txn:
        cursor = txn.cursor()
        key_count = sum(1 for _ in cursor)
    env.close()
    click.echo(f"[Done] Built LMDB at {config['env_path']} with {key_count} keys.")

if __name__ == "__main__":
    cli()
