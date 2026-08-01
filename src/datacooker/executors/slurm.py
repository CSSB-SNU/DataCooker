"""SLURM executor: one job array per size tier, merge + index chained by afterok.

Each tier-shard's exact command is written to a ``.cmds`` manifest (one line per
shard); a small generic sbatch script runs the line for its
``SLURM_ARRAY_TASK_ID`` inside the project's pixi env. Nothing about the recipe
or tier is hard-coded in shell -- the driver builds the commands in Python.
"""
from __future__ import annotations

import shlex
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from datacooker.executors.base import TaskHandle


@dataclass
class SlurmExecutor:
    workdir: Path                       # scratch for manifests, scripts, logs
    repo: Path                          # cwd + pixi manifest root
    # cpu-long-q is the only uncapped cpu QOS (GrpTRES node=0 = unlimited);
    # cpu-standard-q is throttled to ~3-4 nodes, which needlessly serialized an
    # 8-shard array into 2 waves. cpu-long's 30-day limit is a non-issue for builds.
    partition: str = "cpu-long"
    qos: str = "cpu-long-q"
    log_dir: Path | None = None
    pixi_manifest: Path | None = None   # defaults to repo/pyproject.toml
    submit: bool = True                 # False -> dry run (write scripts, no sbatch)
    _dry_seq: int = 0

    def __post_init__(self) -> None:
        self.workdir = Path(self.workdir)
        self.repo = Path(self.repo)
        self.workdir.mkdir(parents=True, exist_ok=True)
        if self.log_dir is None:
            self.log_dir = self.workdir / "logs"
        Path(self.log_dir).mkdir(parents=True, exist_ok=True)
        if self.pixi_manifest is None:
            self.pixi_manifest = self.repo / "pyproject.toml"

    # ---- public API -------------------------------------------------------
    def run_tier_shards(self, *, name, argv_for_shard: Callable[[int], list[str]],
                        n_shards: int, mem_gb: int, cores: int,
                        depends_on: tuple[str, ...] = ()) -> TaskHandle:
        manifest = self.workdir / f"{name}.cmds"
        with manifest.open("w", encoding="utf-8") as f:
            for i in range(n_shards):
                f.write(" ".join(shlex.quote(a) for a in argv_for_shard(i)) + "\n")
        script = self._write_script(name, manifest, array=True)
        job_id = self._sbatch(script, name, mem_gb, cores, depends_on,
                              array=f"0-{n_shards - 1}")
        return TaskHandle(name, job_id, depends_on)

    def run_once(self, *, name, argv: list[str], mem_gb: int, cores: int,
                 depends_on: tuple[str, ...] = ()) -> TaskHandle:
        manifest = self.workdir / f"{name}.cmds"
        manifest.write_text(" ".join(shlex.quote(a) for a in argv) + "\n", encoding="utf-8")
        script = self._write_script(name, manifest, array=False)
        job_id = self._sbatch(script, name, mem_gb, cores, depends_on, array=None)
        return TaskHandle(name, job_id, depends_on)

    # ---- internals --------------------------------------------------------
    def _write_script(self, name: str, manifest: Path, *, array: bool) -> Path:
        out_pat = "%x_%A_%a" if array else "%x_%j"
        line = "$((SLURM_ARRAY_TASK_ID + 1))" if array else "1"
        script = self.workdir / f"{name}.sbatch"
        script.write_text(
            "#!/bin/bash\n"
            f"#SBATCH -J {name}\n"
            f"#SBATCH -p {self.partition}\n"
            f"#SBATCH --qos={self.qos}\n"
            f"#SBATCH -o {self.log_dir}/{out_pat}.out\n"
            f"#SBATCH -e {self.log_dir}/{out_pat}.err\n"
            "set -euo pipefail\n"
            f"cd {shlex.quote(str(self.repo))}\n"
            'export PYTHONPATH="."\n'
            f"CMD=$(sed -n \"{line}p\" {shlex.quote(str(manifest))})\n"
            'echo "host=$(hostname) task=${SLURM_ARRAY_TASK_ID:-0} start=$(date)"\n'
            f"pixi run --manifest-path {shlex.quote(str(self.pixi_manifest))} bash -c \"$CMD\"\n"
            'echo "end=$(date)"\n',
            encoding="utf-8",
        )
        return script

    def _sbatch(self, script: Path, name: str, mem_gb: int, cores: int,
                depends_on: tuple[str, ...], array: str | None) -> str | None:
        cmd = ["sbatch", "--parsable", f"--mem={mem_gb}g", f"-c{cores}"]
        if array is not None:
            cmd.append(f"--array={array}")
        dep = [d for d in depends_on if d]
        if dep:
            cmd.append("--dependency=afterok:" + ":".join(dep))
        cmd.append(str(script))
        if not self.submit:
            self._dry_seq += 1
            fake = f"DRY{self._dry_seq}"
            print(f"[dry-run] {name}: {' '.join(cmd)}  -> {fake}"
                  + (f"  (afterok {dep})" if dep else ""))
            return fake
        out = subprocess.run(cmd, capture_output=True, text=True, check=True)  # noqa: S603
        job_id = out.stdout.strip().split(";")[0]
        print(f"[slurm] {name}: job {job_id}"
              + (f" afterok {dep}" if dep else "")
              + (f" array {array}" if array else ""))
        return job_id
