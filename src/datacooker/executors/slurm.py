"""SLURM executor: one job array per size tier, merge + index chained by afterok.

Each tier-shard's exact command is written to a ``.cmds`` manifest (one line per
shard); a small generic sbatch script runs the line for its
``SLURM_ARRAY_TASK_ID`` inside the project's pixi env. Nothing about the recipe
or tier is hard-coded in shell -- the driver builds the commands in Python.
"""
from __future__ import annotations

import os
import shlex
import subprocess
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from datacooker.executors.base import TaskHandle


@dataclass
class SlurmExecutor:
    """Submit each tier as a SLURM job array, chaining merge + index via afterok."""

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
        """Normalise paths and create the workdir/log directories."""
        self.workdir = Path(self.workdir)
        self.repo = Path(self.repo)
        self.workdir.mkdir(parents=True, exist_ok=True)
        if self.log_dir is None:
            self.log_dir = self.workdir / "logs"
        Path(self.log_dir).mkdir(parents=True, exist_ok=True)
        if self.pixi_manifest is None:
            self.pixi_manifest = self.repo / "pyproject.toml"

    # ---- public API -------------------------------------------------------
    def run_tier_shards(self, *, name: str, argv_for_shard: Callable[[int], list[str]],
                        n_shards: int, mem_gb: int, cores: int,
                        depends_on: tuple[str, ...] = ()) -> TaskHandle:
        """Submit one size tier's shards as a single SLURM job array."""
        manifest = self.workdir / f"{name}.cmds"
        with manifest.open("w", encoding="utf-8") as f:
            for i in range(n_shards):
                f.write(" ".join(shlex.quote(a) for a in argv_for_shard(i)) + "\n")
        script = self._write_script(name, manifest, array=True)
        job_id = self._sbatch(script, name, mem_gb, cores, depends_on,
                              array=f"0-{n_shards - 1}")
        return TaskHandle(name, job_id, depends_on)

    def run_once(self, *, name: str, argv: list[str], mem_gb: int, cores: int,
                 depends_on: tuple[str, ...] = ()) -> TaskHandle:
        """Submit a single command (merge, index) after its dependencies."""
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
            f"job_status_dir={shlex.quote(str(self.log_dir))}\n"
            "record_job_exit() {\n"
            "  job_exit_code=$?\n"
            "  trap - EXIT\n"
            "  set +e\n"
            "  printf '{\"job_id\":\"%s\",\"exit_code\":%d}\\n' \"${SLURM_JOB_ID}\" \"$job_exit_code\" "
            '> "$job_status_dir/${SLURM_JOB_ID}.exit.json.tmp"\n'
            '  mv -- "$job_status_dir/${SLURM_JOB_ID}.exit.json.tmp" "$job_status_dir/${SLURM_JOB_ID}.exit.json"\n'
            '  exit "$job_exit_code"\n'
            "}\n"
            "trap record_job_exit EXIT\n"
            'export PYTHONPATH="${PYTHONPATH:-.}"\n'
            f'CMD=$(sed -n "{line}p" {shlex.quote(str(manifest))})\n'
            'echo "host=$(hostname) task=${SLURM_ARRAY_TASK_ID:-0} start=$(date)"\n'
            f'pixi run --frozen --manifest-path {shlex.quote(str(self.pixi_manifest))} bash -c "$CMD"\n'
            'echo "end=$(date)"\n',
            encoding="utf-8",
        )
        return script

    def _sbatch(self, script: Path, name: str, mem_gb: int, cores: int,
                depends_on: tuple[str, ...], array: str | None) -> str | None:
        cmd = ["sbatch", "--parsable", f"--mem={mem_gb}g", f"-c{cores}"]
        if excluded := os.environ.get("DATACOOKER_EXCLUDE_NODES"):
            cmd.append(f"--exclude={excluded}")
        if array is not None:
            cmd.append(f"--array={array}")
        dep = [d for d in depends_on if d]
        if dep:
            cmd.append("--dependency=afterok:" + ":".join(dep))
        cmd.append(str(script))
        if not self.submit:
            self._dry_seq += 1
            fake = f"DRY{self._dry_seq}"
            print(f"[dry-run] {name}: {' '.join(cmd)}  -> {fake}"  # noqa: T201
                  + (f"  (afterok {dep})" if dep else ""))
            return fake
        out = subprocess.run(cmd, capture_output=True, text=True, check=True)  # noqa: S603
        job_id = out.stdout.strip().split(";")[0]
        script.with_suffix(".job-id").write_text(job_id + "\n", encoding="utf-8")
        print(f"[slurm] {name}: job {job_id}"  # noqa: T201
              + (f" afterok {dep}" if dep else "")
              + (f" array {array}" if array else ""))
        return job_id


def wait_for_jobs(job_ids: tuple[str, ...], *, poll_seconds: float = 30.0) -> None:
    """Wait for upstream terminal jobs before reading their generated inputs.

    ``afterok`` alone cannot order the planner's local filesystem reads. A pipeline
    invocation with dependencies therefore stays alive until those inputs are final.
    Unknown jobs and failed terminal states fail explicitly; neither means success.
    """
    import re
    import time

    pending = set(job_ids)
    active_states = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED", "REQUEUED", "RESIZING"}
    while pending:
        for job_id in sorted(pending):
            result = subprocess.run(  # noqa: S603
                ["scontrol", "show", "job", "-o", job_id],  # noqa: S607
                capture_output=True, text=True, check=False,
            )
            states = re.findall(r"\bJobState=(\S+)", result.stdout)
            reasons = re.findall(r"\bReason=(\S+)", result.stdout)
            if result.returncode or not states:
                msg = f"Cannot confirm upstream job {job_id}: {result.stderr.strip() or result.stdout.strip()}"
                raise RuntimeError(msg)
            if "DependencyNeverSatisfied" in reasons:
                msg = f"Upstream job {job_id} cannot run: DependencyNeverSatisfied."
                raise RuntimeError(msg)
            if all(state == "COMPLETED" for state in states):
                pending.remove(job_id)
                continue
            failed = set(states) - active_states - {"COMPLETED"}
            if failed:
                msg = f"Upstream job {job_id} failed: {', '.join(sorted(failed))}."
                raise RuntimeError(msg)
        if pending:
            print(f"[planner] waiting for upstream jobs: {', '.join(sorted(pending))}", flush=True)  # noqa: T201
            time.sleep(poll_seconds)
