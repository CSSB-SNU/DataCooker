from __future__ import annotations

import subprocess
import unittest
from unittest.mock import patch

from datacooker.executors.slurm import wait_for_jobs


def result(state: str, reason: str = "None") -> subprocess.CompletedProcess:
    return subprocess.CompletedProcess([], 0, f"JobId=42 JobState={state} Reason={reason}", "")


class SlurmWaitTests(unittest.TestCase):
    def test_running_does_not_mean_ready(self) -> None:
        with patch("datacooker.executors.slurm.subprocess.run", side_effect=[result("RUNNING"), result("COMPLETED")]) as run:
            wait_for_jobs(("42",), poll_seconds=0)
            self.assertEqual(run.call_count, 2)

    def test_failed_upstream_aborts_planning(self) -> None:
        with patch("datacooker.executors.slurm.subprocess.run", return_value=result("FAILED")), self.assertRaisesRegex(RuntimeError, "FAILED"):
            wait_for_jobs(("42",), poll_seconds=0)

    def test_impossible_dependency_aborts_planning(self) -> None:
        with patch("datacooker.executors.slurm.subprocess.run", return_value=result("PENDING", "DependencyNeverSatisfied")), self.assertRaisesRegex(RuntimeError, "DependencyNeverSatisfied"):
            wait_for_jobs(("42",), poll_seconds=0)

    def test_unknown_state_is_not_success(self) -> None:
        with patch("datacooker.executors.slurm.subprocess.run", return_value=subprocess.CompletedProcess([], 1, "", "unknown job")), self.assertRaisesRegex(RuntimeError, "Cannot confirm"):
            wait_for_jobs(("42",), poll_seconds=0)
