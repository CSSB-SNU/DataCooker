from __future__ import annotations

import os
import time
import unittest
from typing import Any

from datacooker._ray import iter_ray_chunks

_STATE: Any = None
_INSTALLS = 0


def install(state: Any) -> None:
    global _STATE, _INSTALLS  # noqa: PLW0603
    _STATE = state
    _INSTALLS += 1


def work(key: int) -> tuple:
    start = time.monotonic()
    time.sleep(0.03)
    return key, _STATE, _INSTALLS, os.getpid(), start, time.monotonic()


@unittest.skipUnless(os.environ.get("DATACOOKER_TEST_RAY") == "1", "requires a compute-node Ray session")
class RayExecutionTests(unittest.TestCase):
    def test_cpu_limit_and_state_installation_across_runs(self) -> None:
        for state in ("first", "second"):
            chunks = list(iter_ray_chunks(
                items=list(range(50)), chunk_size=8, n_jobs=2,
                worker_qualname=f"{__name__}.work",
                installer_qualname=f"{__name__}.install",
                state=state, log_label="integration",
            ))
            rows = [row for chunk in chunks for row in chunk]
            self.assertEqual(sorted(row[0] for row in rows), list(range(50)))
            self.assertEqual({row[1] for row in rows}, {state})
            installs: dict[int, set[int]] = {}
            events = []
            for row in rows:
                installs.setdefault(row[3], set()).add(row[2])
                events.extend([(row[4], 1), (row[5], -1)])
            self.assertTrue(all(len(counts) == 1 for counts in installs.values()))
            active = 0
            for _, delta in sorted(events):
                active += delta
                self.assertLessEqual(active, 2)
