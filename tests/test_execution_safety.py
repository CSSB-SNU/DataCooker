import unittest

from datacooker._ray import Admission
from datacooker.cache import ExecutionContext
from datacooker.lmdb.index import ItemSize
from datacooker.planning import NodeSpec, plan, shard_keys


class ExecutionSafetyTests(unittest.TestCase):
    def test_context_lists_declared_names_without_traversing_payload(self) -> None:
        context = ExecutionContext(lambda name: name.split("."))
        payload = {}
        payload["cycle"] = payload
        context.add_data("input.structure", payload)
        context.add_data("input.count", 2)
        self.assertEqual(context.keys(), ["input.structure", "input.count"])
        self.assertIs(context.snapshot()["input.structure"], payload)
        self.assertEqual(context["input.count"], 2)


    def test_planner_distributes_once_with_exact_key_coverage(self) -> None:
        sizes = [ItemSize(str(i), 100-i) for i in range(12)]
        shards = shard_keys(plan(sizes, NodeSpec(n_nodes=3)).tiers[0])
        self.assertEqual(shards, [["0", "5", "6", "11"], ["1", "4", "7", "10"], ["2", "3", "8", "9"]])
        self.assertEqual(sorted(key for shard in shards for key in shard), sorted(str(i) for i in range(12)))


    def test_admission_ramps_and_recovers_without_permanent_ratchet(self) -> None:
        policy = Admission(8)
        self.assertEqual(policy.target, 1)
        for _ in range(40):
            policy.observe(0.2, 1)
        self.assertEqual(policy.target, 8)
        policy.observe(0.8, 0)
        self.assertEqual(policy.target, 4)
        self.assertFalse(policy.admits(1, 0.7))
        for _ in range(40):
            policy.observe(0.2, 1)
        self.assertEqual(policy.target, 8)
        self.assertTrue(policy.admits(0, 0.9))

    def test_pipeline_merges_only_planned_shards(self) -> None:
        import tempfile
        from pathlib import Path

        from datacooker.executors.slurm import SlurmExecutor
        from datacooker.pipeline import PipelineConfig, run
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            cfg = PipelineConfig(op="rebuild", source=root / "source.lmdb",
                                 target=root / "result.lmdb", engine_config=root / "engine.yaml",
                                 node=NodeSpec(n_nodes=2))
            executor = SlurmExecutor(workdir=root, repo=root, submit=False)
            run(cfg, executor, keyfiles_dir=root / "items",
                sizes=[ItemSize("a", 10), ItemSize("b", 20)])
            paths = (root / "items/merge_shards.txt").read_text().splitlines()
            self.assertEqual(len(paths), 2)
            self.assertTrue(all("result_all_shard" in p for p in paths))
            command = (root / "result_publish.cmds").read_text()
            self.assertIn("--shard-list", command)
            self.assertNotIn("*", command)

    def test_nested_external_dependency_validates_without_payload_traversal(self) -> None:
        from datacooker import RecipeBook, execute
        recipe = RecipeBook()
        recipe.step(outputs=(("result", str),), instruction=str, args=[("_entry.id", str)])
        cycle = {}
        cycle["self"] = cycle
        result = execute(recipe, {"_entry": {"id": "5odi"}, "payload": cycle},
                         key_transform=lambda name: name.split("."))
        self.assertEqual(result, {"result": "5odi"})


def test_parallel_run_body_resolves_its_helpers_when_run_as_module(tmp_path):
    # Run through `python -m`, as SLURM jobs do: a helper defined below the module's
    # `if __name__ == "__main__": main()` line is not yet bound there (NameError on every node).
    import subprocess
    import sys

    config = tmp_path / "empty.yaml"
    config.write_text("{}\n")
    result = subprocess.run(
        [sys.executable, "-m", "datacooker.cli.workflow", "parallel-run", str(config)],
        capture_output=True, text=True, check=False,
    )
    assert "NameError" not in result.stderr


def test_node_keys_in_config_do_not_reach_the_runner_twice():
    from datacooker.cli._common import resolve_optional_int

    config = {"node_count": 8, "n_jobs": -1}
    assert resolve_optional_int(3, config, ("node_count", "world_size")) == 3
    assert "node_count" not in config
    config = {"world_size": 4}
    assert resolve_optional_int(None, config, ("node_count", "world_size")) == 4
    assert config == {}
