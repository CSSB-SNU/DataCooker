"""Workflow CLI: runtime_environment is applied to os.environ, never passed to a recipe."""
import os

from click.testing import CliRunner

from datacooker.cli import workflow


def test_run_applies_runtime_environment_and_hides_it_from_the_recipe(tmp_path, monkeypatch):
    recipe = tmp_path / "recipe.py"
    recipe.write_text("")
    config = tmp_path / "config.yaml"
    config.write_text(
        f"recipe_path: {recipe}\n"
        f"output_data_path: {tmp_path / 'out.tsv'}\n"
        "inputs: {a: 1}\n"
        "runtime_environment: {DC_TEST_REFERENCE: pinned}\n",
    )
    seen = {}

    def fake_run_recipe(recipe, **kwargs):
        seen.update(kwargs)
        return {}

    monkeypatch.delenv("DC_TEST_REFERENCE", raising=False)
    monkeypatch.setattr(workflow, "run_recipe", fake_run_recipe)
    result = CliRunner().invoke(workflow.cli, ["run", str(config)])
    assert result.exit_code == 0, result.output
    assert "runtime_environment" not in seen
    assert os.environ["DC_TEST_REFERENCE"] == "pinned"
