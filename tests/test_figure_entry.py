"""Regression checks for the portable public plotting notebook entry.

Real R rendering and input immutability are recorded separately in the freeze
verification. These tests exercise path resolution and subprocess contracts.
"""
import ast
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch


REPO = Path(__file__).resolve().parents[1]


def entry_functions():
    notebook = json.loads((REPO / "pipelines/8-visualization.ipynb").read_text())
    setup = next("".join(c["source"]) for c in notebook["cells"]
                 if c["cell_type"] == "code"
                 and "def locate_figure_layout(" in "".join(c["source"]))
    tree = ast.parse(setup)
    tree.body = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    namespace = dict(Path=Path, os=os, shutil=shutil, subprocess=subprocess)
    exec(compile(tree, "figure_entry", "exec"), namespace)
    return namespace


class FigureEntryTest(unittest.TestCase):
    def setUp(self):
        self.entry = entry_functions()

    def test_clone_and_notebook_directory_resolve_to_same_inputs(self):
        expected = (REPO, REPO / "figures", REPO / "figures/data",
                    REPO / "figures/outputs")
        for start in (REPO, REPO / "pipelines"):
            self.assertEqual(self.entry["locate_figure_layout"](start), expected)

    def test_explicit_runtime_prefers_binary_over_redirecting_wrapper(self):
        with tempfile.TemporaryDirectory() as folder:
            home = Path(folder).resolve()
            (home / "bin/exec").mkdir(parents=True)
            (home / "bin/R").touch()
            (home / "bin/exec/R").touch()
            with patch.dict(os.environ, {"ICCA_R_HOME": str(home)}):
                command, env = self.entry["locate_r_runtime"]()
            self.assertEqual(command, home / "bin/exec/R")
            self.assertEqual(env["R_HOME"], str(home))

    def test_explicit_inputs_and_missing_pdf_recovery(self):
        with tempfile.TemporaryDirectory() as folder:
            output = Path(folder)
            self.entry.update(ROOT=REPO, SCRIPT_DIR=REPO / "figures",
                              DATA_ROOT=REPO / "figures/data", OUTPUT_ROOT=output,
                              R_EXEC=Path("/test/R"), R_ENV={})

            def render(command, **kwargs):
                self.assertEqual(kwargs["cwd"], REPO)
                self.assertEqual(command[command.index("--data-dir") + 1],
                                 str(REPO / "figures/data/figure1"))
                self.assertEqual(command[command.index("-f") + 1],
                                 str(REPO / "figures/plot_figure1.R"))
                dest = Path(command[command.index("--out-dir") + 1])
                dest.mkdir(parents=True, exist_ok=True)
                for suffix in ("pdf", "png"):
                    (dest / f"Figure_1_composite.{suffix}").write_bytes(b"test fixture")
                return subprocess.CompletedProcess(command, 0)

            with patch.object(subprocess, "run", side_effect=render) as run, \
                    contextlib.redirect_stdout(io.StringIO()):
                png, pdf = self.entry["run_figure_script"](1)
                self.entry["run_figure_script"](1)
                self.assertEqual(run.call_count, 1)
                pdf.unlink()
                self.entry["run_figure_script"](1)
                self.assertEqual(run.call_count, 2)
                self.assertTrue(png.is_file() and pdf.is_file())

    def test_renderer_success_requires_both_nonempty_outputs(self):
        with tempfile.TemporaryDirectory() as folder:
            self.entry.update(ROOT=REPO, SCRIPT_DIR=REPO / "figures",
                              DATA_ROOT=REPO / "figures/data", OUTPUT_ROOT=Path(folder),
                              R_EXEC=Path("/test/R"), R_ENV={})
            with patch.object(subprocess, "run", return_value=None), \
                    contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(RuntimeError, "both non-empty"):
                    self.entry["run_figure_script"](1)


if __name__ == "__main__":
    unittest.main()
