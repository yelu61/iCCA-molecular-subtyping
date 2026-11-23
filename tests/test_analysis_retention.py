"""Guard the selected upstream implementation and analysis entry inventory."""
import csv
import hashlib
import json
from pathlib import Path
import unittest

REPO = Path(__file__).resolve().parents[1]


class AnalysisRetentionTest(unittest.TestCase):
    def test_selected_implementation_matches_its_manifest(self):
        with (REPO / "IMPLEMENTATION_MANIFEST.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        self.assertTrue(rows)
        for row in rows:
            with self.subTest(path=row["path"]):
                path = REPO / "_implementation/project" / row["path"]
                self.assertTrue(path.is_file(), "Required analysis source is missing")
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(),
                                 row["packaged_sha256"])

    def test_every_analysis_entry_has_source(self):
        catalog = json.loads((REPO / "config/analysis_stages.json").read_text())
        self.assertTrue(catalog)
        for name, stage in catalog.items():
            with self.subTest(stage=name):
                self.assertTrue((REPO / "_implementation/project" / stage["entry"]).is_file())


if __name__ == "__main__":
    unittest.main()
