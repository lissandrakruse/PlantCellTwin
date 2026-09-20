import ast
import csv
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
EXPECTED = json.loads((ROOT / "tests" / "expected_metrics.json").read_text())

def rows(name):
    with (RESULTS / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))

class ReleaseSmokeTests(unittest.TestCase):
    def test_required_files(self):
        required = ["README.md", "LICENSE", "CITATION.cff", "environment.yml",
                    "requirements.txt", "index.html", "app.js",
                    "results/bioconductor-python-concordance.csv",
                    "results/osd251-gravity-dose-validation.csv"]
        for name in required:
            self.assertTrue((ROOT / name).is_file(), name)

    def test_reference_concordance(self):
        root = next(r for r in rows("bioconductor-python-concordance.csv")
                    if r.get("tissue", "").lower() == "root")
        values = {k.lower(): v for k, v in root.items()}
        self.assertEqual(int(float(values["bioc_fdr_005"])), EXPECTED["root_r_fdr_005"])
        self.assertEqual(int(float(values["python_fdr_005"])), EXPECTED["root_python_fdr_005"])
        self.assertGreaterEqual(float(values["effect_spearman"]), EXPECTED["root_effect_rho_min"])
        self.assertEqual(int(float(values["top300_overlap"])), EXPECTED["root_top300_overlap"])

    def test_osd251_panel(self):
        panel = [r for r in rows("osd251-gravity-dose-validation.csv")
                 if "root" in r.get("panel", "").lower()]
        self.assertEqual(len(panel), EXPECTED["osd251_root_panel_size"])
        self.assertEqual(sum(float(r["gravity_spearman_rho"]) > 0 for r in panel),
                         EXPECTED["osd251_root_panel_positive_trends"])

    def test_no_machine_paths_in_interface(self):
        for name in ("index.html", "app.js", "styles.css", "plant.css"):
            self.assertNotIn("/workspace/scratch/", (ROOT / name).read_text(), name)

    def test_analysis_python_syntax(self):
        for path in (ROOT / "analysis").glob("*.py"):
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

if __name__ == "__main__":
    unittest.main()
