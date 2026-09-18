"""Contracts of the source attribution scripts: the shared result table, the call
rule, the SourceTracker2 converter, cross-validation folds, the summary and the
evaluation. Plots are skipped (--no-plot), so only the standard library is needed."""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

import mstflow_attribution as attribution  # noqa: E402


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(script, *args, cwd, check=True):
    completed = subprocess.run([sys.executable, str(BIN / script), *map(str, args)], cwd=cwd,
                               capture_output=True, text=True)
    if check and completed.returncode != 0:
        raise AssertionError(f"{script} failed:\n{completed.stdout}\n{completed.stderr}")
    return completed


MAP_HEADER = ["#SampleID", "SourceSink", "Env", "project"]


class DecideTests(unittest.TestCase):
    def decide(self, **estimates):
        return attribution.decide(estimates, max_unknown=0.5, min_top_known=0.5, min_margin=0.05)

    def test_clear_sink_is_called(self):
        result = self.decide(oral=0.8, gut=0.1, Unknown=0.1)
        self.assertEqual(result["call"], "oral")
        self.assertEqual(result["reasons"], "NA")
        self.assertAlmostEqual(result["margin"], 0.7)

    def test_every_failed_rule_is_reported(self):
        result = self.decide(oral=0.2, gut=0.18, Unknown=0.62)
        self.assertEqual(result["call"], "ABSTAIN")
        self.assertEqual(result["reasons"], "unknown>=0.5;top_known<0.5;margin<0.05")

    def test_tie_abstains(self):
        result = self.decide(oral=0.5, gut=0.5, Unknown=0.0)
        self.assertEqual(result["call"], "ABSTAIN")
        self.assertIn("tie", result["reasons"])


class ResultTableTests(unittest.TestCase):
    def test_estimates_are_rescaled_and_ordered(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r.tsv"
            attribution.write_results(path, "main", "feast", "proportion",
                                      {"s1": {"Unknown": (0.2, None), "soil": (0.4004, 0.1), "gut": (0.4, None)}})
            rows = read_tsv(path)
            self.assertEqual([row["source"] for row in rows], ["gut", "soil", "Unknown"])
            self.assertAlmostEqual(sum(float(row["estimate"]) for row in rows), 1.0, places=5)
            self.assertEqual(rows[0]["sd"], "NA")

    def test_bad_tables_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "r.tsv"
            cases = {
                "no Unknown": {"s1": {"soil": (1.0, None)}},
                "sum to 0.6": {"s1": {"soil": (0.4, None), "Unknown": (0.2, None)}},
                "different sources": {"s1": {"soil": (1.0, None), "Unknown": (0.0, None)},
                                      "s2": {"gut": (1.0, None), "Unknown": (0.0, None)}},
            }
            for message, estimates in cases.items():
                with self.subTest(message=message), self.assertRaises(SystemExit) as raised:
                    attribution.write_results(path, "main", "feast", "proportion", estimates)
                self.assertIn(message, str(raised.exception))


class SourceTracker2ConverterTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        write_tsv(self.root / "map.tsv", MAP_HEADER,
                  [["a", "source", "soil", "P1"], ["b", "source", "gut", "P2"], ["s1", "sink", "NA", "P3"]])

    def tearDown(self):
        self.directory.cleanup()

    def convert(self, means, stds="SampleID\tgut\tsoil\tUnknown\ns1\t0.01\t0.02\t0.03\n", check=True):
        (self.root / "means.txt").write_text(means)
        (self.root / "stds.txt").write_text(stds)
        return run("sourcetracker2_results.py", "--proportions", "means.txt", "--stds", "stds.txt",
                   "--map", "map.tsv", "--analysis", "main", "--output", "out.tsv", cwd=self.root, check=check)

    def test_converts_means_and_deviations(self):
        self.convert("SampleID\tgut\tsoil\tUnknown\ns1\t0.1\t0.7\t0.2\n")
        rows = {row["source"]: row for row in read_tsv(self.root / "out.tsv")}
        self.assertEqual(rows["soil"]["estimate"], "0.7")
        self.assertEqual(rows["soil"]["sd"], "0.02")
        self.assertEqual({row["method"] for row in rows.values()}, {"sourcetracker2"})

    def test_rejects_duplicate_columns_and_missing_sinks(self):
        duplicate = self.convert("SampleID\tsoil\tsoil\tUnknown\ns1\t0.4\t0.4\t0.2\n", check=False)
        self.assertNotEqual(duplicate.returncode, 0)
        self.assertIn("duplicate columns: soil", duplicate.stderr)
        missing = self.convert("SampleID\tgut\tsoil\tUnknown\n", check=False)
        self.assertNotEqual(missing.returncode, 0)
        self.assertIn("exactly the sinks", missing.stderr)


class FoldTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        units = [["oral1", "source", "oral", "P1"], ["oral2", "source", "oral", "P2"],
                 ["soil1", "source", "soil", "P3"], ["soil2", "source", "soil", "P3"],
                 ["gut1", "source", "gut", ""], ["sink1", "sink", "NA", "P9"]]
        write_tsv(self.root / "map.tsv", MAP_HEADER, units)
        write_tsv(self.root / "counts.tsv", ["#OTU ID", *[u[0] for u in units]],
                  [["t1", 10, 10, 0, 0, 5, 3], ["t2", 0, 0, 10, 10, 5, 3], ["t3", 0, 0, 0, 0, 0, 0]])

    def tearDown(self):
        self.directory.cleanup()

    def split(self, *extra, check=True):
        return run("split_attribution_folds.py", "--counts", "counts.tsv", "--map", "map.tsv",
                   "--group-column", "project", "--out-dir", "folds", "--summary", "folds.tsv", *extra,
                   cwd=self.root, check=check)

    def test_whole_groups_are_held_out_as_sinks_of_known_origin(self):
        self.split()
        summary = {row["group"]: row for row in read_tsv(self.root / "folds.tsv")}
        self.assertEqual(set(summary), {"P1", "P2", "P3", "gut1"})
        self.assertEqual(summary["P3"]["held_out_units"], "soil1,soil2")
        self.assertEqual(summary["P3"]["envs_absent_from_training"], "soil")
        fold = read_tsv(self.root / "folds" / "cv_P1.map.tsv")
        held_out = [row for row in fold if row["SourceSink"] == "sink"]
        self.assertEqual([(row["#SampleID"], row["Env"], row["expected_env"]) for row in held_out],
                         [("oral1", "NA", "oral")])
        self.assertNotIn("sink1", [row["#SampleID"] for row in fold])
        with open(self.root / "folds" / "cv_P1.counts.tsv") as handle:
            taxa = [line.split("\t")[0] for line in handle][1:]
        self.assertEqual(taxa, ["t1", "t2"])

    def test_shallow_units_are_skipped_and_one_group_is_refused(self):
        self.split("--min-sink-count", "15")
        notes = {row["group"]: row["note"] for row in read_tsv(self.root / "folds.tsv")}
        self.assertTrue(notes["gut1"].startswith("skipped"))
        write_tsv(self.root / "map.tsv", MAP_HEADER,
                  [["oral1", "source", "oral", "P1"], ["oral2", "source", "oral", "P1"],
                   ["soil1", "source", "soil", "P1"], ["soil2", "source", "soil", "P1"],
                   ["gut1", "source", "gut", "P1"], ["sink1", "sink", "NA", "P1"]])
        refused = self.split(check=False)
        self.assertNotEqual(refused.returncode, 0)
        self.assertIn("at least two 'project' groups", refused.stderr)


class SummaryAndEvaluationTests(unittest.TestCase):
    def test_methods_are_compared_per_sink(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_tsv(root / "map.tsv", MAP_HEADER,
                      [["a", "source", "oral", "P1"], ["b", "source", "soil", "P2"],
                       ["s1", "sink", "NA", "P3"], ["s2", "sink", "NA", "P3"]])
            for method, s1, s2 in (("feast", (0.9, 0.05, 0.05), (0.1, 0.1, 0.8)),
                                   ("sourcetracker2", (0.8, 0.1, 0.1), (0.2, 0.7, 0.1))):
                attribution.write_results(root / f"{method}.tsv", "main", method, "proportion", {
                    "s1": dict(zip(["oral", "soil", "Unknown"], [(v, None) for v in s1])),
                    "s2": dict(zip(["oral", "soil", "Unknown"], [(v, None) for v in s2])),
                })
            run("summarize_attribution.py", "--no-plot", "--map", "map.tsv", "--max-unknown", "0.5",
                "--min-top-known", "0.5", "--min-margin", "0.05", "feast.tsv", "sourcetracker2.tsv", cwd=root)
            comparison = {row["sink"]: row for row in read_tsv(root / "comparison.tsv")}
            self.assertEqual(comparison["s1"]["agreement"], "agree:oral")
            self.assertEqual(comparison["s2"]["agreement"], "one_method")
            self.assertEqual(comparison["s2"]["feast"], "ABSTAIN")
            self.assertEqual(list(comparison["s1"])[1:3], ["sourcetracker2", "feast"])
            self.assertEqual(len(read_tsv(root / "assignments.tsv")), 4)

            # With one method, every call is that method's own.
            run("summarize_attribution.py", "--no-plot", "--map", "map.tsv", "--max-unknown", "0.5",
                "--min-top-known", "0.5", "--min-margin", "0.05", "feast.tsv", cwd=root)
            self.assertEqual(read_tsv(root / "comparison.tsv")[0]["agreement"], "agree:oral")
            write_tsv(root / "map.tsv", MAP_HEADER,
                      [["a", "source", "oral", "P1"], ["s1", "sink", "NA", "P3"], ["s3", "sink", "NA", "P3"]])
            incomplete = run("summarize_attribution.py", "--no-plot", "--map", "map.tsv", "--max-unknown", "0.5",
                             "--min-top-known", "0.5", "--min-margin", "0.05", "feast.tsv", cwd=root, check=False)
            self.assertNotEqual(incomplete.returncode, 0)
            self.assertIn("feast/s3", incomplete.stderr)

    def test_evaluation_scores_calls_and_abstentions(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write_tsv(root / "folds.tsv",
                      ["fold", "group", "held_out_units", "training_envs", "envs_absent_from_training", "note"],
                      [["cv_A", "A", "o1,o2,s1", "oral,soil", "", "ok"], ["cv_B", "B", "g1", "oral,soil", "gut", "ok"]])
            header = ["#SampleID", "SourceSink", "Env", "expected_env"]
            write_tsv(root / "cv_A.map.tsv", header, [["o1", "sink", "NA", "oral"], ["o2", "sink", "NA", "oral"],
                                                      ["s1", "sink", "NA", "soil"], ["x", "source", "oral", "NA"],
                                                      ["y", "source", "soil", "NA"]])
            write_tsv(root / "cv_B.map.tsv", header, [["g1", "sink", "NA", "gut"], ["x", "source", "oral", "NA"],
                                                      ["y", "source", "soil", "NA"]])
            shares = {"o1": (0.9, 0.0, 0.1), "o2": (0.1, 0.8, 0.1), "s1": (0.1, 0.1, 0.8)}
            attribution.write_results(root / "a.tsv", "cv_A", "feast", "proportion", {
                unit: dict(zip(["oral", "soil", "Unknown"], [(v, None) for v in values]))
                for unit, values in shares.items()})
            attribution.write_results(root / "b.tsv", "cv_B", "feast", "proportion", {
                "g1": {"oral": (0.1, None), "soil": (0.1, None), "Unknown": (0.8, None)}})
            run("evaluate_attribution.py", "--no-plot", "--folds", "folds.tsv", "--maps", "cv_A.map.tsv",
                "cv_B.map.tsv", "--results", "a.tsv", "b.tsv", "--max-unknown", "0.5", "--min-top-known", "0.5",
                "--min-margin", "0.05", cwd=root)
            row = read_tsv(root / "evaluation.tsv")[0]
            # o1 right, o2 wrong (called soil), s1 abstained; g1's environment was unseen and it abstained.
            self.assertEqual((row["held_out"], row["called"], row["coverage"], row["accuracy"]), ("3", "2", "0.6667", "0.5"))
            self.assertEqual(row["correct_overall"], "0.3333")
            self.assertEqual(row["balanced_accuracy"], "0.25")
            self.assertEqual((row["held_out_unseen_env"], row["abstained_on_unseen_env"]), ("1", "1"))


if __name__ == "__main__":
    unittest.main()
