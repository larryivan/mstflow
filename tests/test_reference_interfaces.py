import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
BIN = REPO / "bin"


def write_tsv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_script(script, *args, cwd=None, check=True):
    completed = subprocess.run(
        [sys.executable, str(BIN / script), *map(str, args)],
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
    )
    if check and completed.returncode != 0:
        raise AssertionError(
            f"{script} failed with exit {completed.returncode}:\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    return completed


class ReferenceInterfaceTests(unittest.TestCase):
    def test_build_requires_one_to_one_and_reports_kraken_qc(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata.tsv"
            write_tsv(metadata, ["sample_id"], [{"sample_id": "s1"}, {"sample_id": "s2"}])
            reports = []
            for sample in ("s1", "s2"):
                report = root / f"{sample}.kraken2.report.txt"
                report.write_text(
                    "20.00\t200\t200\tU\t0\tunclassified\n"
                    "80.00\t800\t0\tR\t1\troot\n"
                    "60.00\t600\t600\tG\t1301\t  Streptococcus\n"
                )
                reports.append(report)

            run_script("build_feature_table.py", "--metadata", metadata, "--prefix", "kraken2_G",
                       "--taxon-filter", "none", *reports, cwd=root)
            qc = {row["sample_id"]: row for row in read_tsv(root / "kraken2_G_none_sample_qc.tsv")}
            self.assertEqual(qc["s1"]["kraken_classified_fraction"], "0.8")
            self.assertEqual(qc["s1"]["rank_total_reads"], "600")
            self.assertEqual(qc["s1"]["rank_fraction_of_classified"], "0.75")

            incomplete = run_script("build_feature_table.py", "--metadata", metadata, "--prefix", "bad",
                                    reports[0], cwd=root, check=False)
            self.assertNotEqual(incomplete.returncode, 0)
            self.assertIn("missing Kraken report for: s2", incomplete.stderr)

    def test_build_keeps_samples_without_counts_as_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata.tsv"
            write_tsv(
                metadata,
                ["sample_id"],
                [{"sample_id": "blank1"}, {"sample_id": "blank2"}, {"sample_id": "source1"}],
            )
            # blank1 has classified reads but none at the genus; blank2 has none at all.
            (root / "blank1.kraken2.report.txt").write_text(
                "71.11\t64\t64\tU\t0\tunclassified\n"
                "28.89\t26\t0\tR\t1\troot\n"
                "28.89\t26\t26\tF\t1300\t  Streptococcaceae\n"
            )
            (root / "blank2.kraken2.report.txt").write_text("100.00\t12\t12\tU\t0\tunclassified\n")
            (root / "source1.kraken2.report.txt").write_text(
                "20.00\t2\t2\tU\t0\tunclassified\n"
                "80.00\t8\t0\tR\t1\troot\n"
                "80.00\t8\t8\tG\t1301\t  Streptococcus\n"
            )
            run_script("build_feature_table.py", "--metadata", metadata, "--prefix", "kraken2_G",
                       "--taxon-filter", "none", "blank1.kraken2.report.txt", "blank2.kraken2.report.txt",
                       "source1.kraken2.report.txt", cwd=root)

            with open(root / "kraken2_G_counts.tsv", newline="") as handle:
                self.assertEqual(
                    list(csv.reader(handle, delimiter="\t")),
                    [
                        ["#OTU ID", "blank1", "blank2", "source1"],
                        ["taxid_1301__Streptococcus", "0", "0", "8"],
                    ],
                )
            qc = {row["sample_id"]: row for row in read_tsv(root / "kraken2_G_none_sample_qc.tsv")}
            self.assertEqual(qc["blank1"]["kraken_total_reads"], "90")
            self.assertEqual(qc["blank1"]["kraken_classified_reads"], "26")
            self.assertEqual(qc["blank1"]["kraken_unclassified_reads"], "64")
            self.assertEqual(qc["blank1"]["kraken_classified_fraction"], "0.288888888889")
            self.assertEqual(qc["blank1"]["rank_total_reads"], "0")
            self.assertEqual(qc["blank1"]["rank_fraction_of_classified"], "0")
            self.assertEqual(qc["blank1"]["taxon_filter_retained_reads"], "0")
            self.assertEqual(qc["blank1"]["taxon_filter_retained_fraction"], "NA")
            self.assertEqual(qc["blank1"]["raw_nonzero_features"], "0")
            self.assertEqual(qc["blank2"]["kraken_total_reads"], "12")
            self.assertEqual(qc["blank2"]["kraken_classified_fraction"], "0")
            self.assertEqual(qc["blank2"]["rank_fraction_of_classified"], "NA")
            self.assertEqual(qc["source1"]["rank_total_reads"], "8")
            self.assertEqual(qc["source1"]["raw_nonzero_features"], "1")

    def _integer_count_fixture(self, root, genus_reads):
        metadata = root / "metadata.tsv"
        write_tsv(metadata, ["sample_id"], [{"sample_id": "one"}, {"sample_id": "two"}])
        reports = []
        for sample, reads in zip(("one", "two"), genus_reads):
            report = root / f"{sample}.kraken2.report.txt"
            report.write_text(
                "0.00\t0\t0\tU\t0\tunclassified\n"
                "100.00\t10\t0\tR\t1\troot\n"
                f"100.00\t{reads}\t0\tG\t1301\t  Streptococcus\n"
            )
            reports.append(report)
        return metadata, reports

    def test_counts_are_exact_integers_and_reject_fractions(self):
        """Counts must survive as exact integers, and fractions must fail loudly."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            # 8, 8.0 and 8e0 are all the integer 8; a huge count must not become
            # exponent notation, which %.12g would have produced.
            metadata, reports = self._integer_count_fixture(root, ["8.0", "9007199254740993"])
            run_script("build_feature_table.py", "--metadata", metadata, "--prefix", "kraken2_G",
                       "--taxon-filter", "none", *reports, cwd=root)
            with open(root / "kraken2_G_counts.tsv", newline="") as handle:
                self.assertEqual(
                    list(csv.reader(handle, delimiter="\t")),
                    [
                        ["#OTU ID", "one", "two"],
                        ["taxid_1301__Streptococcus", "8", "9007199254740993"],
                    ],
                )

            metadata, reports = self._integer_count_fixture(root, ["8e0", "3"])
            run_script("build_feature_table.py", "--metadata", metadata, "--prefix", "kraken2_G",
                       "--taxon-filter", "none", *reports, cwd=root)
            self.assertEqual(read_tsv(root / "kraken2_G_counts.tsv")[0]["one"], "8")
            # Relative abundance stays floating point.
            self.assertEqual(read_tsv(root / "kraken2_G_relative_abundance.tsv")[0]["one"], "1")

            # "1e999999999" is finite to Decimal but would expand to a
            # billion-digit integer inside a non-interruptible C call. It must
            # be rejected, not materialised.
            for bad in ("8.5", "-1", "nan", "inf", "eight", "1e999999999"):
                metadata, reports = self._integer_count_fixture(root, [bad, "3"])
                failed = run_script("build_feature_table.py", "--metadata", metadata, "--prefix", "bad",
                                    "--taxon-filter", "none", *reports, cwd=root, check=False)
                self.assertNotEqual(failed.returncode, 0, bad)
                # Failures identify the file, the line and the field.
                self.assertIn("one.kraken2.report.txt:3:clade_reads", failed.stderr, bad)

    def test_feature_table_is_invariant_to_report_order(self):
        """Terminal collect() order must not change any published table."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, reports = self._integer_count_fixture(root, ["8", "3"])
            outputs = [
                "kraken2_G_counts.tsv",
                "kraken2_G_relative_abundance.tsv",
                "kraken2_G_features.tsv",
                "kraken2_G_none_counts.tsv",
                "kraken2_G_none_sample_qc.tsv",
            ]
            digests = []
            for order in (reports, list(reversed(reports))):
                run_script("build_feature_table.py", "--metadata", metadata, "--prefix", "kraken2_G",
                           "--taxon-filter", "none", *order, cwd=root)
                digests.append({name: (root / name).read_bytes() for name in outputs})
            self.assertEqual(digests[0], digests[1])
            self.assertEqual(read_tsv(root / "kraken2_G_counts.tsv")[0]["one"], "8")

    def _source_tracker_fixture(self, root):
        metadata = root / "metadata.tsv"
        # The last five columns are reserved by the map; their values must not leak into it.
        fields = [
            "sample_id",
            "unit_id",
            "source_sink",
            "model_source_env",
            "project",
            "host_action",
            "host_group",
            "SourceSink",
            "Env",
            "aggregation_n",
            "aggregation_members",
            "#SampleID",
        ]
        rows = [
            ["bone1", "bone_unit", "source", "bone", "P1", "remove", "human", "bad", "bad", "bad", "bad", "bad"],
            ["bone2", "bone_unit", "source", "bone", "P1", "remove", "human", "bad", "bad", "bad", "bad", "bad"],
            ["soil1", "soil_unit", "source", "soil", "P2", "skip", "NA", "bad", "bad", "bad", "bad", "bad"],
            ["soil2", "soil_unit", "source", "soil", "P2", "skip", "NA", "bad", "bad", "bad", "bad", "bad"],
            ["sink1", "sink_unit", "sink", "NA", "P3", "remove", "human", "bad", "bad", "bad", "bad", "bad"],
            ["blank1", "blank_unit", "sink", "NA", "P3", "skip", "NA", "bad", "bad", "bad", "bad", "bad"],
        ]
        write_tsv(metadata, fields, [dict(zip(fields, row)) for row in rows])
        counts = root / "counts.tsv"
        with open(counts, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["#OTU ID", *[row[0] for row in rows]])
            writer.writerow(["taxid_2__Bacteria", 600, 600, 700, 800, 1500, 2000])
        feature_qc = root / "sample_qc.tsv"
        qc_fields = [
            "sample_id",
            "kraken_total_reads",
            "kraken_classified_reads",
            "kraken_classified_fraction",
            "rank_total_reads",
            "taxon_filter_retained_reads",
            "taxon_filter_retained_fraction",
        ]
        write_tsv(
            feature_qc,
            qc_fields,
            [
                {
                    "sample_id": row[0],
                    "kraken_total_reads": 10000,
                    "kraken_classified_reads": 2000,
                    "kraken_classified_fraction": 0.2,
                    "rank_total_reads": 2500,
                    "taxon_filter_retained_reads": 1200 if row[0].startswith("bone") else 1500,
                    "taxon_filter_retained_fraction": 0.8,
                }
                for row in rows
            ],
        )
        return metadata, counts, feature_qc

    def _prepare(self, root, metadata, counts, feature_qc, name, *extra, check=True):
        return run_script(
            "prepare_attribution_inputs.py",
            "--counts", counts,
            "--metadata", metadata,
            "--feature-sample-qc", feature_qc,
            "--output-counts", root / f"{name}_counts.tsv",
            "--output-map", root / f"{name}_map.tsv",
            "--summary", root / f"{name}_summary.tsv",
            "--sample-qc", root / f"{name}_qc.tsv",
            "--aggregate-by", "unit_id",
            *extra,
            check=check,
        )

    def test_prepare_aggregates_units(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            self._prepare(root, metadata, counts, feature_qc, "st",
                          "--min-source-count", "1000", "--min-sink-count", "1000")
            mapping = {row["#SampleID"]: row for row in read_tsv(root / "st_map.tsv")}
            self.assertEqual(set(mapping), {"bone_unit", "soil_unit", "sink_unit", "blank_unit"})
            self.assertEqual(mapping["soil_unit"]["aggregation_n"], "2")
            self.assertEqual(mapping["soil_unit"]["Env"], "soil")
            self.assertEqual(mapping["soil_unit"]["project"], "P2")
            self.assertEqual(mapping["soil_unit"]["host_action"], "skip")
            qc = {row["input_sample_id"]: row for row in read_tsv(root / "st_qc.tsv")}
            self.assertEqual(qc["soil1"]["analysis_sample_id"], "soil_unit")
            self.assertEqual(qc["soil1"]["included"], "true")
            self.assertEqual(qc["soil1"]["origin"], "samplesheet")
            with open(root / "st_counts.tsv", newline="") as handle:
                header = next(csv.reader(handle, delimiter="\t"))
            self.assertEqual(header, ["#OTU ID", "bone_unit", "soil_unit", "sink_unit", "blank_unit"])

            rows = read_tsv(metadata)
            for row in rows[:2]:
                row["source_sink"] = "NA"
            write_tsv(metadata, list(rows[0]), rows)
            unlabelled = self._prepare(root, metadata, counts, feature_qc, "bad", check=False)
            self.assertIn("neither source nor sink", unlabelled.stderr)

    def test_prepare_names_the_units_it_had_to_leave_out(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            failed = self._prepare(root, metadata, counts, feature_qc, "deep", "--min-sink-count", "5000", check=False)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("no sample to trace is deep enough for source attribution: "
                          "sink_unit (below_min_sink_count:5000), blank_unit (below_min_sink_count:5000)", failed.stderr)
            # A run of sources only is prepared for cross-validation and the panel.
            self._prepare(root, metadata, counts, feature_qc, "sources", "--min-sink-count", "5000", "--no-sinks")
            self.assertEqual({row["SourceSink"] for row in read_tsv(root / "sources_map.tsv")}, {"source"})

    def _panel_fixture(self, root):
        """Two sinks and two sources over three taxa, one only in the sources and one
        only in the sinks; samples to trace come first, as in the samplesheet MSTFlow writes."""
        fields = ["sample_id", "unit_id", "source_sink", "model_source_env", "project", "host_action"]
        rows = [
            ["calculus", "calculus", "sink", "NA", "NA", "remove"],
            ["sediment", "sediment", "sink", "NA", "NA", "skip"],
            ["plaque", "plaque", "source", "oral", "P1", "remove"],
            ["soil", "soil", "source", "soil", "P2", "skip"],
        ]
        taxa = {"taxid_10__Alpha": [0, 0, 30, 0], "taxid_20__Beta": [10, 20, 5, 40], "taxid_5__Gamma": [7, 1, 0, 0]}
        write_tsv(root / "metadata.tsv", fields, [dict(zip(fields, row)) for row in rows])
        with open(root / "counts.tsv", "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["#OTU ID", *[row[0] for row in rows]])
            for taxon, values in taxa.items():
                writer.writerow([taxon, *values])
        write_tsv(root / "qc.tsv", ["sample_id", "kraken_total_reads", "kraken_classified_reads",
                                    "rank_total_reads", "taxon_filter_retained_reads"],
                  [{"sample_id": row[0], "kraken_total_reads": 100, "kraken_classified_reads": 50,
                    "rank_total_reads": 50, "taxon_filter_retained_reads": 50} for row in rows])
        return rows, taxa

    def _prepare_panel(self, root, name, *extra, check=True):
        (root / name).mkdir(exist_ok=True)
        return self._prepare(root, root / "metadata.tsv", root / "counts.tsv", root / "qc.tsv", name,
                             "--min-taxon-filter-retained-reads", "0", "--write-panel", root / name, *extra,
                             check=check)

    def test_a_panel_reproduces_the_run_that_made_it(self):
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "first", Path(directory) / "second"
            first.mkdir()
            second.mkdir()
            rows, taxa = self._panel_fixture(first)
            self._prepare_panel(first, "full")
            self.assertEqual(read_tsv(first / "full/sources.tsv"), [
                {"source_id": "plaque", "environment": "oral", "project": "P1", "libraries": "plaque"},
                {"source_id": "soil", "environment": "soil", "project": "P2", "libraries": "soil"},
            ])
            # Gamma is in the sinks only, so it is not part of the panel.
            self.assertEqual([row["#OTU ID"] for row in read_tsv(first / "full/counts.tsv")],
                             ["taxid_10__Alpha", "taxid_20__Beta"])

            # The same sinks against the panel instead of the sources' reads.
            sinks = [row for row in rows if row[2] == "sink"]
            fields = ["sample_id", "unit_id", "source_sink", "model_source_env", "project", "host_action"]
            write_tsv(second / "metadata.tsv", fields, [dict(zip(fields, row)) for row in sinks])
            with open(second / "counts.tsv", "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["#OTU ID", "calculus", "sediment"])
                for taxon, values in taxa.items():
                    if any(values[:2]):
                        writer.writerow([taxon, *values[:2]])
            qc = [row for row in read_tsv(first / "qc.tsv") if row["sample_id"] in ("calculus", "sediment")]
            write_tsv(second / "qc.tsv", list(qc[0]), qc)
            self._prepare_panel(second, "reused", "--panel-counts", first / "full/counts.tsv",
                                "--panel-sources", first / "full/sources.tsv")

            self.assertEqual((first / "full_counts.tsv").read_bytes(), (second / "reused_counts.tsv").read_bytes())
            labels = ("#SampleID", "SourceSink", "Env", "project", "aggregation_members")
            self.assertEqual([tuple(row[field] for field in labels) for row in read_tsv(first / "full_map.tsv")],
                             [tuple(row[field] for field in labels) for row in read_tsv(second / "reused_map.tsv")])
            qc = {row["input_sample_id"]: row for row in read_tsv(second / "reused_qc.tsv")}
            self.assertEqual((qc["plaque"]["origin"], qc["plaque"]["included"]), ("panel", "true"))
            self.assertEqual((qc["calculus"]["origin"], qc["calculus"]["included"]), ("samplesheet", "true"))

    def test_panel_sources_are_depth_checked_and_named_apart(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self._panel_fixture(root)
            self._prepare_panel(root, "full")
            panel = ["--panel-counts", root / "full/counts.tsv", "--panel-sources", root / "full/sources.tsv"]

            clash = self._prepare_panel(root, "clash", *panel, check=False)
            self.assertIn("the panel has source(s) named like samples of this run: plaque, soil", clash.stderr)

            rows = read_tsv(root / "full/sources.tsv")
            for row in rows:
                row["source_id"] = "panel_" + row["source_id"]
                row["environment"] = " " + row["environment"].upper()
            write_tsv(root / "full/sources.tsv", list(rows[0]), rows)
            with open(root / "full/counts.tsv") as handle:
                text = handle.read().replace("\tplaque\tsoil", "\tpanel_plaque\tpanel_soil", 1)
            (root / "full/counts.tsv").write_text(text)
            self._prepare_panel(root, "deep", *panel, "--min-source-count", "40")
            qc = {row["input_sample_id"]: row for row in read_tsv(root / "deep_qc.tsv")}
            self.assertEqual(qc["panel_plaque"]["exclusion_reason"], "below_min_source_count:40")
            self.assertEqual(qc["panel_soil"]["included"], "true")
            mapping = {row["#SampleID"]: row for row in read_tsv(root / "deep_map.tsv")}
            self.assertEqual(mapping["panel_soil"]["Env"], "soil")
            # This run's sources and the panel's together make the next panel, without
            # those too shallow to use.
            self.assertEqual([row["source_id"] for row in read_tsv(root / "deep/sources.tsv")],
                             ["soil", "panel_soil"])

    def test_prepare_reports_undefined_qc_as_na_and_rejects_fractional_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)

            # A sample with no counts at the rank has no defined retained
            # fraction. That must be reported as NA, not as an empty cell.
            qc_rows = read_tsv(feature_qc)
            for row in qc_rows:
                if row["sample_id"] == "sink1":
                    row["rank_total_reads"] = "0"
                    row["taxon_filter_retained_reads"] = "0"
                    row["taxon_filter_retained_fraction"] = "NA"
            write_tsv(feature_qc, list(qc_rows[0]), qc_rows)
            # Promote the blank to the surviving sink so the run still completes
            # while the zero-estimate sink is excluded on undefined QC.
            metadata_rows = read_tsv(metadata)
            for row in metadata_rows:
                if row["sample_id"] == "blank1":
                    row["source_sink"] = "sink"
            write_tsv(metadata, list(metadata_rows[0]), metadata_rows)
            self._prepare(root, metadata, counts, feature_qc, "st",
                          "--min-taxon-filter-retained-reads", "0",
                          "--min-taxon-filter-retained-fraction", "0")
            qc = {row["input_sample_id"]: row for row in read_tsv(root / "st_qc.tsv")}
            self.assertEqual(qc["sink1"]["taxon_filter_retained_fraction"], "NA")
            self.assertEqual(qc["sink1"]["rank_total_reads"], "0")
            # Counts stay exact integers, never 1500.0 or 1.5e+03.
            self.assertEqual(qc["soil1"]["taxon_filter_retained_reads"], "3000")
            self.assertEqual(qc["soil1"]["input_total_count"], "700")

            fractional = root / "fractional_counts.tsv"
            with open(fractional, "w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    ["#OTU ID", "bone1", "bone2", "soil1", "soil2", "sink1", "blank1"]
                )
                writer.writerow(["taxid_2__Bacteria", 600, 600, "700.5", 800, 1500, 2000])
            failed = self._prepare(root, metadata, fractional, feature_qc, "bad", check=False)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("fractional_counts.tsv:2", failed.stderr)
            self.assertIn("taxid_2__Bacteria", failed.stderr)
            self.assertIn("soil1", failed.stderr)
            self.assertIn("whole read count", failed.stderr)

    def test_prepare_validates_and_filters_complete_units(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            rows = read_tsv(metadata)
            for row in rows:
                if row["sample_id"] == "bone1":
                    row["model_source_env"] = "soil"
            write_tsv(metadata, list(rows[0]), rows)
            conflicting = self._prepare(root, metadata, counts, feature_qc, "conflict", check=False)
            self.assertNotEqual(conflicting.returncode, 0)
            self.assertIn("bone_unit", conflicting.stderr)
            self.assertIn("model_source_env", conflicting.stderr)

            rows = read_tsv(metadata)
            for row in rows:
                if row["sample_id"] == "bone1":
                    row["model_source_env"] = "bone"
            write_tsv(metadata, list(rows[0]), rows)
            qc_rows = read_tsv(feature_qc)
            for row in qc_rows:
                if row["sample_id"] == "bone1":
                    row["kraken_classified_reads"] = "0"
                    row["kraken_classified_fraction"] = "0"
            write_tsv(feature_qc, list(qc_rows[0]), qc_rows)
            self._prepare(root, metadata, counts, feature_qc, "filtered",
                          "--min-kraken-classified-fraction", "0.15")
            mapping = {row["#SampleID"] for row in read_tsv(root / "filtered_map.tsv")}
            self.assertNotIn("bone_unit", mapping)
            excluded = {
                row["input_sample_id"]: row["exclusion_reason"]
                for row in read_tsv(root / "filtered_qc.tsv")
                if row["included"] == "false"
            }
            self.assertTrue(excluded["bone1"].startswith("below_min_kraken"))
            self.assertEqual(excluded["bone1"], excluded["bone2"])

    def test_prepare_rejects_reserved_source_environment_names(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            rows = read_tsv(metadata)
            for row in rows:
                if row["sample_id"].startswith("soil"):
                    row["model_source_env"] = "Unknown"
            write_tsv(metadata, list(rows[0]), rows)
            failed = self._prepare(root, metadata, counts, feature_qc, "reserved", check=False)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("reserved result columns Unknown or SampleID", failed.stderr)

    def test_all_taxon_scope_uses_rank_yield_not_filtered_yield(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            qc_rows = read_tsv(feature_qc)
            for row in qc_rows:
                row["taxon_filter_retained_reads"] = "0"
                row["taxon_filter_retained_fraction"] = "0"
            write_tsv(feature_qc, list(qc_rows[0]), qc_rows)
            self._prepare(root, metadata, counts, feature_qc, "all", "--taxon-scope", "all")
            mapping = {row["#SampleID"] for row in read_tsv(root / "all_map.tsv")}
            self.assertEqual(mapping, {"bone_unit", "soil_unit", "sink_unit", "blank_unit"})

    def test_seeded_sourcetracker_wrapper_validates_and_forwards_seed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "numpy.py").write_text(
                "class Random:\n"
                "    def seed(self, value):\n"
                "        print(f'numpy_seed={value}')\n"
                "random = Random()\n"
            )
            cli_dir = root / "sourcetracker" / "_cli"
            cli_dir.mkdir(parents=True)
            (root / "sourcetracker" / "__init__.py").write_text("")
            (cli_dir / "__init__.py").write_text(
                "class CLI:\n"
                "    def main(self, args, prog_name, standalone_mode):\n"
                "        print('forwarded=' + ' '.join(args))\n"
                "cli = CLI()\n"
            )
            env = os.environ.copy()
            env["PYTHONPATH"] = str(root)
            env["PYTHONHASHSEED"] = "17"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(BIN / "run_sourcetracker2_seeded.py"),
                    "--seed",
                    "17",
                    "gibbs",
                    "--jobs",
                    "1",
                ],
                text=True,
                capture_output=True,
                env=env,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertIn("numpy_seed=17", completed.stdout)
            self.assertIn("forwarded=gibbs --jobs 1", completed.stdout)

            env["PYTHONHASHSEED"] = "18"
            mismatch = subprocess.run(
                [
                    sys.executable,
                    str(BIN / "run_sourcetracker2_seeded.py"),
                    "--seed",
                    "17",
                    "gibbs",
                ],
                text=True,
                capture_output=True,
                env=env,
                check=False,
            )
            self.assertNotEqual(mismatch.returncode, 0)
            self.assertIn("PYTHONHASHSEED must be exported", mismatch.stderr)



if __name__ == "__main__":
    unittest.main()
