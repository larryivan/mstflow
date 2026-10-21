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
    def test_validate_host_action_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fields = [
                "sample_id",
                "fastq_1",
                "fastq_2",
                "layout",
                "age_group",
                "host_action",
                "host_group",
                "host_ref",
            ]
            base = {
                "sample_id": "sample_1",
                "fastq_1": "reads.fq.gz",
                "fastq_2": "NA",
                "layout": "SE",
                "age_group": "ancient",
                "host_action": "skip",
                "host_group": "none",
                "host_ref": "NA",
            }
            metadata = root / "metadata.tsv"
            output = root / "validated.tsv"

            write_tsv(metadata, fields, [base])
            run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                "--require-explicit-host-action",
            )
            normalized = read_tsv(output)[0]
            self.assertEqual(normalized["host_action"], "skip")
            self.assertEqual(normalized["host_group"], "NA")
            self.assertEqual(normalized["host_ref"], "NA")

            missing_action = dict(base, host_action="NA")
            write_tsv(metadata, fields, [missing_action])
            missing = run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                "--require-explicit-host-action",
                check=False,
            )
            self.assertNotEqual(missing.returncode, 0)
            self.assertIn("host_action must be explicitly set", missing.stderr)

            skip_with_reference = dict(base, host_ref="default")
            write_tsv(metadata, fields, [skip_with_reference])
            conflicted = run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                check=False,
            )
            self.assertNotEqual(conflicted.returncode, 0)
            self.assertIn("host_action=skip conflicts with host_ref=default", conflicted.stderr)

            # Without a host_ref column, skip needs nothing else.
            no_ref_fields = [f for f in fields if f != "host_ref"]
            no_ref = {k: v for k, v in base.items() if k != "host_ref"}
            write_tsv(metadata, no_ref_fields, [no_ref])
            run_script("validate_metadata.py", "--input", metadata, "--output", output)
            self.assertEqual(read_tsv(output)[0]["host_ref"], "NA")

            remove = dict(
                base,
                host_action="remove",
                host_group="human",
                host_ref="default",
            )
            write_tsv(metadata, fields, [remove])
            run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                "--require-explicit-host-action",
            )
            self.assertEqual(read_tsv(output)[0]["host_action"], "remove")

    def test_validate_keeps_extra_columns_and_guards_host_group(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata.tsv"
            output = root / "validated.tsv"
            fields = [
                "sample_id",
                "fastq_1",
                "fastq_2",
                "layout",
                "age_group",
                "host_ref",
                "host_group",
                "source_sink",
                "model_source_env",
                "unit_id",
                "custom_panel_note",
            ]
            row = {
                "sample_id": "dog_1",
                "fastq_1": "reads.fq.gz",
                "fastq_2": "NA",
                "layout": "SE",
                "age_group": "ancient",
                "host_ref": "default",
                "host_group": "Canid",
                "source_sink": "source",
                "model_source_env": "canid_paleofeces",
                "unit_id": "dog_unit",
                "custom_panel_note": "keep me",
            }
            write_tsv(metadata, fields, [row])

            failed = run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                "--require-explicit-host-for-nondefault",
                check=False,
            )
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("requires an explicit host_ref", failed.stderr)

            row["host_ref"] = "dog.fa"
            write_tsv(metadata, fields, [row])
            run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                "--require-explicit-host-for-nondefault",
            )
            validated = read_tsv(output)[0]
            self.assertEqual(validated["custom_panel_note"], "keep me")
            self.assertEqual(validated["model_source_env"], "canid_paleofeces")
            self.assertEqual(validated["host_group"], "canid")
            self.assertEqual(validated["unit_id"], "dog_unit")
            self.assertEqual(validated["preprocess_mode"], "adapterremoval_se")

            na_row = dict(row)
            na_row.update(
                {
                    "sample_id": "human_na",
                    "age_group": "modern",
                    "host_ref": "NA",
                    "host_group": "NA",
                    "source_sink": "sink",
                    "model_source_env": "NA",
                    "unit_id": "NA",
                }
            )
            write_tsv(metadata, fields, [na_row])
            run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                "--require-explicit-host-for-nondefault",
            )
            normalized = read_tsv(output)[0]
            self.assertEqual(normalized["preprocess_mode"], "fastp_se")
            self.assertEqual(normalized["host_group"], "human")
            self.assertEqual(normalized["model_source_env"], "NA")
            self.assertEqual(normalized["unit_id"], "human_na")

            na_row["source_sink"] = "source"
            write_tsv(metadata, fields, [na_row])
            no_env = run_script("validate_metadata.py", "--input", metadata, "--output", output, check=False)
            self.assertNotEqual(no_env.returncode, 0)
            self.assertIn("source rows require model_source_env", no_env.stderr)

    def test_validate_rejects_unsafe_ids_and_empty_manifests(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fq = "reads.fq.gz"
            fields = ["sample_id", "fastq_1", "fastq_2", "layout", "age_group", "unit_id"]
            metadata = root / "metadata.tsv"
            output = root / "validated.tsv"
            write_tsv(
                metadata,
                fields,
                [
                    {
                        "sample_id": "unsafe/name",
                        "fastq_1": fq,
                        "fastq_2": "NA",
                        "layout": "SE",
                        "age_group": "modern",
                        "unit_id": "unit with spaces",
                    }
                ],
            )
            failed = run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                check=False,
            )
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("sample_id 'unsafe/name' is not filesystem-safe", failed.stderr)
            self.assertIn("unit_id 'unit with spaces' is not filesystem-safe", failed.stderr)

            write_tsv(metadata, fields, [])
            empty = run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                check=False,
            )
            self.assertNotEqual(empty.returncode, 0)
            self.assertIn("metadata contains no sample rows", empty.stderr)

            write_tsv(
                metadata,
                fields,
                [
                    {
                        "sample_id": "NA",
                        "fastq_1": fq,
                        "fastq_2": "NA",
                        "layout": "SE",
                        "age_group": "modern",
                        "unit_id": "NA",
                    }
                ],
            )
            missing_id = run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                check=False,
            )
            self.assertNotEqual(missing_id.returncode, 0)
            self.assertIn("reserved missing-value token", missing_id.stderr)

    def test_build_requires_one_to_one_and_reports_kraken_qc(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata.tsv"
            write_tsv(metadata, ["sample_id"], [{"sample_id": "s1"}, {"sample_id": "s2"}])
            bracken_header = [
                "name",
                "taxonomy_id",
                "taxonomy_lvl",
                "kraken_assigned_reads",
                "added_reads",
                "new_est_reads",
                "fraction_total_reads",
            ]
            bracken = []
            reports = []
            for sample, count in (("s1", 800), ("s2", 400)):
                table = root / f"{sample}.bracken.tsv"
                write_tsv(
                    table,
                    bracken_header,
                    [
                        {
                            "name": "Bacteria",
                            "taxonomy_id": "2",
                            "taxonomy_lvl": "D",
                            "kraken_assigned_reads": count,
                            "added_reads": 0,
                            "new_est_reads": count,
                            "fraction_total_reads": 1,
                        }
                    ],
                )
                bracken.append(table)
                report = root / f"{sample}.kraken2.report.txt"
                report.write_text(
                    "20.00\t200\t200\tU\t0\tunclassified\n"
                    "80.00\t800\t0\tR\t1\troot\n"
                )
                reports.append(report)

            run_script(
                "build_feature_table.py",
                "--metadata",
                metadata,
                "--prefix",
                "bracken_G",
                "--taxon-filter",
                "none",
                *bracken,
                "--kraken-reports",
                *reports,
                cwd=root,
            )
            qc = {row["sample_id"]: row for row in read_tsv(root / "bracken_G_none_sample_qc.tsv")}
            self.assertEqual(qc["s1"]["kraken_classified_fraction"], "0.8")
            self.assertEqual(qc["s1"]["bracken_total_reads"], "800")

            incomplete = run_script(
                "build_feature_table.py",
                "--metadata",
                metadata,
                "--prefix",
                "bad",
                bracken[0],
                cwd=root,
                check=False,
            )
            self.assertNotEqual(incomplete.returncode, 0)
            self.assertIn("missing Bracken table for: s2", incomplete.stderr)

    def test_build_preserves_header_only_bracken_sample_as_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata.tsv"
            write_tsv(
                metadata,
                ["sample_id"],
                [{"sample_id": "blank1"}, {"sample_id": "blank2"}, {"sample_id": "source1"}],
            )
            bracken_header = [
                "name",
                "taxonomy_id",
                "taxonomy_lvl",
                "kraken_assigned_reads",
                "added_reads",
                "new_est_reads",
                "fraction_total_reads",
            ]
            bracken_files = []
            for sample in ("blank1", "blank2"):
                table = root / f"{sample}.bracken.tsv"
                write_tsv(table, bracken_header, [])
                bracken_files.append(table)
            source_table = root / "source1.bracken.tsv"
            write_tsv(
                source_table,
                bracken_header,
                [
                    {
                        "name": "Bacteria",
                        "taxonomy_id": "2",
                        "taxonomy_lvl": "D",
                        "kraken_assigned_reads": 8,
                        "added_reads": 0,
                        "new_est_reads": 8,
                        "fraction_total_reads": 1,
                    }
                ],
            )
            bracken_files.append(source_table)

            blank1_report = root / "blank1.kraken2.report.txt"
            blank1_report.write_text(
                "71.11\t64\t64\tU\t0\tunclassified\n"
                "28.89\t26\t0\tR\t1\troot\n"
                "7.78\t7\t7\tG\t101\t  Genus_one\n"
                "4.44\t4\t4\tG\t102\t  Genus_two\n"
                "2.22\t2\t2\tG\t103\t  Genus_three\n"
                "2.22\t2\t2\tG\t104\t  Genus_four\n"
                "2.22\t2\t2\tG\t105\t  Genus_five\n"
            )
            blank2_report = root / "blank2.kraken2.report.txt"
            blank2_report.write_text("100.00\t12\t12\tU\t0\tunclassified\n")
            source_report = root / "source1.kraken2.report.txt"
            source_report.write_text(
                "20.00\t2\t2\tU\t0\tunclassified\n"
                "80.00\t8\t0\tR\t1\troot\n"
            )

            run_script(
                "build_feature_table.py",
                "--metadata",
                metadata,
                "--prefix",
                "bracken_G",
                "--taxon-filter",
                "none",
                *bracken_files,
                "--kraken-reports",
                blank1_report,
                blank2_report,
                source_report,
                cwd=root,
            )

            with open(root / "bracken_G_counts.tsv", newline="") as handle:
                self.assertEqual(
                    list(csv.reader(handle, delimiter="\t")),
                    [
                        ["#OTU ID", "blank1", "blank2", "source1"],
                        ["taxid_2__Bacteria", "0", "0", "8"],
                    ],
                )
            qc = {
                row["sample_id"]: row
                for row in read_tsv(root / "bracken_G_none_sample_qc.tsv")
            }
            self.assertEqual(qc["blank1"]["kraken_total_reads"], "90")
            self.assertEqual(qc["blank1"]["kraken_classified_reads"], "26")
            self.assertEqual(qc["blank1"]["kraken_unclassified_reads"], "64")
            self.assertEqual(qc["blank1"]["kraken_classified_fraction"], "0.288888888889")
            self.assertEqual(qc["blank1"]["bracken_total_reads"], "0")
            self.assertEqual(qc["blank1"]["bracken_fraction_of_kraken_classified"], "0")
            self.assertEqual(qc["blank1"]["taxon_filter_retained_reads"], "0")
            self.assertEqual(qc["blank1"]["taxon_filter_retained_fraction"], "NA")
            self.assertEqual(qc["blank1"]["raw_nonzero_features"], "0")
            self.assertEqual(qc["blank1"]["filtered_nonzero_features"], "0")
            self.assertEqual(qc["blank2"]["kraken_total_reads"], "12")
            self.assertEqual(qc["blank2"]["kraken_classified_reads"], "0")
            self.assertEqual(qc["blank2"]["kraken_classified_fraction"], "0")
            self.assertEqual(qc["blank2"]["bracken_fraction_of_kraken_classified"], "NA")
            self.assertEqual(qc["source1"]["bracken_total_reads"], "8")
            self.assertEqual(qc["source1"]["raw_nonzero_features"], "1")

    def _integer_count_fixture(self, root, new_est_reads):
        metadata = root / "metadata.tsv"
        write_tsv(metadata, ["sample_id"], [{"sample_id": "one"}, {"sample_id": "two"}])
        bracken_header = [
            "name",
            "taxonomy_id",
            "taxonomy_lvl",
            "kraken_assigned_reads",
            "added_reads",
            "new_est_reads",
            "fraction_total_reads",
        ]
        tables = []
        for index, sample in enumerate(("one", "two")):
            table = root / f"{sample}.bracken.tsv"
            write_tsv(
                table,
                bracken_header,
                [
                    {
                        "name": "Bacteria",
                        "taxonomy_id": "2",
                        "taxonomy_lvl": "D",
                        "kraken_assigned_reads": 5,
                        "added_reads": 0,
                        "new_est_reads": new_est_reads[index],
                        "fraction_total_reads": 1,
                    }
                ],
            )
            tables.append(table)
        return metadata, tables

    def test_bracken_counts_are_exact_integers_and_reject_fractions(self):
        """Counts must survive as exact integers, and fractions must fail loudly."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            # 8, 8.0 and 8e0 are all the integer 8; a huge count must not become
            # exponent notation, which %.12g would have produced.
            metadata, tables = self._integer_count_fixture(root, ["8.0", "9007199254740993"])
            run_script(
                "build_feature_table.py",
                "--metadata", metadata,
                "--prefix", "bracken_G",
                "--taxon-filter", "none",
                *tables,
                cwd=root,
            )
            with open(root / "bracken_G_counts.tsv", newline="") as handle:
                self.assertEqual(
                    list(csv.reader(handle, delimiter="\t")),
                    [
                        ["#OTU ID", "one", "two"],
                        ["taxid_2__Bacteria", "8", "9007199254740993"],
                    ],
                )

            metadata, tables = self._integer_count_fixture(root, ["8e0", "3"])
            run_script(
                "build_feature_table.py",
                "--metadata", metadata,
                "--prefix", "bracken_G",
                "--taxon-filter", "none",
                *tables,
                cwd=root,
            )
            counts = read_tsv(root / "bracken_G_counts.tsv")[0]
            self.assertEqual(counts["one"], "8")

            # Relative abundance stays floating point.
            relative = read_tsv(root / "bracken_G_relative_abundance.tsv")[0]
            self.assertEqual(relative["one"], "1")

            # "1e999999999" is finite to Decimal but would expand to a
            # billion-digit integer inside a non-interruptible C call. It must
            # be rejected, not materialised.
            for bad in ("", "8.5", "-1", "nan", "inf", "eight", "1e999999999"):
                metadata, tables = self._integer_count_fixture(root, [bad, "3"])
                failed = run_script(
                    "build_feature_table.py",
                    "--metadata", metadata,
                    "--prefix", "bad",
                    "--taxon-filter", "none",
                    *tables,
                    cwd=root,
                    check=False,
                )
                self.assertNotEqual(failed.returncode, 0, bad)
                # Failures identify the file, the line, the field and the sample.
                self.assertIn("one.bracken.tsv:2:new_est_reads", failed.stderr, bad)
                self.assertIn("sample one", failed.stderr, bad)

    def test_feature_table_is_invariant_to_bracken_file_order(self):
        """Terminal collect() order must not change any published table."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, tables = self._integer_count_fixture(root, ["8", "3"])
            outputs = [
                "bracken_G_counts.tsv",
                "bracken_G_relative_abundance.tsv",
                "bracken_G_features.tsv",
                "bracken_G_none_counts.tsv",
                "bracken_G_none_sample_qc.tsv",
            ]
            digests = []
            for order in (tables, list(reversed(tables))):
                run_script(
                    "build_feature_table.py",
                    "--metadata", metadata,
                    "--prefix", "bracken_G",
                    "--taxon-filter", "none",
                    *order,
                    cwd=root,
                )
                digests.append({name: (root / name).read_bytes() for name in outputs})
            self.assertEqual(digests[0], digests[1])
            self.assertEqual(read_tsv(root / "bracken_G_counts.tsv")[0]["one"], "8")

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
            ["blank1", "blank_unit", "", "blank", "P3", "skip", "NA", "bad", "bad", "bad", "bad", "bad"],
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
            "bracken_total_reads",
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
                    "bracken_total_reads": 2500,
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

    def test_prepare_aggregates_units_and_leaves_out_unselected_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            rows = read_tsv(metadata)
            for row in rows:
                if row["sample_id"].startswith("bone"):
                    row["source_sink"] = "NA"
            write_tsv(metadata, list(rows[0]), rows)
            self._prepare(root, metadata, counts, feature_qc, "st",
                          "--min-source-count", "1000", "--min-sink-count", "1000")
            mapping = {row["#SampleID"]: row for row in read_tsv(root / "st_map.tsv")}
            self.assertEqual(set(mapping), {"soil_unit", "sink_unit"})
            self.assertEqual(mapping["soil_unit"]["aggregation_n"], "2")
            self.assertEqual(mapping["soil_unit"]["Env"], "soil")
            self.assertEqual(mapping["soil_unit"]["project"], "P2")
            self.assertEqual(mapping["soil_unit"]["host_action"], "skip")
            qc = {row["input_sample_id"]: row for row in read_tsv(root / "st_qc.tsv")}
            self.assertEqual(qc["soil1"]["analysis_sample_id"], "soil_unit")
            self.assertEqual(qc["soil1"]["included"], "true")
            self.assertEqual(qc["bone1"]["exclusion_reason"], "source_sink_not_selected")
            self.assertEqual(qc["blank1"]["exclusion_reason"], "source_sink_not_selected")
            with open(root / "st_counts.tsv", newline="") as handle:
                header = next(csv.reader(handle, delimiter="\t"))
            self.assertEqual(header, ["#OTU ID", "soil_unit", "sink_unit"])

    def test_prepare_reports_undefined_qc_as_na_and_rejects_fractional_counts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)

            # A preserved zero-estimate Bracken sample has no defined retained
            # fraction. That must be reported as NA, not as an empty cell.
            qc_rows = read_tsv(feature_qc)
            for row in qc_rows:
                if row["sample_id"] == "sink1":
                    row["bracken_total_reads"] = "0"
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
            self.assertEqual(qc["sink1"]["bracken_total_reads"], "0")
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

    def test_all_taxon_scope_uses_bracken_yield_not_filtered_yield(self):
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
            self.assertEqual(mapping, {"bone_unit", "soil_unit", "sink_unit"})

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
