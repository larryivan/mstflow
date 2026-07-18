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
            fq = root / "reads.fq.gz"
            fq.touch()
            host = root / "human.fa"
            host.write_text(">chr1\nACGT\n")
            for suffix in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
                (root / f"human.fa{suffix}").write_text("index\n")

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
                "fastq_1": str(fq),
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
                "--default-host-fasta",
                host,
                "--require-explicit-host-action",
            )
            normalized = read_tsv(output)[0]
            self.assertEqual(normalized["host_action"], "skip")
            self.assertEqual(normalized["host_group"], "none")
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
            self.assertIn("host_action=skip requires host_ref=NA", conflicted.stderr)

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
                "--default-host-fasta",
                host,
                "--require-explicit-host-action",
            )
            self.assertEqual(read_tsv(output)[0]["host_action"], "remove")

            (root / "human.fa.sa").unlink()
            missing_index = run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                "--default-host-fasta",
                host,
                check=False,
            )
            self.assertNotEqual(missing_index.returncode, 0)
            self.assertIn("BWA index sidecar does not exist", missing_index.stderr)

    def test_validate_preserves_panel_fields_and_guards_host_group(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fq = root / "reads.fq.gz"
            host = root / "dog.fa"
            fq.touch()
            host.touch()
            metadata = root / "metadata.tsv"
            output = root / "validated.tsv"
            fields = [
                "sample_id",
                "fastq_1",
                "fastq_2",
                "layout",
                "age_group",
                "preprocess_mode",
                "host_ref",
                "host_group",
                "source_sink",
                "model_source_env",
                "category",
                "unit_id",
                "panel_role",
                "panel_layer",
                "validation_group",
                "project_group",
                "custom_panel_note",
            ]
            row = {
                "sample_id": "dog_1",
                "fastq_1": str(fq),
                "fastq_2": "NA",
                "layout": "SE",
                "age_group": "ancient",
                "host_ref": "default",
                "host_group": "Canid",
                "source_sink": "source",
                "category": "canid_paleofeces",
                "unit_id": "dog_unit",
                "panel_role": "source",
                "panel_layer": "reference",
                "validation_group": "VG1",
                "project_group": "PG1",
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

            row["host_ref"] = str(host)
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

            na_row = dict(row)
            na_row.update(
                {
                    "sample_id": "human_na",
                    "age_group": "modern",
                    "preprocess_mode": "NA",
                    "host_ref": "NA",
                    "host_group": "NA",
                    "model_source_env": "NA",
                    "unit_id": "NA",
                    "category": "soil",
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
            self.assertEqual(normalized["model_source_env"], "soil")
            self.assertEqual(normalized["unit_id"], "human_na")

    def test_validate_rejects_unsafe_ids_and_empty_manifests(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fq = root / "reads.fq.gz"
            fq.touch()
            fields = ["sample_id", "fastq_1", "fastq_2", "layout", "age_group", "unit_id"]
            metadata = root / "metadata.tsv"
            output = root / "validated.tsv"
            write_tsv(
                metadata,
                fields,
                [
                    {
                        "sample_id": "unsafe/name",
                        "fastq_1": str(fq),
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
                        "fastq_1": str(fq),
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

            bracken_fields = fields + ["bracken_read_len"]
            write_tsv(
                metadata,
                bracken_fields,
                [
                    {
                        "sample_id": "safe_id",
                        "fastq_1": str(fq),
                        "fastq_2": "NA",
                        "layout": "SE",
                        "age_group": "modern",
                        "unit_id": "safe_unit",
                        "bracken_read_len": '75"; touch injected; echo "',
                    }
                ],
            )
            unsafe_length = run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                check=False,
            )
            self.assertNotEqual(unsafe_length.returncode, 0)
            self.assertIn("bracken_read_len must be auto or a positive integer", unsafe_length.stderr)

    def test_build_requires_one_to_one_and_reports_kraken_qc(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata.tsv"
            fields = ["sample_id", "unit_id", "panel_role", "category"]
            write_tsv(
                metadata,
                fields,
                [
                    {"sample_id": "s1", "unit_id": "u1", "panel_role": "source", "category": "soil"},
                    {"sample_id": "s2", "unit_id": "u2", "panel_role": "sink", "category": "bone"},
                ],
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
            self.assertEqual(qc["s1"]["unit_id"], "u1")

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
            metadata_fields = ["sample_id", "unit_id", "panel_role", "category"]
            write_tsv(
                metadata,
                metadata_fields,
                [
                    {
                        "sample_id": "blank1",
                        "unit_id": "blank_unit",
                        "panel_role": "control",
                        "category": "extraction_library_blank",
                    },
                    {
                        "sample_id": "blank2",
                        "unit_id": "blank_unit_2",
                        "panel_role": "control",
                        "category": "extraction_library_blank",
                    },
                    {
                        "sample_id": "source1",
                        "unit_id": "source_unit",
                        "panel_role": "source",
                        "category": "soil",
                    },
                ],
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

    def _source_tracker_fixture(self, root):
        metadata = root / "metadata.tsv"
        fields = [
            "sample_id",
            "unit_id",
            "source_sink",
            "model_source_env",
            "source_class",
            "panel_role",
            "panel_layer",
            "category",
            "validation_group",
            "project_group",
            "host_action",
            "host_group",
            "SourceSink",
            "Env",
            "aggregation_n",
            "aggregation_members",
            "#SampleID",
        ]
        rows = [
            ["bone1", "bone_unit", "source", "bone", "bone", "source", "reference", "bone", "VG", "PG", "remove", "human", "bad", "bad", "bad", "bad", "bad"],
            ["bone2", "bone_unit", "source", "bone", "bone", "source", "reference", "bone", "VG", "PG", "remove", "human", "bad", "bad", "bad", "bad", "bad"],
            ["soil1", "soil_unit", "source", "soil", "soil", "source", "reference", "soil", "VG", "PG", "skip", "none", "bad", "bad", "bad", "bad", "bad"],
            ["soil2", "soil_unit", "source", "soil", "soil", "source", "reference", "soil", "VG", "PG", "skip", "none", "bad", "bad", "bad", "bad", "bad"],
            ["sink1", "sink_unit", "sink", "bone", "bone", "sink", "target", "bone", "VG", "PG", "remove", "human", "bad", "bad", "bad", "bad", "bad"],
            ["blank1", "blank_unit", "source", "blank", "blank", "control", "qc", "blank", "VG", "PG", "skip", "none", "bad", "bad", "bad", "bad", "bad"],
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

    def test_prepare_aggregates_units_and_reproduces_missing_source_experiment(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            run_script(
                "prepare_sourcetracker2_inputs.py",
                "--counts",
                counts,
                "--metadata",
                metadata,
                "--feature-sample-qc",
                feature_qc,
                "--output-counts",
                root / "st_counts.tsv",
                "--output-map",
                root / "st_map.tsv",
                "--summary",
                root / "summary.tsv",
                "--exclusions",
                root / "excluded.tsv",
                "--sample-qc",
                root / "st_qc.tsv",
                "--aggregate-by",
                "unit_id",
                "--exclude-source-env",
                "bone",
                "--exclude-panel-role",
                "control,holdout",
                "--source-rarefaction-depth",
                "1000",
                "--sink-rarefaction-depth",
                "1000",
                "--analysis-label",
                "NO_BONE",
            )
            mapping = {row["#SampleID"]: row for row in read_tsv(root / "st_map.tsv")}
            self.assertEqual(set(mapping), {"soil_unit", "sink_unit"})
            self.assertEqual(mapping["soil_unit"]["aggregation_n"], "2")
            self.assertEqual(mapping["soil_unit"]["validation_group"], "VG")
            self.assertEqual(mapping["soil_unit"]["host_action"], "skip")
            self.assertEqual(mapping["soil_unit"]["host_group"], "none")
            qc = {row["input_sample_id"]: row for row in read_tsv(root / "st_qc.tsv")}
            self.assertEqual(qc["soil1"]["host_action"], "skip")
            self.assertEqual(qc["soil1"]["host_group"], "none")
            excluded = read_tsv(root / "excluded.tsv")
            reasons = {row["input_sample_id"]: row["exclusion_reason"] for row in excluded}
            self.assertTrue(reasons["bone1"].startswith("source_environment_excluded"))
            self.assertTrue(reasons["blank1"].startswith("panel_role_excluded"))
            with open(root / "st_counts.tsv", newline="") as handle:
                header = next(csv.reader(handle, delimiter="\t"))
            self.assertEqual(header, ["#OTU ID", "soil_unit", "sink_unit"])

    def test_prepare_loo_does_not_require_sink(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            run_script(
                "prepare_sourcetracker2_inputs.py",
                "--counts",
                counts,
                "--metadata",
                metadata,
                "--feature-sample-qc",
                feature_qc,
                "--output-counts",
                root / "loo_counts.tsv",
                "--output-map",
                root / "loo_map.tsv",
                "--summary",
                root / "loo_summary.tsv",
                "--exclusions",
                root / "loo_excluded.tsv",
                "--sample-qc",
                root / "loo_qc.tsv",
                "--aggregate-by",
                "sample_id",
                "--exclude-panel-role",
                "control,holdout",
                "--sink-rarefaction-depth",
                "500",
                "--random-seed",
                "17",
                "--loo",
            )
            mapping = read_tsv(root / "loo_map.tsv")
            self.assertTrue(mapping)
            self.assertEqual({row["SourceSink"] for row in mapping}, {"source"})
            with open(root / "loo_counts.tsv", newline="") as handle:
                count_rows = list(csv.reader(handle, delimiter="\t"))
            column_totals = [
                sum(float(row[index]) for row in count_rows[1:])
                for index in range(1, len(count_rows[0]))
            ]
            self.assertEqual(column_totals, [500.0] * len(column_totals))
            excluded = {row["input_sample_id"]: row for row in read_tsv(root / "loo_excluded.tsv")}
            self.assertEqual(excluded["sink1"]["exclusion_reason"], "loo_sources_only")

    def test_prepare_loo_rejects_rarefied_sources_and_singleton_environments(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            common = [
                "--counts",
                counts,
                "--metadata",
                metadata,
                "--feature-sample-qc",
                feature_qc,
                "--output-counts",
                root / "loo_counts.tsv",
                "--output-map",
                root / "loo_map.tsv",
                "--summary",
                root / "loo_summary.tsv",
                "--exclude-panel-role",
                "control,holdout",
                "--loo",
            ]
            rarefied = run_script(
                "prepare_sourcetracker2_inputs.py",
                *common,
                "--source-rarefaction-depth",
                "1000",
                check=False,
            )
            self.assertNotEqual(rarefied.returncode, 0)
            self.assertIn("LOO requires --source-rarefaction-depth 0", rarefied.stderr)

            singleton = run_script(
                "prepare_sourcetracker2_inputs.py",
                *common,
                "--aggregate-by",
                "unit_id",
                check=False,
            )
            self.assertNotEqual(singleton.returncode, 0)
            self.assertIn("at least two independent analysis units", singleton.stderr)
            self.assertIn("bone", singleton.stderr)
            self.assertIn("soil", singleton.stderr)

    def test_prepare_validates_and_filters_complete_units(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            rows = read_tsv(metadata)
            for row in rows:
                if row["sample_id"] == "bone1":
                    row["panel_role"] = "holdout"
            write_tsv(metadata, list(rows[0]), rows)
            conflicting = run_script(
                "prepare_sourcetracker2_inputs.py",
                "--counts",
                counts,
                "--metadata",
                metadata,
                "--feature-sample-qc",
                feature_qc,
                "--output-counts",
                root / "conflict_counts.tsv",
                "--output-map",
                root / "conflict_map.tsv",
                "--summary",
                root / "conflict_summary.tsv",
                "--aggregate-by",
                "unit_id",
                "--exclude-panel-role",
                "control,holdout",
                check=False,
            )
            self.assertNotEqual(conflicting.returncode, 0)
            self.assertIn("bone_unit", conflicting.stderr)
            self.assertIn("panel_role", conflicting.stderr)

            rows = read_tsv(metadata)
            for row in rows:
                if row["sample_id"] == "bone1":
                    row["panel_role"] = "source"
            write_tsv(metadata, list(rows[0]), rows)
            qc_rows = read_tsv(feature_qc)
            for row in qc_rows:
                if row["sample_id"] == "bone1":
                    row["kraken_classified_reads"] = "0"
                    row["kraken_classified_fraction"] = "0"
            write_tsv(feature_qc, list(qc_rows[0]), qc_rows)
            run_script(
                "prepare_sourcetracker2_inputs.py",
                "--counts",
                counts,
                "--metadata",
                metadata,
                "--feature-sample-qc",
                feature_qc,
                "--output-counts",
                root / "filtered_counts.tsv",
                "--output-map",
                root / "filtered_map.tsv",
                "--summary",
                root / "filtered_summary.tsv",
                "--exclusions",
                root / "filtered_exclusions.tsv",
                "--aggregate-by",
                "unit_id",
                "--exclude-panel-role",
                "control,holdout",
                "--min-kraken-classified-fraction",
                "0.15",
            )
            mapping = {row["#SampleID"] for row in read_tsv(root / "filtered_map.tsv")}
            self.assertNotIn("bone_unit", mapping)
            excluded = {
                row["input_sample_id"]: row["exclusion_reason"]
                for row in read_tsv(root / "filtered_exclusions.tsv")
            }
            self.assertTrue(excluded["bone1"].startswith("below_min_kraken"))
            self.assertEqual(excluded["bone1"], excluded["bone2"])

    def test_prepare_allows_explicit_holdout_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            rows = read_tsv(metadata)
            for row in rows:
                if row["sample_id"] == "blank1":
                    row["panel_role"] = "holdout"
            write_tsv(metadata, list(rows[0]), rows)
            run_script(
                "prepare_sourcetracker2_inputs.py",
                "--counts",
                counts,
                "--metadata",
                metadata,
                "--feature-sample-qc",
                feature_qc,
                "--output-counts",
                root / "holdout_counts.tsv",
                "--output-map",
                root / "holdout_map.tsv",
                "--summary",
                root / "holdout_summary.tsv",
                "--aggregate-by",
                "unit_id",
            )
            mapping = {row["#SampleID"] for row in read_tsv(root / "holdout_map.tsv")}
            self.assertIn("blank_unit", mapping)

    def test_prepare_rejects_reserved_source_environment_names(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata, counts, feature_qc = self._source_tracker_fixture(root)
            rows = read_tsv(metadata)
            for row in rows:
                if row["sample_id"].startswith("soil"):
                    row["model_source_env"] = "Unknown"
            write_tsv(metadata, list(rows[0]), rows)
            failed = run_script(
                "prepare_sourcetracker2_inputs.py",
                "--counts",
                counts,
                "--metadata",
                metadata,
                "--feature-sample-qc",
                feature_qc,
                "--output-counts",
                root / "reserved_counts.tsv",
                "--output-map",
                root / "reserved_map.tsv",
                "--summary",
                root / "reserved_summary.tsv",
                "--aggregate-by",
                "unit_id",
                "--exclude-panel-role",
                "control,holdout",
                check=False,
            )
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
            run_script(
                "prepare_sourcetracker2_inputs.py",
                "--counts",
                counts,
                "--metadata",
                metadata,
                "--feature-sample-qc",
                feature_qc,
                "--output-counts",
                root / "all_counts.tsv",
                "--output-map",
                root / "all_map.tsv",
                "--summary",
                root / "all_summary.tsv",
                "--exclusions",
                root / "all_excluded.tsv",
                "--sample-qc",
                root / "all_qc.tsv",
                "--taxon-scope",
                "all",
                "--aggregate-by",
                "unit_id",
                "--exclude-source-env",
                "bone",
                "--exclude-panel-role",
                "control,holdout",
            )
            mapping = {row["#SampleID"] for row in read_tsv(root / "all_map.tsv")}
            self.assertEqual(mapping, {"soil_unit", "sink_unit"})
            qc = read_tsv(root / "all_qc.tsv")
            self.assertEqual({row["taxon_scope"] for row in qc}, {"all"})

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

    def test_plot_uses_uncertainty_and_abstains_without_hard_call(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            means = root / "mixing_proportions.txt"
            stds = root / "mixing_proportions_stds.txt"
            mapping = root / "map.tsv"
            means.write_text("SampleID\tsoil\tbone\tUnknown\nsink1\t0.2\t0.1\t0.7\n")
            stds.write_text("SampleID\tsoil\tbone\tUnknown\nsink1\t0.05\t0.02\t0.08\n")
            write_tsv(
                mapping,
                ["#SampleID", "SourceSink", "Env"],
                [{"#SampleID": "sink1", "SourceSink": "sink", "Env": "bone"}],
            )
            assignments = root / "assignments.tsv"
            image = root / "plot.png"
            run_script(
                "plot_sourcetracker2.py",
                "--proportions",
                means,
                "--uncertainty",
                stds,
                "--mapping",
                mapping,
                "--assignments",
                assignments,
                "--plot",
                image,
                "--max-unknown-prop",
                "0.5",
                "--min-top-known-prop",
                "0.5",
            )
            result = read_tsv(assignments)[0]
            self.assertEqual(result["assignment"], "ABSTAIN")
            self.assertEqual(result["abstained"], "true")
            self.assertEqual(result["top_known_std"], "0.05")
            self.assertEqual(result["top_known_correct"], "NA")
            self.assertGreater(image.stat().st_size, 100)

            tie_mapping = root / "tie_map.tsv"
            write_tsv(
                tie_mapping,
                ["#SampleID", "SourceSink", "Env"],
                [{"#SampleID": "tie_sink", "SourceSink": "sink", "Env": "soil"}],
            )
            tie_means = root / "tie_means.txt"
            tie_means.write_text(
                "SampleID\tsoil\tbone\tUnknown\ntie_sink\t0.5\t0.5\t0\n"
            )
            tie_assignments = root / "tie_assignments.tsv"
            run_script(
                "plot_sourcetracker2.py",
                "--proportions",
                tie_means,
                "--mapping",
                tie_mapping,
                "--assignments",
                tie_assignments,
                "--plot",
                root / "tie.png",
            )
            tie_result = read_tsv(tie_assignments)[0]
            self.assertEqual(tie_result["assignment"], "ABSTAIN")
            self.assertIn("top_known_tie", tie_result["abstention_reason"])

            duplicate = root / "duplicate_columns.txt"
            duplicate.write_text(
                "SampleID\tsoil\tsoil\tUnknown\nsink1\t0.4\t0.4\t0.2\n"
            )
            duplicate_result = run_script(
                "plot_sourcetracker2.py",
                "--proportions",
                duplicate,
                "--mapping",
                mapping,
                "--assignments",
                root / "duplicate_assignments.tsv",
                "--plot",
                root / "duplicate.png",
                check=False,
            )
            self.assertNotEqual(duplicate_result.returncode, 0)
            self.assertIn("duplicate result column(s): soil", duplicate_result.stderr)

            empty = root / "empty.txt"
            empty.write_text("SampleID\tsoil\tUnknown\n")
            empty_image = root / "empty.png"
            missing_sink = run_script(
                "plot_sourcetracker2.py",
                "--proportions",
                empty,
                "--mapping",
                mapping,
                "--assignments",
                root / "empty_assignments.tsv",
                "--plot",
                empty_image,
                check=False,
            )
            self.assertNotEqual(missing_sink.returncode, 0)
            self.assertIn("missing sink result(s): sink1", missing_sink.stderr)

            malformed = root / "malformed.txt"
            malformed.write_text("SampleID\tsoil\tUnknown\nsink1\t0.2\t0.2\n")
            bad_sum = run_script(
                "plot_sourcetracker2.py",
                "--proportions",
                malformed,
                "--mapping",
                mapping,
                "--assignments",
                root / "malformed_assignments.tsv",
                "--plot",
                root / "malformed.png",
                check=False,
            )
            self.assertNotEqual(bad_sum.returncode, 0)
            self.assertIn("sum to 0.4", bad_sum.stderr)

            loo_mapping = root / "loo_map.tsv"
            write_tsv(
                loo_mapping,
                ["#SampleID", "SourceSink", "Env"],
                [{"#SampleID": "source1", "SourceSink": "source", "Env": "soil"}],
            )
            loo_means = root / "loo_means.txt"
            loo_means.write_text("SampleID\tsoil\tUnknown\nsource1\t1\t0\n")
            loo_image = root / "loo.png"
            run_script(
                "plot_sourcetracker2.py",
                "--proportions",
                loo_means,
                "--mapping",
                loo_mapping,
                "--assignments",
                root / "loo_assignments.tsv",
                "--plot",
                loo_image,
                "--loo",
            )
            self.assertGreater(loo_image.stat().st_size, 100)


if __name__ == "__main__":
    unittest.main()
