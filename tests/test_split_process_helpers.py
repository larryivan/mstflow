import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
FRAGMENT_HELPER = BIN / "mstflow_fragment_tools.py"


def write_fastq(path, records):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def read_fastq(path):
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        lines = handle.read().splitlines()
    return [tuple(lines[index:index + 4]) for index in range(0, len(lines), 4)]


def read_metrics(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def run_script(script, *arguments, cwd):
    completed = subprocess.run(
        [sys.executable, str(BIN / script), *map(str, arguments)],
        cwd=cwd,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(
            f"{script} failed ({completed.returncode})\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    return completed


class SplitProcessHelperTests(unittest.TestCase):
    def test_metadata_container_mode_leaves_paths_for_nextflow_staging(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metadata = root / "metadata.tsv"
            output = root / "validated.tsv"
            metadata.write_text(
                "sample_id\tfastq_1\tfastq_2\tlayout\tage_group\thost_action\t"
                "host_group\thost_ref\n"
                "sample\t/outside/not-mounted.fq.gz\tNA\tSE\tmodern\tremove\t"
                "dog\t/outside/not-mounted-dog.fa\n",
                encoding="utf-8",
            )
            run_script(
                "validate_metadata.py",
                "--input",
                metadata,
                "--output",
                output,
                "--default-host-group",
                "human",
                "--require-explicit-host-for-nondefault",
                "--skip-filesystem-checks",
                cwd=root,
            )
            with output.open(encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["fastq_1"], "/outside/not-mounted.fq.gz")
            self.assertEqual(row["host_ref"], "/outside/not-mounted-dog.fa")

    def test_adapterremoval_collapsed_policy_is_normalized_separately(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw"
            raw.mkdir()
            sample = "ancient"
            write_fastq(raw / f"{sample}.ar.collapsed.gz", [("collapsed", "A" * 40)])
            write_fastq(raw / f"{sample}.ar.singleton.truncated.gz", [("singleton", "C" * 40)])
            write_fastq(raw / f"{sample}.ar.pair1.truncated.gz", [("pair/1", "G" * 40)])
            write_fastq(raw / f"{sample}.ar.pair2.truncated.gz", [("pair/2", "T" * 40)])
            run_script(
                "normalize_adapterremoval.py",
                "--sample-id",
                sample,
                "--raw-dir",
                raw,
                "--mode",
                "collapsed_only",
                "--kmer-length",
                "35",
                "--gzip-level",
                "1",
                "--fragment-helper",
                FRAGMENT_HELPER,
                cwd=root,
            )
            normalized = root / f"{sample}.prehost.fq.gz"
            self.assertEqual(
                [record[0] for record in read_fastq(normalized)],
                ["@collapsed"],
            )
            first_bytes = normalized.read_bytes()
            run_script(
                "normalize_adapterremoval.py",
                "--sample-id",
                sample,
                "--raw-dir",
                raw,
                "--mode",
                "collapsed_only",
                "--kmer-length",
                "35",
                "--gzip-level",
                "1",
                "--fragment-helper",
                FRAGMENT_HELPER,
                cwd=root,
            )
            self.assertEqual(first_bytes, normalized.read_bytes())
            metrics = read_metrics(root / f"{sample}.adapterremoval.fragment_qc.tsv")
            self.assertEqual(metrics["fastq_gzip_level"], "1")
            self.assertEqual(metrics["retained_collapsed_fragments"], "1")
            self.assertEqual(metrics["discarded_singleton_fragments"], "1")
            self.assertEqual(metrics["discarded_unmerged_fragments"], "1")
            self.assertEqual(metrics["retained_fragments"], "1")

    def test_adapterremoval_all_fragments_preserves_three_file_contract(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = root / "raw"
            raw.mkdir()
            sample = "mixed"
            write_fastq(raw / f"{sample}.ar.collapsed.gz", [("collapsed", "A" * 40)])
            write_fastq(raw / f"{sample}.ar.pair1.truncated.gz", [("pair/1", "C" * 40)])
            write_fastq(raw / f"{sample}.ar.pair2.truncated.gz", [("pair/2", "G" * 40)])
            run_script(
                "normalize_adapterremoval.py",
                "--sample-id",
                sample,
                "--raw-dir",
                raw,
                "--mode",
                "all_fragments",
                "--kmer-length",
                "35",
                "--gzip-level",
                "6",
                "--fragment-helper",
                FRAGMENT_HELPER,
                cwd=root,
            )
            outputs = [
                root / f"{sample}.prehost_se.fq.gz",
                root / f"{sample}.prehost_R1.fq.gz",
                root / f"{sample}.prehost_R2.fq.gz",
            ]
            self.assertTrue(all(path.is_file() for path in outputs))
            self.assertEqual([len(read_fastq(path)) for path in outputs], [1, 1, 1])
            metrics = read_metrics(root / f"{sample}.adapterremoval.fragment_qc.tsv")
            self.assertEqual(metrics["retained_single_fragments"], "1")
            self.assertEqual(metrics["retained_paired_fragments"], "1")
            self.assertEqual(metrics["retained_fragments"], "2")

    def test_compressed_bwa_sam_bundle_retains_or_removes_whole_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            alignments = root / "alignments"
            alignments.mkdir()
            (alignments / "mode.txt").write_text("modern_pe\n", encoding="utf-8")
            sequence1 = "A" * 40
            sequence2 = "C" * 40
            quality = "I" * 40
            with gzip.open(alignments / "pair.sam.gz", "wt", encoding="utf-8", newline="") as sam:
                sam.write("@HD\tVN:1.6\tSO:queryname\n")
                sam.write(f"host\t65\thost\t1\t0\t40M\t=\t1\t0\t{sequence1}\t{quality}\n")
                sam.write(f"host\t129\thost\t1\t0\t40M\t=\t1\t0\t{sequence2}\t{quality}\n")
                sam.write(f"keep\t77\t*\t0\t0\t*\t*\t0\t0\t{sequence1}\t{quality}\n")
                sam.write(f"keep\t141\t*\t0\t0\t*\t*\t0\t0\t{sequence2}\t{quality}\n")
            run_script(
                "filter_host_alignments.py",
                "--sample-id",
                "sample",
                "--alignment-dir",
                alignments,
                "--minimum-mapq",
                "0",
                "--gzip-level",
                "1",
                "--fragment-helper",
                FRAGMENT_HELPER,
                cwd=root,
            )
            self.assertEqual(
                [record[0] for record in read_fastq(root / "sample.dehosted_R1.fq.gz")],
                ["@keep/1"],
            )
            self.assertEqual(
                [record[0] for record in read_fastq(root / "sample.dehosted_R2.fq.gz")],
                ["@keep/2"],
            )
            metrics = read_metrics(root / "sample.host_removal.fragment_qc.tsv")
            self.assertEqual(metrics["input_fragments"], "2")
            self.assertEqual(metrics["retained_fragments"], "1")
            self.assertEqual(metrics["removed_host_fragments"], "1")

    def test_kraken_preparation_filters_complete_modern_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            r1 = root / "R1.fq.gz"
            r2 = root / "R2.fq.gz"
            write_fastq(r1, [("keep/1", "A" * 40), ("short/1", "A" * 20)])
            write_fastq(r2, [("keep/2", "C" * 40), ("short/2", "C" * 20)])
            run_script(
                "prepare_kraken_input.py",
                "--sample-id",
                "modern",
                "--age-group",
                "modern",
                "--read",
                r2,
                "--read",
                r1,
                "--minimum-usable-span",
                "35",
                "--kmer-length",
                "35",
                "--gzip-level",
                "6",
                "--fragment-helper",
                FRAGMENT_HELPER,
                cwd=root,
            )
            self.assertEqual((root / "modern.kraken_input.mode.txt").read_text(), "paired\n")
            self.assertEqual(
                [record[0] for record in read_fastq(root / "modern.kraken_input_R1.fq.gz")],
                ["@keep/1"],
            )
            self.assertEqual(
                [record[0] for record in read_fastq(root / "modern.kraken_input_R2.fq.gz")],
                ["@keep/2"],
            )
            metrics = read_metrics(root / "modern.kraken_input.read_length.tsv")
            self.assertEqual(metrics["retained_fragments"], "1")
            self.assertEqual(metrics["dropped_below_minimum_usable_span"], "1")

    def test_ancient_mixed_preparation_keeps_native_pair_semantics(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            collapsed = root / "sample.prehost_se.fq.gz"
            r1 = root / "sample.prehost_R1.fq.gz"
            r2 = root / "sample.prehost_R2.fq.gz"
            write_fastq(collapsed, [("collapsed", "A" * 60)])
            write_fastq(r1, [("pair/1", "C" * 40)])
            write_fastq(r2, [("pair/2", "G" * 40)])
            run_script(
                "prepare_kraken_input.py",
                "--sample-id",
                "ancient",
                "--age-group",
                "ancient",
                "--read",
                collapsed,
                "--read",
                r1,
                "--read",
                r2,
                "--minimum-usable-span",
                "35",
                "--kmer-length",
                "35",
                "--gzip-level",
                "6",
                "--fragment-helper",
                FRAGMENT_HELPER,
                cwd=root,
            )
            self.assertEqual((root / "ancient.kraken_input.mode.txt").read_text(), "paired\n")
            emitted_r1 = read_fastq(root / "ancient.kraken_input_R1.fq.gz")
            emitted_r2 = read_fastq(root / "ancient.kraken_input_R2.fq.gz")
            self.assertEqual(len(emitted_r1), 2)
            self.assertEqual(len(emitted_r2), 2)
            self.assertEqual(emitted_r1[0][1], "A" * 60)
            self.assertEqual(emitted_r2[0][1], "N")
            self.assertIn("mstflow_internal_dummy_mate=1", emitted_r2[0][0])
            self.assertEqual(emitted_r1[1][1], "C" * 40)
            self.assertEqual(emitted_r2[1][1], "G" * 40)
            metrics = read_metrics(root / "ancient.kraken_input.read_length.tsv")
            self.assertEqual(metrics["retained_single_fragments"], "1")
            self.assertEqual(metrics["retained_paired_fragments"], "1")
            self.assertEqual(metrics["retained_fragments"], "2")
            self.assertEqual(metrics["median_kraken_confidence_denominator"], "19")


if __name__ == "__main__":
    unittest.main()
