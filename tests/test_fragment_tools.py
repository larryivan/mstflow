import csv
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HELPER = ROOT / "bin" / "mstflow_fragment_tools.py"


def write_fastq(path, records):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def read_metrics(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return {
            row["metric"]: row["value"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def run_helper(*args, check=True):
    return subprocess.run(
        [sys.executable, str(HELPER), *map(str, args)],
        check=check,
        text=True,
        capture_output=True,
    )


class FragmentToolTests(unittest.TestCase):
    def test_paired_bracken_length_is_one_mate_not_combined_span(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            r1 = tmp / "R1.fq.gz"
            r2 = tmp / "R2.fq.gz"
            stats = tmp / "stats.tsv"
            write_fastq(r1, [("pair/1", "A" * 75)])
            write_fastq(r2, [("pair/2", "C" * 75)])

            run_helper(
                "stats",
                "--r1",
                r1,
                "--r2",
                r2,
                "--output",
                stats,
                "--minimum-usable-span",
                "35",
                "--kmer-length",
                "35",
            )
            metrics = read_metrics(stats)
            self.assertEqual(metrics["median_bracken_read_length"], "75")
            self.assertEqual(metrics["median_kmer_equivalent_read_length"], "116")

            database = tmp / "db"
            database.mkdir()
            (database / "database75mers.kmer_distrib").touch()
            (database / "database100mers.kmer_distrib").touch()
            selection = tmp / "selection.tsv"
            completed = run_helper(
                "select-bracken-length",
                "--stats",
                stats,
                "--database",
                database,
                "--requested",
                "auto",
                "--sample-id",
                "sample",
                "--output",
                selection,
            )
            self.assertEqual(completed.stdout.strip(), "75")
            with open(selection, encoding="utf-8", newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["observed_median_bracken_read_length"], "75")
            self.assertEqual(row["selected_bracken_read_length"], "75")

            tie_stats = tmp / "tie.tsv"
            tie_stats.write_text(
                "metric\tvalue\nmedian_bracken_read_length\t87.5\n",
                encoding="utf-8",
            )
            tie_selection = tmp / "tie_selection.tsv"
            tied = run_helper(
                "select-bracken-length",
                "--stats",
                tie_stats,
                "--database",
                database,
                "--requested",
                "auto",
                "--sample-id",
                "tie",
                "--output",
                tie_selection,
            )
            self.assertEqual(tied.stdout.strip(), "75")

    def test_combine_emits_one_record_per_pair_with_kmer_blocking_gap(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            r1 = tmp / "R1.fq.gz"
            r2 = tmp / "R2.fq.gz"
            output = tmp / "combined.fq.gz"
            stats = tmp / "stats.tsv"
            write_fastq(r1, [("pair/1", "A" * 40)])
            write_fastq(r2, [("pair/2", "C" * 40)])
            run_helper(
                "combine",
                "--r1",
                r1,
                "--r2",
                r2,
                "--output",
                output,
                "--stats",
                stats,
                "--gap-length",
                "100",
                "--kmer-length",
                "35",
                "--minimum-usable-span",
                "35",
            )
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                lines = handle.read().splitlines()
            self.assertEqual(len(lines), 4)
            self.assertEqual(lines[1], f"{'A' * 40}{'N' * 100}{'C' * 40}")
            self.assertEqual(read_metrics(stats)["retained_paired_fragments"], "1")

            rescanned_stats = tmp / "rescanned.tsv"
            run_helper(
                "stats",
                "--single",
                output,
                "--output",
                rescanned_stats,
                "--kmer-length",
                "35",
                "--minimum-usable-span",
                "35",
            )
            rescanned = read_metrics(rescanned_stats)
            self.assertEqual(rescanned["median_bracken_read_length"], "40")
            self.assertEqual(rescanned["median_kmer_equivalent_read_length"], "46")

    def test_pair_name_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            r1 = tmp / "R1.fq.gz"
            r2 = tmp / "R2.fq.gz"
            write_fastq(r1, [("first/1", "A" * 40)])
            write_fastq(r2, [("different/2", "C" * 40)])
            completed = run_helper(
                "stats",
                "--r1",
                r1,
                "--r2",
                r2,
                "--output",
                tmp / "stats.tsv",
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("paired FASTQ name mismatch", completed.stderr)

    def test_host_mapq_threshold_removes_or_retains_complete_pair(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            sequence1 = "A" * 40
            sequence2 = "C" * 40
            quality = "I" * 40
            sam = tmp / "pair.sam"
            sam.write_text(
                "@HD\tVN:1.6\tSO:queryname\n"
                f"pair\t65\thost\t1\t0\t40M\t=\t1\t0\t{sequence1}\t{quality}\n"
                f"pair\t129\thost\t1\t0\t40M\t=\t1\t0\t{sequence2}\t{quality}\n",
                encoding="utf-8",
            )

            removed_r1 = tmp / "removed_R1.fq.gz"
            removed_r2 = tmp / "removed_R2.fq.gz"
            removed_stats = tmp / "removed.tsv"
            run_helper(
                "filter-host-sam",
                "--sam",
                sam,
                "--output-r1",
                removed_r1,
                "--output-r2",
                removed_r2,
                "--minimum-mapq",
                "0",
                "--stats",
                removed_stats,
            )
            self.assertEqual(read_metrics(removed_stats)["removed_host_fragments"], "1")
            with gzip.open(removed_r1, "rt", encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "")

            retained_r1 = tmp / "retained_R1.fq.gz"
            retained_r2 = tmp / "retained_R2.fq.gz"
            retained_stats = tmp / "retained.tsv"
            run_helper(
                "filter-host-sam",
                "--sam",
                sam,
                "--output-r1",
                retained_r1,
                "--output-r2",
                retained_r2,
                "--minimum-mapq",
                "20",
                "--stats",
                retained_stats,
            )
            self.assertEqual(read_metrics(retained_stats)["retained_fragments"], "1")
            with gzip.open(retained_r1, "rt", encoding="utf-8") as handle:
                self.assertTrue(handle.read().startswith("@pair/1\n"))


if __name__ == "__main__":
    unittest.main()
