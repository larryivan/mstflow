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


def read_fastq(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    return [tuple(lines[index:index + 4]) for index in range(0, len(lines), 4)]


class FragmentToolTests(unittest.TestCase):
    def test_pair_fragments_is_kraken_equivalent_to_native_submission(self):
        """A sub-k dummy mate must not perturb any reported statistic.

        Kraken2 2.1.6 scans max(0, length - k + 1) minimizer positions from a
        sequence (src/classify.cc ClassifySequence, src/mmscanner.cc), so a mate
        shorter than k contributes exactly zero positions and cannot change
        total_kmers, required_score, or the call.
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            single = tmp / "collapsed.fq.gz"
            r1 = tmp / "R1.fq.gz"
            r2 = tmp / "R2.fq.gz"
            write_fastq(single, [("collapsed1", "A" * 60), ("collapsed2", "G" * 40)])
            write_fastq(r1, [("pair/1", "A" * 40)])
            write_fastq(r2, [("pair/2", "C" * 40)])

            out_r1 = tmp / "paired_R1.fq.gz"
            out_r2 = tmp / "paired_R2.fq.gz"
            stats = tmp / "paired.stats.tsv"
            run_helper(
                "pair-fragments",
                "--single", single,
                "--r1", r1,
                "--r2", r2,
                "--output-r1", out_r1,
                "--output-r2", out_r2,
                "--stats", stats,
                "--minimum-usable-span", "35",
                "--kmer-length", "35",
            )

            emitted_r1 = read_fastq(out_r1)
            emitted_r2 = read_fastq(out_r2)
            # One fragment in, one Kraken2 classification record out.
            self.assertEqual(len(emitted_r1), 3)
            self.assertEqual(len(emitted_r2), 3)
            self.assertEqual([record[0] for record in emitted_r1[:2]], ["@collapsed1", "@collapsed2"])
            for record in emitted_r2[:2]:
                self.assertIn("mstflow_internal_dummy_mate=1", record[0])
                self.assertEqual(record[1], "N")
                self.assertEqual(record[3], "!")
            # A real unmerged mate is passed through untouched.
            self.assertEqual(emitted_r2[2][1], "C" * 40)

            metrics = read_metrics(stats)
            self.assertEqual(metrics["retained_single_fragments"], "2")
            self.assertEqual(metrics["retained_paired_fragments"], "1")
            # Bracken observations: 60, 40 (collapsed) and 40, 40 (real mates);
            # the placeholder is excluded, so the median stays 40.
            self.assertEqual(metrics["median_bracken_read_length"], "40")
            self.assertEqual(metrics["mean_bracken_read_length"], "45.0")
            # Confidence denominators: 26, 6 (native SE) and 6+6 (native paired).
            self.assertEqual(metrics["median_kraken_confidence_denominator"], "12")

            # KRAKEN2_PROFILE re-scans the emitted pair; the statistics it
            # derives must match those recorded when the pair was written.
            rescanned = tmp / "rescanned.tsv"
            run_helper(
                "stats",
                "--r1", out_r1,
                "--r2", out_r2,
                "--output", rescanned,
                "--minimum-usable-span", "35",
                "--kmer-length", "35",
            )
            rescan = read_metrics(rescanned)
            for metric in (
                "median_bracken_read_length",
                "mean_bracken_read_length",
                "median_kraken_confidence_denominator",
                "mean_kraken_confidence_denominator",
                "median_informative_bases",
                "median_kmer_equivalent_read_length",
            ):
                self.assertEqual(rescan[metric], metrics[metric], metric)

    def test_n_spacer_concatenation_inflates_the_confidence_denominator(self):
        """Regression guard for the corrected all_fragments representation.

        Two 40 bp mates give Kraken2 12 minimizer positions natively but 146
        when concatenated with a 100 N spacer, because ambiguous positions are
        counted in total_kmers. required_score = ceil(confidence * total_kmers).
        """
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            r1 = tmp / "R1.fq.gz"
            r2 = tmp / "R2.fq.gz"
            write_fastq(r1, [("pair/1", "A" * 40)])
            write_fastq(r2, [("pair/2", "C" * 40)])

            spacer_stats = tmp / "spacer.tsv"
            run_helper(
                "combine",
                "--r1", r1,
                "--r2", r2,
                "--output", tmp / "combined.fq.gz",
                "--stats", spacer_stats,
                "--gap-length", "100",
                "--kmer-length", "35",
                "--minimum-usable-span", "35",
            )
            paired_stats = tmp / "paired.tsv"
            run_helper(
                "pair-fragments",
                "--r1", r1,
                "--r2", r2,
                "--output-r1", tmp / "p_R1.fq.gz",
                "--output-r2", tmp / "p_R2.fq.gz",
                "--stats", paired_stats,
                "--kmer-length", "35",
                "--minimum-usable-span", "35",
            )
            self.assertEqual(
                read_metrics(spacer_stats)["median_kraken_confidence_denominator"], "146"
            )
            self.assertEqual(
                read_metrics(paired_stats)["median_kraken_confidence_denominator"], "12"
            )
            # Fragment membership is identical between the two representations.
            self.assertEqual(
                read_metrics(spacer_stats)["retained_fragments"],
                read_metrics(paired_stats)["retained_fragments"],
            )

    def test_dummy_mate_must_be_shorter_than_the_kmer_length(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            single = tmp / "collapsed.fq.gz"
            write_fastq(single, [("collapsed1", "A" * 60)])
            completed = run_helper(
                "pair-fragments",
                "--single", single,
                "--output-r1", tmp / "R1.fq.gz",
                "--output-r2", tmp / "R2.fq.gz",
                "--stats", tmp / "stats.tsv",
                "--kmer-length", "35",
                "--dummy-mate-length", "35",
                check=False,
            )
            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("must be shorter than the k-mer length", completed.stderr)

    def test_gzip_output_is_reproducible_and_level_independent(self):
        """Compression level changes bytes on disk only, never the FASTQ."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp / "in.fq.gz"
            write_fastq(source, [(f"read{index}", "ACGT" * 20) for index in range(200)])

            def prepare(destination, level):
                run_helper(
                    "prepare-single",
                    "--single", source,
                    "--output", destination,
                    "--stats", tmp / f"stats{level}.tsv",
                    "--minimum-usable-span", "35",
                    "--kmer-length", "35",
                    "--gzip-level", str(level),
                )

            # Repeated writes of the same content are byte-identical: the GZIP
            # header mtime is pinned so checksum manifests stay stable.
            first = tmp / "a.fq.gz"
            prepare(first, 6)
            first_bytes = first.read_bytes()
            first.unlink()
            prepare(first, 6)
            self.assertEqual(first_bytes, first.read_bytes())

            reference = None
            for level in (0, 1, 6, 9):
                target = tmp / f"level{level}" / "out.fq.gz"
                target.parent.mkdir()
                prepare(target, level)
                with gzip.open(target, "rt", encoding="utf-8") as handle:
                    decompressed = handle.read()
                if reference is None:
                    reference = decompressed
                self.assertEqual(decompressed, reference, f"level {level}")
                self.assertEqual(
                    read_metrics(tmp / f"stats{level}.tsv")["retained_fragments"], "200"
                )
            # Level 0 must be honoured, not replaced by a truthiness default.
            self.assertGreater(
                (tmp / "level0" / "out.fq.gz").stat().st_size,
                (tmp / "level9" / "out.fq.gz").stat().st_size,
            )
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
