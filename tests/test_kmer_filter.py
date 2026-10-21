"""The distinct k-mer filter in build_feature_table.py removes, sample by sample, taxa
whose reads cover too few distinct minimizers (Kraken2) or unique k-mers (KrakenUniq),
and reports what it removed. Without Bracken, counts come from either report's clade
counts at the chosen rank."""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin"
BRACKEN_HEADER = ["name", "taxonomy_id", "taxonomy_lvl", "kraken_assigned_reads", "added_reads",
                  "new_est_reads", "fraction_total_reads"]


def write(path, rows):
    with open(path, "w", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Kraken2MinimizerFilterTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        write(self.root / "metadata.tsv", [["sample_id"], ["s1"], ["s2"]])
        # Both samples have 500 Streptococcus reads; in s2 they pile on 40 minimizers.
        for sample, strep_distinct in (("s1", 5000), ("s2", 40)):
            write(self.root / f"{sample}.bracken.tsv", [
                BRACKEN_HEADER,
                ["Streptococcus", "1301", "G", "500", "0", "500", "0.5"],
                ["Bacillus", "1386", "G", "500", "0", "500", "0.5"],
            ])
            (self.root / f"{sample}.kraken2.report.txt").write_text(
                "0.00\t0\t0\tU\t0\tunclassified\n100.00\t1000\t0\tR\t1\troot\n")
            (self.root / f"{sample}.kraken2.minimizer.report.txt").write_text(
                "100.00\t1000\t0\t90000\t9000\tR\t1\troot\n"
                f"50.00\t500\t500\t8000\t{strep_distinct}\tG\t1301\t  Streptococcus\n"
                "50.00\t500\t500\t8000\t4000\tG\t1386\t  Bacillus\n")

    def tearDown(self):
        self.directory.cleanup()

    def build(self, *extra):
        names = ("s1", "s2")
        return subprocess.run(
            [sys.executable, str(BIN / "build_feature_table.py"), "--metadata", "metadata.tsv", "--prefix", "t",
             *[f"{s}.bracken.tsv" for s in names], "--kraken-reports",
             *[f"{s}.kraken2.report.txt" for s in names], *extra],
            cwd=self.root, capture_output=True, text=True)

    def test_off_by_default(self):
        self.assertEqual(self.build().returncode, 0)
        self.assertEqual(read(self.root / "t_kmer_filter.tsv"), [])
        qc = {row["sample_id"]: row for row in read(self.root / "t_none_sample_qc.tsv")}
        self.assertEqual(qc["s2"]["kmer_filter_removed_reads"], "NA")

    def test_poorly_covered_taxa_are_removed_per_sample(self):
        done = self.build("--min-distinct-kmers", "100", "--minimizer-reports",
                          "s1.kraken2.minimizer.report.txt", "s2.kraken2.minimizer.report.txt")
        self.assertEqual(done.returncode, 0, done.stderr)
        removed = read(self.root / "t_kmer_filter.tsv")
        self.assertEqual([(r["sample_id"], r["name"], r["removed_reads"], r["distinct_kmers"]) for r in removed],
                         [("s2", "Streptococcus", "500", "40")])
        counts = {row["#OTU ID"]: row for row in read(self.root / "t_counts.tsv")}
        strep = counts["taxid_1301__Streptococcus"]
        self.assertEqual((strep["s1"], strep["s2"]), ("500", "0"))
        qc = {row["sample_id"]: row for row in read(self.root / "t_none_sample_qc.tsv")}
        self.assertEqual((qc["s1"]["kmer_filter_removed_reads"], qc["s2"]["kmer_filter_removed_reads"]),
                         ("0", "500"))

    def test_reports_are_required(self):
        failed = self.build("--min-distinct-kmers", "100")
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("or, for Kraken2, --minimizer-reports", failed.stderr)



class KrakenUniqTests(unittest.TestCase):
    REPORT = (
        "# KrakenUniq v1.0.4 DATE:2026-10-22 DB:db\n"
        "%\treads\ttaxReads\tkmers\tdup\tcov\ttaxID\trank\ttaxName\n"
        "20\t200\t200\t0\t0\tNA\t0\tno rank\tunclassified\n"
        "80\t800\t0\t90000\t1.2\tNA\t1\tno rank\troot\n"
        "50\t500\t20\t{strep}\t1.1\tNA\t1301\tgenus\t      Streptococcus\n"
        "48\t480\t480\t{strep}\t1.1\tNA\t1305\tspecies\t        Streptococcus sanguinis\n"
        "15\t150\t150\t3000\t1.0\tNA\t1386\tgenus\t      Bacillus\n"
    )

    def build(self, root, *args):
        return subprocess.run([sys.executable, str(BIN / "build_feature_table.py"), "--metadata", "metadata.tsv",
                               "--prefix", "u", *args], cwd=root, capture_output=True, text=True)

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        write(self.root / "metadata.tsv", [["sample_id"], ["s1"], ["s2"]])
        for sample, strep in (("s1", 20000), ("s2", 300)):
            (self.root / f"{sample}.krakenuniq.report.txt").write_text(self.REPORT.format(strep=strep))
        self.reports = ["--kraken-reports", "s1.krakenuniq.report.txt", "s2.krakenuniq.report.txt"]

    def tearDown(self):
        self.directory.cleanup()

    def test_clade_counts_without_bracken(self):
        done = self.build(self.root, "--rank", "G", "--min-distinct-kmers", "1000", "--min-taxon-reads", "200",
                          *self.reports)
        self.assertEqual(done.returncode, 0, done.stderr)
        removed = {(r["sample_id"], r["name"]) for r in read(self.root / "u_kmer_filter.tsv")}
        # Bacillus has too few reads everywhere; Streptococcus too few k-mers in s2.
        self.assertEqual(removed, {("s1", "Bacillus"), ("s2", "Bacillus"), ("s2", "Streptococcus")})
        qc = {row["sample_id"]: row for row in read(self.root / "u_none_sample_qc.tsv")}
        self.assertEqual((qc["s1"]["kraken_classified_reads"], qc["s1"]["kraken_unclassified_reads"]), ("800", "200"))
        self.assertEqual(qc["s1"]["bracken_total_reads"], "500")

    def test_bracken_estimates_filtered_by_krakenuniq_kmers(self):
        # Bracken moves reads classified above genus down to it; the unique k-mers
        # still come from the KrakenUniq reports.
        for sample in ("s1", "s2"):
            write(self.root / f"{sample}.bracken.tsv", [
                BRACKEN_HEADER,
                ["Streptococcus", "1301", "G", "500", "120", "620", "0.8"],
                ["Bacillus", "1386", "G", "150", "30", "180", "0.2"],
            ])
        done = self.build(self.root, "--min-distinct-kmers", "1000", "s1.bracken.tsv", "s2.bracken.tsv", *self.reports)
        self.assertEqual(done.returncode, 0, done.stderr)
        removed = [(r["sample_id"], r["name"], r["removed_reads"]) for r in read(self.root / "u_kmer_filter.tsv")]
        self.assertEqual(removed, [("s2", "Streptococcus", "620")])
        counts = {row["#OTU ID"]: row for row in read(self.root / "u_counts.tsv")}
        self.assertEqual(counts["taxid_1301__Streptococcus"]["s1"], "620")

    def test_kraken2_reports_give_clade_counts(self):
        for sample in ("s1", "s2"):
            (self.root / f"{sample}.kraken2.report.txt").write_text(
                " 25.00\t250\t250\tU\t0\tunclassified\n"
                " 75.00\t750\t10\tR\t1\troot\n"
                " 60.00\t600\t40\tG\t1301\t            Streptococcus\n"
                " 56.00\t560\t560\tS\t1305\t              Streptococcus sanguinis\n")
        done = self.build(self.root, "--kraken-reports", "s1.kraken2.report.txt", "s2.kraken2.report.txt")
        self.assertEqual(done.returncode, 0, done.stderr)
        counts = {row["#OTU ID"]: row for row in read(self.root / "u_counts.tsv")}
        self.assertEqual(list(counts), ["taxid_1301__Streptococcus"])
        self.assertEqual(counts["taxid_1301__Streptococcus"]["s2"], "600")

    def test_needs_tables_or_reports(self):
        failed = self.build(self.root)
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("give Bracken tables, or --kraken-reports", failed.stderr)

if __name__ == "__main__":
    unittest.main()
