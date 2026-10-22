"""build_feature_table.py counts, for each taxon at the chosen rank, the reads a Kraken2
or KrakenUniq report places in it or below it, and removes, sample by sample, taxa whose
reads cover too few distinct minimizers (Kraken2) or unique k-mers (KrakenUniq),
reporting what it removed."""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin"


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
            (self.root / f"{sample}.kraken2.report.txt").write_text(
                "0.00\t0\t0\tU\t0\tunclassified\n100.00\t1000\t0\tR\t1\troot\n"
                "50.00\t500\t500\tG\t1301\t  Streptococcus\n50.00\t500\t500\tG\t1386\t  Bacillus\n")
            (self.root / f"{sample}.kraken2.minimizer.report.txt").write_text(
                "100.00\t1000\t0\t90000\t9000\tR\t1\troot\n"
                f"50.00\t500\t500\t8000\t{strep_distinct}\tG\t1301\t  Streptococcus\n"
                "50.00\t500\t500\t8000\t4000\tG\t1386\t  Bacillus\n")

    def tearDown(self):
        self.directory.cleanup()

    def build(self, *extra):
        return subprocess.run(
            [sys.executable, str(BIN / "build_feature_table.py"), "--metadata", "metadata.tsv", "--prefix", "t",
             "s1.kraken2.report.txt", "s2.kraken2.report.txt", *extra],
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

    def test_minimizer_reports_are_required(self):
        failed = self.build("--min-distinct-kmers", "100")
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("or, for Kraken2, --minimizer-reports", failed.stderr)


class CladeCountTests(unittest.TestCase):
    KRAKENUNIQ = (
        "# KrakenUniq v1.0.4 DATE:2026-10-22 DB:db\n"
        "%\treads\ttaxReads\tkmers\tdup\tcov\ttaxID\trank\ttaxName\n"
        "20\t200\t200\t0\t0\tNA\t0\tno rank\tunclassified\n"
        "80\t800\t0\t90000\t1.2\tNA\t1\tno rank\troot\n"
        "50\t500\t20\t{strep}\t1.1\tNA\t1301\tgenus\t      Streptococcus\n"
        "48\t480\t480\t{strep}\t1.1\tNA\t1305\tspecies\t        Streptococcus sanguinis\n"
        "15\t150\t150\t3000\t1.0\tNA\t1386\tgenus\t      Bacillus\n"
    )
    # Reads left at the family stay out of the genus table; species reads count for
    # their genus.
    KRAKEN2 = (
        " 25.00\t250\t250\tU\t0\tunclassified\n"
        " 75.00\t750\t10\tR\t1\troot\n"
        " 70.00\t700\t95\tF\t1300\t          Streptococcaceae\n"
        " 60.00\t600\t40\tG\t1301\t            Streptococcus\n"
        " 56.00\t560\t560\tS\t1305\t              Streptococcus sanguinis\n"
        " {lacto:.2f}\t{lacto}\t{lacto}\tG\t1357\t            Lactococcus\n"
    )

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        write(self.root / "metadata.tsv", [["sample_id"], ["s1"], ["s2"]])

    def tearDown(self):
        self.directory.cleanup()

    def build(self, *args):
        return subprocess.run([sys.executable, str(BIN / "build_feature_table.py"), "--metadata", "metadata.tsv",
                               "--prefix", "u", *args], cwd=self.root, capture_output=True, text=True)

    def test_krakenuniq_reads_and_unique_kmers(self):
        for sample, strep in (("s1", 20000), ("s2", 300)):
            (self.root / f"{sample}.krakenuniq.report.txt").write_text(self.KRAKENUNIQ.format(strep=strep))
        done = self.build("--rank", "G", "--min-distinct-kmers", "1000",
                          "s1.krakenuniq.report.txt", "s2.krakenuniq.report.txt")
        self.assertEqual(done.returncode, 0, done.stderr)
        removed = {(r["sample_id"], r["name"]) for r in read(self.root / "u_kmer_filter.tsv")}
        # Streptococcus covers too few unique k-mers in s2.
        self.assertEqual(removed, {("s2", "Streptococcus")})
        qc = {row["sample_id"]: row for row in read(self.root / "u_none_sample_qc.tsv")}
        self.assertEqual((qc["s1"]["kraken_classified_reads"], qc["s1"]["kraken_unclassified_reads"]), ("800", "200"))
        self.assertEqual(qc["s1"]["rank_total_reads"], "650")

    def test_kraken2_clade_counts_add_up(self):
        for sample, lacto in (("s1", 5), ("s2", 50)):
            (self.root / f"{sample}.kraken2.report.txt").write_text(self.KRAKEN2.format(lacto=lacto))
        done = self.build("s1.kraken2.report.txt", "s2.kraken2.report.txt")
        self.assertEqual(done.returncode, 0, done.stderr)
        counts = {row["#OTU ID"]: row for row in read(self.root / "u_counts.tsv")}
        self.assertEqual(list(counts), ["taxid_1301__Streptococcus", "taxid_1357__Lactococcus"])
        self.assertEqual(counts["taxid_1301__Streptococcus"]["s2"], "600")
        # Every read counts, however few: a mixture's counts are the sum of its parts'.
        self.assertEqual((counts["taxid_1357__Lactococcus"]["s1"], counts["taxid_1357__Lactococcus"]["s2"]), ("5", "50"))

    def test_reports_are_required(self):
        failed = self.build()
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("reports", failed.stderr)


if __name__ == "__main__":
    unittest.main()
