"""The distinct k-mer filter in build_feature_table.py removes, sample by sample, taxa
whose reads cover too few distinct minimizers (Kraken2) or unique k-mers (KrakenUniq),
and reports what it removed. KrakenUniq reports can replace Bracken tables."""
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
        counts = {row["#OTU ID"] if "#OTU ID" in row else next(iter(row.values())): row
                  for row in read(self.root / "t_counts.tsv")}
        strep = next(row for key, row in counts.items() if "Streptococcus" in key)
        self.assertEqual((strep["s1"], strep["s2"]), ("500", "0"))
        qc = {row["sample_id"]: row for row in read(self.root / "t_none_sample_qc.tsv")}
        self.assertEqual((qc["s1"]["kmer_filter_removed_reads"], qc["s2"]["kmer_filter_removed_reads"]),
                         ("0", "500"))

    def test_reports_are_required(self):
        failed = self.build("--min-distinct-kmers", "100")
        self.assertNotEqual(failed.returncode, 0)
        self.assertIn("needs --minimizer-reports", failed.stderr)



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

    def test_genus_counts_totals_and_filters(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write(root / "metadata.tsv", [["sample_id"], ["s1"], ["s2"]])
            for sample, strep in (("s1", 20000), ("s2", 300)):
                (root / f"{sample}.krakenuniq.report.txt").write_text(self.REPORT.format(strep=strep))
            done = subprocess.run(
                [sys.executable, str(BIN / "build_feature_table.py"), "--metadata", "metadata.tsv", "--prefix", "u",
                 "--rank", "G", "--min-distinct-kmers", "1000", "--min-taxon-reads", "200",
                 "--krakenuniq-reports", "s1.krakenuniq.report.txt", "s2.krakenuniq.report.txt"],
                cwd=root, capture_output=True, text=True)
            self.assertEqual(done.returncode, 0, done.stderr)
            removed = {(r["sample_id"], r["name"]) for r in read(root / "u_kmer_filter.tsv")}
            # Bacillus has too few reads everywhere; Streptococcus too few k-mers in s2.
            self.assertEqual(removed, {("s1", "Bacillus"), ("s2", "Bacillus"), ("s2", "Streptococcus")})
            qc = {row["sample_id"]: row for row in read(root / "u_none_sample_qc.tsv")}
            self.assertEqual((qc["s1"]["kraken_classified_reads"], qc["s1"]["kraken_unclassified_reads"]), ("800", "200"))
            self.assertEqual(qc["s1"]["bracken_total_reads"], "500")

    def test_tables_and_reports_are_exclusive(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            write(root / "metadata.tsv", [["sample_id"], ["s1"]])
            failed = subprocess.run([sys.executable, str(BIN / "build_feature_table.py"), "--metadata", "metadata.tsv",
                                     "--prefix", "u"], cwd=root, capture_output=True, text=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("either Bracken tables or --krakenuniq-reports", failed.stderr)


if __name__ == "__main__":
    unittest.main()
