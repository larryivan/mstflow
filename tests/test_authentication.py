"""Authentication chooses the most abundant taxa of each ancient library, finds the
reference sequences below each in the taxonomy, and reports per taxon how many reads
aligned and how much of its genomes they cover, before testing it for damage."""
import csv
import gzip
import importlib.util
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


class SelectTargetsTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        write(self.root / "counts.tsv", [
            ["feature", "ancient", "modern"],
            ["taxid_1301__Streptococcus", "500", "900"],
            ["taxid_1386__Bacillus", "300", "0"],
            ["taxid_816__Bacteroides", "150", "0"],
            ["taxid_561__Escherichia", "50", "0"],
        ])
        write(self.root / "metadata.tsv", [["sample_id", "age_group"], ["ancient", "ancient"], ["modern", "modern"]])
        # Two Streptococcus species and one Bacillus species; 9606 is not in the taxonomy.
        (self.root / "nodes.dmp").write_text("".join(
            f"{child}\t|\t{parent}\t|\trank\t|\n"
            for child, parent in [(1, 1), (2, 1), (1301, 2), (1302, 1301), (1313, 1301), (1386, 2), (1423, 1386)]))
        with gzip.open(self.root / "genomes.fna.gz", "wt") as handle:
            for seqid, taxid in [("a", 1302), ("b", 1313), ("c", 1423), ("d", 9606)]:
                handle.write(f">{seqid}|kraken:taxid|{taxid} description\nACGT\n")

    def tearDown(self):
        self.directory.cleanup()

    def select(self, fasta="genomes.fna.gz"):
        return subprocess.run(
            [sys.executable, BIN / "select_authentication_targets.py", "--counts", "counts.tsv",
             "--metadata", "metadata.tsv", "--fasta", fasta, "--nodes", "nodes.dmp", "--top", "3", "--min-reads", "100"],
            cwd=self.root, capture_output=True, text=True)

    def test_top_taxa_of_ancient_libraries_with_their_genomes(self):
        self.assertEqual(self.select().returncode, 0)
        targets = [(row["sample_id"], row["taxid"], row["table_reads"], row["has_reference"])
                   for row in read(self.root / "targets.tsv")]
        self.assertEqual(targets, [
            ("ancient", "1301", "500", "true"),
            ("ancient", "1386", "300", "true"),
            ("ancient", "816", "150", "false"),
        ])
        references = {row["seqid"]: row["taxid"] for row in read(self.root / "references.tsv")}
        self.assertEqual(references, {"a|kraken:taxid|1302": "1301", "b|kraken:taxid|1313": "1301",
                                      "c|kraken:taxid|1423": "1386"})

    def test_header_without_a_taxon_is_rejected(self):
        (self.root / "plain.fa").write_text(">a\nACGT\n")
        result = self.select("plain.fa")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("kraken:taxid", result.stderr)


@unittest.skipUnless(importlib.util.find_spec("pysam") and importlib.util.find_spec("numpy"), "needs pysam and numpy")
class AuthenticateTests(unittest.TestCase):
    def test_reads_and_breadth_per_taxon(self):
        import pysam

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "a", "LN": 1000}, {"SN": "b", "LN": 1000}, {"SN": "c", "LN": 500}]}
            with pysam.AlignmentFile(str(root / "sample.bam"), "wb", header=header) as bam:
                for number, (reference, start) in enumerate([(0, 0), (0, 25), (0, 900), (2, 100)]):
                    segment = pysam.AlignedSegment(bam.header)
                    segment.query_name, segment.reference_id, segment.reference_start = f"r{number}", reference, start
                    segment.cigarstring, segment.query_sequence, segment.mapping_quality = "50M", "A" * 50, 60
                    segment.query_qualities = pysam.qualitystring_to_array("I" * 50)
                    bam.write(segment)
            write(root / "targets.tsv", [
                ["sample_id", "taxid", "name", "table_reads", "has_reference"],
                ["sample", "1301", "Streptococcus", "500", "true"],
                ["sample", "1386", "Bacillus", "300", "true"],
                ["sample", "816", "Bacteroides", "150", "false"],
                ["other", "1301", "Streptococcus", "900", "true"],
            ])
            write(root / "references.tsv", [["seqid", "taxid"], ["a", "1301"], ["b", "1301"], ["c", "1386"]])
            subprocess.run(
                [sys.executable, BIN / "authenticate_taxa.py", "--bam", "sample.bam", "--sample", "sample",
                 "--targets", "targets.tsv", "--references", "references.tsv", "--min-reads", "5", "--output", "out.tsv"],
                cwd=root, check=True)
            rows = {row["taxid"]: row for row in read(root / "out.tsv")}

        self.assertEqual(list(rows), ["1301", "1386", "816"])
        # Streptococcus: 75 + 50 bases of a's 1,000 covered, none of b's.
        self.assertEqual((rows["1301"]["aligned_reads"], rows["1301"]["breadth"]), ("3", "0.0625"))
        self.assertEqual((rows["1386"]["aligned_reads"], rows["1386"]["breadth"]), ("1", "0.1"))
        self.assertEqual({row["verdict"] for row in (rows["1301"], rows["1386"])}, {"too few reads"})
        self.assertEqual((rows["816"]["aligned_reads"], rows["816"]["verdict"]), ("NA", "no reference genome"))


if __name__ == "__main__":
    unittest.main()
