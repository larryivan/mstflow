"""Contract tests against a real Kraken2 installation.

Skipped unless a Kraken2 database is provided:

    MSTFLOW_TEST_KRAKEN2_DB=/path/to/db  python3 -m unittest tests.test_kraken2_real_contract

Optionally set ``MSTFLOW_TEST_KRAKEN2_FASTA`` to a FASTA whose sequences are in
the database; otherwise a ``library/`` FASTA inside the database is used.  These
tests assert behaviour that cannot be established from a mock: that Kraken2
silently ignores unknown long options, and that pairing a single-end fragment
with a sub-k-mer placeholder mate reproduces native single-end classification
exactly.
"""

import gzip
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
HELPER = REPO / "bin" / "mstflow_fragment_tools.py"
KRAKEN2 = shutil.which("kraken2")
DATABASE = os.environ.get("MSTFLOW_TEST_KRAKEN2_DB")
KMER_LENGTH = 35
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reference_sequences():
    explicit = os.environ.get("MSTFLOW_TEST_KRAKEN2_FASTA")
    candidates = [Path(explicit)] if explicit else []
    if not candidates and DATABASE:
        candidates = sorted(Path(DATABASE).glob("library/**/*.fna"))
    for candidate in candidates:
        if not candidate.is_file():
            continue
        sequences = []
        current = []
        for line in candidate.read_text().splitlines():
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                current = []
            else:
                current.append(line.strip())
        if current:
            sequences.append("".join(current))
        longest = max(sequences, key=len, default="")
        if len(longest) >= 4000:
            return longest
    return ""


def write_fastq(path, records):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


@unittest.skipUnless(KRAKEN2, "kraken2 is required for real Kraken2 contract tests")
@unittest.skipUnless(DATABASE, "set MSTFLOW_TEST_KRAKEN2_DB to run real Kraken2 contract tests")
class RealKraken2ContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reference = reference_sequences()
        if not cls.reference:
            raise unittest.SkipTest(
                "no database library FASTA found; set MSTFLOW_TEST_KRAKEN2_FASTA"
            )

    def classify(self, *arguments, confidence="0.05", hit_groups="2"):
        """Return {read_id: (status, taxid, kmer_string)} from a real run."""
        completed = subprocess.run(
            [
                KRAKEN2,
                "--db",
                DATABASE,
                "--threads",
                "1",
                "--confidence",
                confidence,
                "--minimum-hit-groups",
                hit_groups,
                "--report",
                os.devnull,
                *map(str, arguments),
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        calls = {}
        for line in completed.stdout.splitlines():
            fields = line.split("\t")
            if len(fields) >= 5 and fields[0] in {"C", "U"}:
                calls[fields[1]] = (fields[0], fields[2], fields[4])
        return completed, calls

    def test_kraken2_silently_ignores_unknown_long_options(self):
        """Documents why MSTFlow must not rely on --check-names."""
        with tempfile.TemporaryDirectory() as tmp:
            reads = Path(tmp) / "reads.fq"
            write_fastq(reads, [("read0", self.reference[0:80])])
            completed, calls = self.classify("--check-names", reads)
            self.assertEqual(
                completed.returncode,
                0,
                "Kraken2 is expected to keep classifying after an unknown option",
            )
            self.assertIn("Unknown option", completed.stderr)
            self.assertEqual(len(calls), 1, "classification still happened")

    def test_dummy_mate_reproduces_native_single_end_classification(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            records = [
                (f"collapsed{index}", self.reference[200 + index * 90:200 + index * 90 + 40 + index])
                for index in range(30)
            ]
            single = tmp / "se.fq.gz"
            write_fastq(single, records)

            _, native = self.classify("--gzip-compressed", single)
            self.assertEqual(len(native), len(records))
            if not any(status == "C" for status, _, _ in native.values()):
                self.skipTest("no read classified against this database; cannot compare calls")

            subprocess.run(
                [
                    sys.executable,
                    str(HELPER),
                    "pair-fragments",
                    "--single",
                    str(single),
                    "--output-r1",
                    str(tmp / "R1.fq.gz"),
                    "--output-r2",
                    str(tmp / "R2.fq.gz"),
                    "--stats",
                    str(tmp / "stats.tsv"),
                    "--minimum-usable-span",
                    str(KMER_LENGTH),
                    "--kmer-length",
                    str(KMER_LENGTH),
                ],
                check=True,
            )
            _, paired = self.classify(
                "--gzip-compressed", "--paired", tmp / "R1.fq.gz", tmp / "R2.fq.gz"
            )

            self.assertEqual(sorted(paired), sorted(native))
            for read_id, (status, taxid, kmers) in native.items():
                dummy_status, dummy_taxid, dummy_kmers = paired[read_id]
                self.assertEqual((status, taxid), (dummy_status, dummy_taxid), read_id)
                # The placeholder mate contributes no minimizer positions, so the
                # hit string is the native one followed by an empty mate block.
                self.assertEqual(dummy_kmers.rstrip(" |:"), kmers.rstrip(" |:"), read_id)

    def test_pair_fragments_preserves_native_paired_classification(self):
        """Unmerged PE fragments pass through the mixed-input helper unchanged."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            r1_records = []
            r2_records = []
            for index in range(20):
                length = 40 + index % 5
                r1_records.append(
                    (f"paired{index}/1", self.reference[200 + index * 80:200 + index * 80 + length])
                )
                r2_records.append(
                    (f"paired{index}/2", self.reference[2100 + index * 80:2100 + index * 80 + length])
                )
            r1 = tmp / "native_R1.fq.gz"
            r2 = tmp / "native_R2.fq.gz"
            write_fastq(r1, r1_records)
            write_fastq(r2, r2_records)

            native_run, native = self.classify(
                "--gzip-compressed", "--paired", r1, r2
            )
            self.assertEqual(native_run.returncode, 0, native_run.stderr)
            self.assertEqual(len(native), len(r1_records))

            subprocess.run(
                [
                    sys.executable,
                    str(HELPER),
                    "pair-fragments",
                    "--r1",
                    str(r1),
                    "--r2",
                    str(r2),
                    "--output-r1",
                    str(tmp / "helper_R1.fq.gz"),
                    "--output-r2",
                    str(tmp / "helper_R2.fq.gz"),
                    "--stats",
                    str(tmp / "stats.tsv"),
                    "--minimum-usable-span",
                    str(KMER_LENGTH),
                    "--kmer-length",
                    str(KMER_LENGTH),
                ],
                check=True,
            )
            helper_run, rebuilt = self.classify(
                "--gzip-compressed",
                "--paired",
                tmp / "helper_R1.fq.gz",
                tmp / "helper_R2.fq.gz",
            )
            self.assertEqual(helper_run.returncode, 0, helper_run.stderr)
            self.assertEqual(rebuilt, native)

    def test_n_spacer_concatenation_changes_paired_classification(self):
        """The representation MSTFlow replaced is not Kraken2-equivalent."""
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            mate1 = self.reference[1000:1037]
            mate2 = self.reference[1100:1137].translate(COMPLEMENT)[::-1]
            r1 = tmp / "R1.fq.gz"
            r2 = tmp / "R2.fq.gz"
            write_fastq(r1, [("frag", mate1)])
            write_fastq(r2, [("frag", mate2)])
            spacer = tmp / "spacer.fq.gz"
            write_fastq(spacer, [("frag", f"{mate1}{'N' * 100}{mate2}")])

            _, native = self.classify("--gzip-compressed", "--paired", r1, r2)
            _, concatenated = self.classify("--gzip-compressed", spacer)
            self.assertEqual(sorted(native), ["frag"])
            self.assertEqual(sorted(concatenated), ["frag"])
            # Only the concatenated form reports ambiguous minimizer spans, and
            # those inflate total_kmers and therefore required_score.
            self.assertNotIn("A:", native["frag"][2])
            self.assertIn("A:", concatenated["frag"][2])


if __name__ == "__main__":
    unittest.main()
