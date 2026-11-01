"""How MSTFlow reads --input, --sources and --panel: every problem is reported at once,
before any task starts, and a good sheet is written out with its defaults filled in."""
import csv
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


def references(root):
    """A Kraken2 database and a BWA-indexed host as far as the start-up checks look."""
    database = root / "db"
    database.mkdir()
    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (database / name).write_text("stub\n")
    host = root / "host.fa"
    host.write_text(">host\nACGT\n")
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        Path(f"{host}{suffix}").write_text("stub\n")
    reads = root / "reads"
    reads.mkdir()
    for name in ("a_1.fq.gz", "a_2.fq.gz", "b.fq.gz", "c.fq.gz"):
        (reads / name).write_text("")
    return ["--kraken2_db", database, "--host_fasta", host]


def nextflow(root, *arguments):
    completed = subprocess.run(
        [NEXTFLOW, "run", ROOT, *map(str, arguments)],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )
    return completed.returncode, completed.stdout + completed.stderr


@unittest.skipUnless(NEXTFLOW, "nextflow is required for input tests")
class InputTests(unittest.TestCase):
    def test_every_problem_is_reported_before_the_run_starts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = references(root)
            (root / "samples.csv").write_text(
                "sample_id,fastq_1,fastq_2,age_group,host_action,host_ref,host_group\n"
                "bad id,reads/a_1.fq.gz,reads/a_2.fq.gz,old,,,\n"
                "missing,reads/none.fq.gz,,ancient,,,\n"
                "skipped,reads/b.fq.gz,,modern,skip,host.fa,\n"
                "dog,reads/b.fq.gz,,ancient,remove,,dog\n"
                "shared,reads/b.fq.gz,,modern,,,\n"
            )
            (root / "sources.tsv").write_text(
                "sample_id\tfastq_1\tlayout\tage_group\tenvironment\tproject\n"
                "shared\treads/c.fq.gz\tPE\tmodern\toral\tP1\n"
                "blank\treads/c.fq.gz\t\tmodern\t\tP1\n"
                "unknown\treads/c.fq.gz\t\tmodern\tUnknown\tP2\n"
            )
            status, output = nextflow(root, "-preview", "--input", "samples.csv", "--sources", "sources.tsv", *options)
            self.assertNotEqual(status, 0)
            for problem in [
                "MSTFlow cannot start:",
                "--input samples.csv line 2: sample_id 'bad id' may use only",
                "--input samples.csv line 2: age_group must be ancient or modern, not 'old'",
                "--input samples.csv line 3: reads/none.fq.gz was not found",
                "--input samples.csv line 4: host_action is skip, so host_ref must be empty",
                "--input samples.csv line 5: host_group 'dog' needs its own host_ref",
                "--sources sources.tsv line 2: layout is PE but fastq_2 is empty",
                "--sources sources.tsv line 3: environment is empty",
                "--sources sources.tsv line 4: environment 'unknown' is reserved",
            ]:
                self.assertIn(problem, output)
            # No task started: the problems were found while reading the sheets.
            self.assertNotIn("[PROCESS", output)

    def test_sheets_in_the_old_form_or_without_their_columns_are_explained(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = references(root)
            (root / "old.tsv").write_text(
                "sample_id\tfastq_1\tage_group\tsource_sink\tmodel_source_env\n"
                "a\treads/b.fq.gz\tancient\tsink\tNA\n"
            )
            (root / "sources.tsv").write_text("sample_id\tfastq_1\tage_group\nb\treads/c.fq.gz\tmodern\n")
            status, output = nextflow(root, "-preview", "--input", "old.tsv", "--sources", "sources.tsv", *options)
            self.assertNotEqual(status, 0)
            self.assertIn("give the reference samples in their own sheet with --sources", output)
            self.assertIn("--sources sources.tsv: missing column(s) environment", output)

    def test_a_mistyped_option_stops_the_run(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = references(root)
            (root / "samples.tsv").write_text("sample_id\tfastq_1\tage_group\na\treads/b.fq.gz\tancient\n")
            status, output = nextflow(root, "-preview", "--input", "samples.tsv", "--attribution_method", "feast", *options)
            self.assertNotEqual(status, 0)
            self.assertIn("--attribution_method", output)

    def test_a_panel_must_hold_both_tables_and_not_reuse_names(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = references(root)
            (root / "samples.tsv").write_text("sample_id\tfastq_1\tage_group\nsoil_1\treads/b.fq.gz\tancient\n")
            panel = root / "panel"
            panel.mkdir()
            (panel / "sources.tsv").write_text("source_id\tenvironment\tproject\tlibraries\nsoil_1\tsoil\tP1\tsoil_1\n")
            status, output = nextflow(root, "-preview", "--input", "samples.tsv", "--panel", panel, *options)
            self.assertNotEqual(status, 0)
            self.assertIn("has no counts.tsv; a panel is the panel/ directory of an earlier run", output)

            (panel / "counts.tsv").write_text("#OTU ID\tsoil_1\ntaxid_1__A\t5\n")
            status, output = nextflow(root, "-preview", "--input", "samples.tsv", "--panel", panel, *options)
            self.assertNotEqual(status, 0)
            self.assertIn("the panel already has source(s) named soil_1", output)

            # With --sources, a run writes <outdir>/panel, so it cannot also read it.
            (root / "samples.tsv").write_text("sample_id\tfastq_1\tage_group\nbone_1\treads/b.fq.gz\tancient\n")
            (root / "sources.tsv").write_text("sample_id\tfastq_1\tage_group\tenvironment\ngut_1\treads/c.fq.gz\tmodern\tgut\n")
            status, output = nextflow(root, "-preview", "--input", "samples.tsv", "--sources", "sources.tsv",
                                      "--panel", "./panel/", "--outdir", ".", *options)
            self.assertNotEqual(status, 0)
            self.assertIn("would be overwritten by this run's panel", output)

    def test_good_sheets_are_summarised_and_written_with_their_defaults(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            options = references(root)
            # As a spreadsheet program saves CSV: a byte-order mark, CRLF line ends,
            # quotes around a cell holding a comma, and a trailing blank line.
            (root / "samples.csv").write_bytes(
                "﻿sample_id,fastq_1,fastq_2,age_group,notes\r\n"
                'calculus,reads/a_1.fq.gz,reads/a_2.fq.gz,Ancient,"burial 3, layer 2"\r\n'
                "sediment,reads/b.fq.gz,,ancient,\r\n\r\n".encode("utf-8")
            )
            (root / "sources.tsv").write_text(
                "sample_id\tfastq_1\tage_group\thost_action\tenvironment\tproject\n"
                "plaque\treads/c.fq.gz\tmodern\tremove\t Oral\tP1\n"
                "soil\treads/c.fq.gz\tmodern\tskip\tsoil\tP2\n"
            )
            status, output = nextflow(root, "-preview", "--input", "samples.csv", "--sources", "sources.tsv", *options)
            self.assertEqual(status, 0, output)
            self.assertIn("Samples to trace  2 (2 ancient, 0 modern)", output)
            self.assertRegex(output, r"oral\s+1\s+1")
            self.assertIn("Source attribution with sourcetracker2.", output)
            self.assertIn("for reuse with --panel", output)

            # Without sources the run only profiles, and says so.
            status, output = nextflow(
                root, "-stub-run", "--input", "samples.csv", "--run_kraken2", "false", "--remove_host", "false",
                "--outdir", "results", *options,
            )
            self.assertEqual(status, 0, output)
            self.assertIn("Source attribution is off", output)
            with open(root / "results/pipeline_info/samplesheet.tsv", newline="") as handle:
                rows = {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["calculus"]["age_group"], "ancient")
            self.assertEqual(rows["calculus"]["layout"], "PE")
            self.assertEqual(rows["calculus"]["preprocess_mode"], "adapterremoval_pe_collapse")
            self.assertEqual(rows["calculus"]["notes"], "burial 3, layer 2")
            self.assertEqual(rows["calculus"]["source_sink"], "sink")
            self.assertEqual(rows["sediment"]["fastq_2"], "NA")
            self.assertEqual(rows["sediment"]["host_action"], "remove")
            self.assertEqual(rows["sediment"]["host_ref"], "default")
            self.assertEqual(rows["sediment"]["unit_id"], "sediment")


if __name__ == "__main__":
    unittest.main()
