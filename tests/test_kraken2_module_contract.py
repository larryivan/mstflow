"""Kraken2 module contract tests.

Kraken2 2.1.6 parses long options in its Perl wrapper with Getopt::Long and
never inspects the parse result (scripts/kraken2), and its classifier binary
only understands short options (src/classify.cc uses plain getopt).  An
unrecognised long option such as ``--check-names`` therefore prints
``Unknown option: check-names`` on stderr and classification still succeeds with
exit status 0.  These tests pin the two consequences for MSTFlow:

1. the module must not pass ``--check-names``, and must submit paired data with
   native ``--paired``;
2. if Kraken2 ever reports an ignored option, the task must fail closed instead
   of publishing a profile produced under unintended settings.
"""

import gzip
import json
import os
import shutil
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
KRAKEN2_MODULE = REPO / "modules" / "kraken2.nf"
PREPARE_KRAKEN_MODULE = REPO / "modules" / "prepare_kraken_input.nf"
PREPARE_KRAKEN_SCRIPT = REPO / "bin" / "prepare_kraken_input.py"
FRAGMENT_HELPER = REPO / "bin" / "mstflow_fragment_tools.py"
NEXTFLOW = shutil.which("nextflow")


def groovy_string(value):
    return json.dumps(str(value))


def write_fastq(path, records):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


@unittest.skipUnless(NEXTFLOW, "nextflow is required for Kraken2 module contract tests")
class Kraken2ModuleContractTests(unittest.TestCase):
    def _run_case(self, mode, layout="PE"):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        root = Path(temporary_directory.name)

        fake_bin = root / "fake_bin"
        fake_bin.mkdir()
        fake_kraken2 = fake_bin / "kraken2"
        fake_kraken2.write_text(
            textwrap.dedent(
                """\
                #!/usr/bin/env python3
                import json
                import os
                import sys
                from pathlib import Path


                arguments = sys.argv[1:]
                if arguments and arguments[0] == "--version":
                    print("Kraken version 2.1.6")
                    raise SystemExit(0)

                Path(os.environ["FAKE_KRAKEN2_ARGS_LOG"]).write_text(
                    json.dumps(arguments), encoding="utf-8"
                )

                mode = os.environ["FAKE_KRAKEN2_MODE"]
                if mode == "unknown_option":
                    # Reproduce the Getopt::Long behaviour of the real wrapper:
                    # warn about the option, then classify anyway and exit 0.
                    print("Unknown option: check-names", file=sys.stderr)
                if mode == "fatal_error":
                    print("classify: Error: unable to open database", file=sys.stderr)
                    raise SystemExit(2)

                print("Loading database information... done.", file=sys.stderr)


                def option(flag):
                    return arguments[arguments.index(flag) + 1]


                report = Path(option("--report"))
                if "--report-minimizer-data" in arguments:
                    report.write_text(
                        "  0.00\\t0\\t0\\t0\\t0\\tU\\t0\\tunclassified\\n"
                        "100.00\\t1\\t0\\t1\\t1\\tR\\t1\\troot\\n",
                        encoding="utf-8",
                    )
                else:
                    report.write_text(
                        "  0.00\\t0\\t0\\tU\\t0\\tunclassified\\n"
                        "100.00\\t1\\t0\\tR\\t1\\troot\\n",
                        encoding="utf-8",
                    )
                print("1 sequences (0.00 Mbp) processed", file=sys.stderr)
                raise SystemExit(0)
                """
            ),
            encoding="utf-8",
        )
        fake_kraken2.chmod(0o755)

        database = root / "db"
        database.mkdir()
        for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
            (database / name).touch()

        reads = root / "reads"
        reads.mkdir()
        if layout == "PE":
            write_fastq(reads / "sample.prehost_R1.fq.gz", [("frag", "ACGT" * 15)])
            write_fastq(reads / "sample.prehost_R2.fq.gz", [("frag", "TGCA" * 15)])
            read_files = [
                reads / "sample.prehost_R1.fq.gz",
                reads / "sample.prehost_R2.fq.gz",
            ]
            age_group = "modern"
        elif layout in {"MIXED", "MIXED_DEHOSTED"}:
            stage = "prehost" if layout == "MIXED" else "dehosted"
            write_fastq(reads / f"sample.{stage}_se.fq.gz", [("collapsed", "ACGT" * 15)])
            write_fastq(reads / f"sample.{stage}_R1.fq.gz", [("unmerged/1", "ACGT" * 15)])
            write_fastq(reads / f"sample.{stage}_R2.fq.gz", [("unmerged/2", "TGCA" * 15)])
            read_files = [
                reads / f"sample.{stage}_se.fq.gz",
                reads / f"sample.{stage}_R1.fq.gz",
                reads / f"sample.{stage}_R2.fq.gz",
            ]
            age_group = "ancient"
        else:
            write_fastq(reads / "sample.prehost.fq.gz", [("frag", "ACGT" * 15)])
            read_files = [reads / "sample.prehost.fq.gz"]
            age_group = "ancient"

        results = root / "results"
        read_list = ", ".join(f"file({groovy_string(path)})" for path in read_files)
        workflow = root / "main.nf"
        workflow.write_text(
            textwrap.dedent(
                f"""\
                nextflow.enable.dsl = 2

                params.outdir = {groovy_string(results)}
                params.kraken2_extra_args = ''
                params.kraken2_kmer_len = 35
                params.kraken2_min_read_len = 35
                params.kraken2_confidence = 0.05
                params.kraken2_minimum_hit_groups = 2
                params.kraken2_use_names = true
                params.kraken2_report_minimizer_data = true
                params.kraken2_save_read_classifications = false
                params.fastq_gzip_level = 6
                params.container_python_tools = ''
                params.container_kraken2 = ''

                include {{ PREPARE_KRAKEN_INPUT }} from {groovy_string(PREPARE_KRAKEN_MODULE)}
                include {{ KRAKEN2_PROFILE }} from {groovy_string(KRAKEN2_MODULE)}

                workflow {{
                    biological_reads = Channel.of(
                        tuple(
                            [id: 'sample', age_group: {groovy_string(age_group)}],
                            [{read_list}]
                        )
                    )
                    prepare_script = Channel.value(file({groovy_string(PREPARE_KRAKEN_SCRIPT)}))
                    fragment_helper = Channel.value(file({groovy_string(FRAGMENT_HELPER)}))
                    PREPARE_KRAKEN_INPUT(biological_reads, prepare_script, fragment_helper)
                    kraken_input = PREPARE_KRAKEN_INPUT.out.prepared.map {{
                        meta, prepared_reads, read_length, input_mode ->
                            tuple(
                                meta,
                                prepared_reads,
                                read_length,
                                input_mode,
                                file({groovy_string(database)})
                            )
                    }}
                    KRAKEN2_PROFILE(kraken_input)
                }}
                """
            ),
            encoding="utf-8",
        )
        (root / "nextflow.config").write_text(
            textwrap.dedent(
                """\
                process.executor = 'local'
                process.shell = ['/bin/bash', '-euo', 'pipefail']
                process.cpus = 1
                """
            ),
            encoding="utf-8",
        )

        args_log = root / "kraken2.args.json"
        environment = os.environ.copy()
        environment["PATH"] = os.pathsep.join([str(fake_bin), environment.get("PATH", "")])
        environment["FAKE_KRAKEN2_MODE"] = mode
        environment["FAKE_KRAKEN2_ARGS_LOG"] = str(args_log)
        environment["NXF_ANSI_LOG"] = "false"
        completed = subprocess.run(
            [
                NEXTFLOW,
                "run",
                str(workflow),
                "-ansi-log",
                "false",
                "-work-dir",
                str(root / "work"),
            ],
            cwd=root,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )
        return root, completed, args_log

    def test_paired_input_uses_native_paired_without_check_names(self):
        root, completed, args_log = self._run_case("normal", layout="PE")
        self.assertEqual(
            completed.returncode,
            0,
            f"Nextflow failed:\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}",
        )
        arguments = json.loads(args_log.read_text(encoding="utf-8"))
        self.assertNotIn(
            "--check-names",
            arguments,
            "--check-names is not a Kraken2 2.1.6 option; it is silently ignored",
        )
        self.assertIn("--paired", arguments)
        # Exactly two positional read operands follow the flags.
        operands = [
            argument
            for argument in arguments
            if argument.endswith((".fq.gz", ".fastq.gz"))
        ]
        self.assertEqual(len(operands), 2, arguments)
        self.assertIn(".kraken_input_R1.", Path(operands[0]).name)
        self.assertIn(".kraken_input_R2.", Path(operands[1]).name)
        self.assertTrue((root / "results" / "kraken2" / "sample.kraken2.report.txt").is_file())

    def test_single_end_input_is_not_submitted_as_paired(self):
        _, completed, args_log = self._run_case("normal", layout="SE")
        self.assertEqual(completed.returncode, 0, completed.stderr)
        arguments = json.loads(args_log.read_text(encoding="utf-8"))
        self.assertNotIn("--paired", arguments)
        self.assertNotIn("--check-names", arguments)

    def test_ancient_mixed_input_uses_native_paired_preparation(self):
        for layout in ("MIXED", "MIXED_DEHOSTED"):
            with self.subTest(layout=layout):
                root, completed, args_log = self._run_case("normal", layout=layout)
                self.assertEqual(completed.returncode, 0, completed.stderr)
                arguments = json.loads(args_log.read_text(encoding="utf-8"))
                self.assertIn("--paired", arguments)
                self.assertNotIn("--check-names", arguments)

                qc_path = root / "results" / "kraken2" / "sample.kraken2.read_length.tsv"
                with open(qc_path, encoding="utf-8") as handle:
                    qc = dict(
                        line.rstrip("\n").split("\t", 1)
                        for line in handle
                        if not line.startswith("metric\t")
                    )
                self.assertEqual(qc["input_single_fragments"], "1")
                self.assertEqual(qc["input_paired_fragments"], "1")
                self.assertEqual(qc["retained_single_fragments"], "1")
                self.assertEqual(qc["retained_paired_fragments"], "1")
                self.assertEqual(qc["retained_fragments"], "2")

    def test_ignored_unknown_option_fails_the_task(self):
        root, completed, _ = self._run_case("unknown_option", layout="PE")
        self.assertNotEqual(
            completed.returncode,
            0,
            "Kraken2 exits 0 after ignoring an unknown option; MSTFlow must fail closed",
        )
        self.assertIn("Kraken2 ignored an unrecognised option", completed.stdout + completed.stderr)
        publish_dir = root / "results" / "kraken2"
        published = list(publish_dir.glob("sample.kraken2*")) if publish_dir.exists() else []
        self.assertEqual(published, [])

    def test_real_kraken2_failure_still_reports_its_stderr(self):
        """Capturing stderr for inspection must not hide a genuine failure."""
        root, completed, _ = self._run_case("fatal_error", layout="PE")
        self.assertNotEqual(completed.returncode, 0)
        command_stderr = "\n".join(
            path.read_text(encoding="utf-8", errors="replace")
            for path in (root / "work").rglob(".command.err")
        )
        self.assertIn("unable to open database", command_stderr)
        self.assertNotIn("ignored an unrecognised option", command_stderr)

    def test_mate_identity_is_verified_before_kraken2_runs(self):
        """Kraken2 2.1.6 cannot check mate names, so the helper must."""
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        root = Path(temporary_directory.name)
        r1 = root / "R1.fq.gz"
        r2 = root / "R2.fq.gz"
        write_fastq(r1, [("first/1", "ACGT" * 15)])
        write_fastq(r2, [("different/2", "TGCA" * 15)])
        completed = subprocess.run(
            [
                "python3",
                str(FRAGMENT_HELPER),
                "stats",
                "--r1",
                str(r1),
                "--r2",
                str(r2),
                "--output",
                str(root / "stats.tsv"),
                "--minimum-usable-span",
                "35",
                "--kmer-length",
                "35",
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("paired FASTQ name mismatch", completed.stderr)


if __name__ == "__main__":
    unittest.main()
