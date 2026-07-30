import csv
import json
import os
import shutil
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
BRACKEN_MODULE = REPO / "modules" / "bracken.nf"
SELECT_BRACKEN_MODULE = REPO / "modules" / "select_bracken_length.nf"
FRAGMENT_HELPER = REPO / "bin" / "mstflow_fragment_tools.py"
NEXTFLOW = shutil.which("nextflow")

NO_READS_ERROR = "Error: no reads found. Please check your Kraken report"
BRACKEN_HEADER = (
    "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\t"
    "new_est_reads\tfraction_total_reads\n"
)


def groovy_string(value):
    """Return a safely quoted Groovy string for a generated test workflow."""
    return json.dumps(str(value))


@unittest.skipUnless(NEXTFLOW, "nextflow is required for Bracken module contract tests")
class BrackenModuleContractTests(unittest.TestCase):
    def _run_case(self, mode, genus_reads):
        temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(temporary_directory.cleanup)
        root = Path(temporary_directory.name)

        fake_bin = root / "fake_bin"
        fake_bin.mkdir()
        fake_bracken = fake_bin / "bracken"
        fake_bracken.write_text(
            textwrap.dedent(
                f"""\
                #!/usr/bin/env python3
                import json
                import os
                import sys
                from pathlib import Path


                NO_READS_ERROR = {NO_READS_ERROR!r}


                def option(arguments, flag):
                    try:
                        return arguments[arguments.index(flag) + 1]
                    except (ValueError, IndexError):
                        print(f"fake bracken: missing {{flag}}", file=sys.stderr)
                        raise SystemExit(64)


                arguments = sys.argv[1:]
                if arguments in (["-v"], ["--version"]):
                    print("Bracken version 2.9")
                    raise SystemExit(0)
                if arguments in (["-h"], ["--help"]):
                    print("Bracken version 2.9")
                    raise SystemExit(0)

                Path(os.environ["FAKE_BRACKEN_ARGS_LOG"]).write_text(
                    json.dumps(arguments), encoding="utf-8"
                )
                mode = os.environ["FAKE_BRACKEN_MODE"]

                if mode == "ordinary_error":
                    print("Error: database corrupt", file=sys.stderr)
                    raise SystemExit(42)
                if mode == "no_reads_wrong_exit":
                    print(NO_READS_ERROR, file=sys.stderr)
                    raise SystemExit(42)
                if mode == "mixed_error":
                    print(NO_READS_ERROR, file=sys.stderr)
                    print("Error: database corrupt", file=sys.stderr)
                    raise SystemExit(1)
                if mode != "threshold_model":
                    print(f"fake bracken: unknown mode {{mode}}", file=sys.stderr)
                    raise SystemExit(64)

                report = Path(option(arguments, "-i"))
                threshold = int(option(arguments, "-t"))
                genus_counts = []
                for line in report.read_text(encoding="utf-8").splitlines():
                    fields = line.split("\\t")
                    if len(fields) >= 6 and fields[3] == "G":
                        genus_counts.append(int(fields[1]))
                genus_count = max(genus_counts, default=0)

                if genus_count < threshold:
                    print(f">> Checking report file: {{report}}", file=sys.stderr)
                    print(NO_READS_ERROR, file=sys.stderr)
                    raise SystemExit(1)

                output = Path(option(arguments, "-o"))
                output.write_text(
                    {BRACKEN_HEADER!r}
                    + f"Testgenus\\t123\\tG\\t{{genus_count}}\\t0\\t{{genus_count}}\\t1.0\\n",
                    encoding="utf-8",
                )
                bracken_report = Path(option(arguments, "-w"))
                bracken_report.write_text(
                    f"100.00\\t{{genus_count}}\\t{{genus_count}}\\tG\\t123\\t  Testgenus\\n",
                    encoding="utf-8",
                )
                """
            ),
            encoding="utf-8",
        )
        fake_bracken.chmod(0o755)

        database = root / "db"
        database.mkdir()
        (database / "database75mers.kmer_distrib").touch()

        stats = root / "sample.fragment_length.tsv"
        stats.write_text(
            "metric\tvalue\nmedian_bracken_read_length\t75\n",
            encoding="utf-8",
        )
        report = root / "sample.kraken2.report.txt"
        report_text = (
            f"0.00\t0\t0\tU\t0\tunclassified\n"
            f"100.00\t{genus_reads}\t{genus_reads}\tG\t123\t  Testgenus\n"
        )
        report.write_text(report_text, encoding="utf-8")

        results = root / "results"
        workflow = root / "main.nf"
        workflow.write_text(
            textwrap.dedent(
                f"""\
                nextflow.enable.dsl = 2

                params.outdir = {groovy_string(results)}
                params.bracken_extra_args = ''
                params.bracken_rank = 'G'
                params.bracken_threshold = 10
                params.bracken_read_len = 'auto'
                params.container_python_tools = ''
                params.container_bracken = ''

                include {{ SELECT_BRACKEN_LENGTH }} from {groovy_string(SELECT_BRACKEN_MODULE)}
                include {{ BRACKEN_PROFILE }} from {groovy_string(BRACKEN_MODULE)}

                workflow {{
                    selection_input = Channel.of(
                        tuple(
                            [id: 'sample', bracken_read_len: 'auto'],
                            file({groovy_string(stats)}),
                            file({groovy_string(database)})
                        )
                    )
                    fragment_helper = Channel.value(file({groovy_string(FRAGMENT_HELPER)}))
                    SELECT_BRACKEN_LENGTH(selection_input, fragment_helper)
                    bracken_input = SELECT_BRACKEN_LENGTH.out.selections.map {{
                        meta, selection ->
                            tuple(
                                meta,
                                file({groovy_string(report)}),
                                selection,
                                file({groovy_string(database)})
                            )
                    }}
                    BRACKEN_PROFILE(bracken_input)
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
                """
            ),
            encoding="utf-8",
        )

        args_log = root / "bracken.args.json"
        environment = os.environ.copy()
        environment["PATH"] = os.pathsep.join([str(fake_bin), environment.get("PATH", "")])
        environment["FAKE_BRACKEN_MODE"] = mode
        environment["FAKE_BRACKEN_ARGS_LOG"] = str(args_log)
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
        self.assertEqual(report.read_text(encoding="utf-8"), report_text)
        return root, completed, args_log

    def _assert_success(self, completed):
        self.assertEqual(
            completed.returncode,
            0,
            f"Nextflow failed:\nSTDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}",
        )

    def _read_status(self, root):
        status_path = root / "results" / "bracken" / "sample.bracken.status.tsv"
        with open(status_path, encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        return rows[0]

    def _read_invocation(self, args_log):
        arguments = json.loads(args_log.read_text(encoding="utf-8"))

        def value(flag):
            return arguments[arguments.index(flag) + 1]

        return arguments, value

    def _assert_no_published_bracken_result(self, root):
        publish_dir = root / "results" / "bracken"
        published = (
            [
                path
                for path in publish_dir.glob("sample.bracken*")
                if not path.name.endswith(".bracken.read_length_selection.tsv")
            ]
            if publish_dir.exists()
            else []
        )
        self.assertEqual(published, [])

    def _command_stderr(self, root):
        return "\n".join(
            path.read_text(encoding="utf-8", errors="replace")
            for path in (root / "work").rglob(".command.err")
        )

    def test_below_threshold_exact_error_becomes_explicit_zero_result(self):
        root, completed, args_log = self._run_case("threshold_model", genus_reads=9)
        self._assert_success(completed)

        output = root / "results" / "bracken" / "sample.bracken.tsv"
        self.assertEqual(output.read_text(encoding="utf-8"), BRACKEN_HEADER)
        report = root / "results" / "bracken" / "sample.bracken.report.txt"
        self.assertEqual(report.read_text(encoding="utf-8"), "0.00\t0\t0\tR\t1\troot\n")

        status = self._read_status(root)
        self.assertEqual(status["sample_id"], "sample")
        self.assertEqual(status["status"], "zero_estimates")
        self.assertEqual(status["reason"], "no_estimates_at_requested_rank_threshold")
        self.assertEqual(status["bracken_rank"], "G")
        self.assertEqual(status["threshold"], "10")
        self.assertEqual(status["selected_read_length"], "75")
        self.assertEqual(status["bracken_exit_status"], "1")

        _, argument = self._read_invocation(args_log)
        self.assertEqual(argument("-t"), "10")
        self.assertEqual(argument("-l"), "G")
        self.assertEqual(argument("-r"), "75")

    def test_threshold_ten_is_inclusive(self):
        root, completed, args_log = self._run_case("threshold_model", genus_reads=10)
        self._assert_success(completed)

        output = root / "results" / "bracken" / "sample.bracken.tsv"
        self.assertEqual(
            output.read_text(encoding="utf-8"),
            BRACKEN_HEADER + "Testgenus\t123\tG\t10\t0\t10\t1.0\n",
        )
        status = self._read_status(root)
        self.assertEqual(status["status"], "estimated")
        self.assertEqual(status["reason"], "abundance_estimated")
        self.assertEqual(status["threshold"], "10")
        self.assertEqual(status["selected_read_length"], "75")
        self.assertEqual(status["bracken_exit_status"], "0")

        _, argument = self._read_invocation(args_log)
        self.assertEqual(argument("-t"), "10")
        self.assertEqual(argument("-l"), "G")
        self.assertEqual(argument("-r"), "75")

    def test_unrelated_bracken_error_remains_fatal(self):
        root, completed, _ = self._run_case("ordinary_error", genus_reads=9)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("Error: database corrupt", self._command_stderr(root))
        self._assert_no_published_bracken_result(root)

    def test_exact_no_reads_error_with_unexpected_exit_remains_fatal(self):
        root, completed, _ = self._run_case("no_reads_wrong_exit", genus_reads=9)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn(NO_READS_ERROR, self._command_stderr(root))
        self._assert_no_published_bracken_result(root)

    def test_exact_no_reads_line_plus_another_error_remains_fatal(self):
        root, completed, _ = self._run_case("mixed_error", genus_reads=9)
        self.assertNotEqual(
            completed.returncode,
            0,
            "A no-reads line must not hide an additional fatal Bracken error",
        )
        command_stderr = self._command_stderr(root)
        self.assertIn(NO_READS_ERROR, command_stderr)
        self.assertIn("Error: database corrupt", command_stderr)
        self._assert_no_published_bracken_result(root)


if __name__ == "__main__":
    unittest.main()
