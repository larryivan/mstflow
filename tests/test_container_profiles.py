import shutil
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = shutil.which("nextflow")


def flat_config(profile):
    completed = subprocess.run(
        [NEXTFLOW, "config", "-profile", profile, "-flat"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise AssertionError(
            f"nextflow config failed for {profile}\n"
            f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
        )
    values = {}
    for line in completed.stdout.splitlines():
        if " = " not in line:
            continue
        key, value = line.split(" = ", 1)
        values[key] = value.strip("'")
    return values


class ContainerDirectiveTests(unittest.TestCase):
    def test_every_process_declares_one_container(self):
        for module in sorted((ROOT / "modules").glob("*.nf")):
            text = module.read_text(encoding="utf-8")
            with self.subTest(module=module.name):
                self.assertEqual(
                    text.count("process "),
                    text.count('container "${params.container_'),
                    f"{module.name} has a process without an explicit logical container",
                )

    def test_single_tool_processes_use_matching_images(self):
        expected = {
            "ADAPTERREMOVAL_RAW": ("preprocess_adapterremoval.nf", "container_adapterremoval"),
            "NORMALIZE_ADAPTERREMOVAL": ("normalize_adapterremoval.nf", "container_python_tools"),
            "HOST_ALIGN_BWA": ("host_remove_bwa.nf", "container_bwa"),
            "FILTER_HOST_ALIGNMENTS": ("filter_host_alignments.nf", "container_python_tools"),
            "PREPARE_KRAKEN_INPUT": ("prepare_kraken_input.nf", "container_python_tools"),
            "KRAKEN2_PROFILE": ("kraken2.nf", "container_kraken2"),
            "SELECT_BRACKEN_LENGTH": ("select_bracken_length.nf", "container_python_tools"),
            "BRACKEN_PROFILE": ("bracken.nf", "container_bracken"),
        }
        for process, (filename, parameter) in expected.items():
            text = (ROOT / "modules" / filename).read_text(encoding="utf-8")
            with self.subTest(process=process):
                self.assertIn(f"process {process} {{", text)
                self.assertIn(f'container "${{params.{parameter}}}"', text)

    def test_container_defaults_are_digest_pinned(self):
        config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
        for parameter in (
            "container_fastqc",
            "container_fastp",
            "container_adapterremoval",
            "container_bwa",
            "container_kraken2",
            "container_bracken",
            "container_multiqc",
            "container_python_tools",
            "container_sourcetracker2",
            "container_sourcepredict",
            "container_krakenuniq",
        ):
            with self.subTest(parameter=parameter):
                line = next(
                    candidate
                    for candidate in config.splitlines()
                    if candidate.strip().startswith(f"{parameter} = ")
                )
                self.assertIn("@sha256:", line)


@unittest.skipUnless(NEXTFLOW, "nextflow is required for profile tests")
class ExecutionProfileTests(unittest.TestCase):
    def test_profiles_enable_exactly_one_dependency_engine(self):
        expectations = {
            "standard": ("local", None),
            "conda": ("local", "conda"),
            "docker": ("local", "docker"),
            "apptainer": ("local", "apptainer"),
            "singularity": ("local", "singularity"),
        }
        engines = ("conda", "docker", "apptainer", "singularity")
        for profile, (executor, enabled_engine) in expectations.items():
            with self.subTest(profile=profile):
                values = flat_config(profile)
                self.assertEqual(values["process.executor"], executor)
                for engine in engines:
                    expected = "true" if engine == enabled_engine else "false"
                    self.assertEqual(values[f"{engine}.enabled"], expected)

    def test_no_scheduler_profiles_are_shipped(self):
        """Executors are site policy and are supplied with -c, never as a profile."""
        for profile in ("slurm", "slurm_conda", "slurm_apptainer", "slurm_singularity"):
            with self.subTest(profile=profile):
                result = subprocess.run(
                    [NEXTFLOW, "config", "-profile", profile, "-flat"],
                    cwd=ROOT, capture_output=True, text=True,
                )
                self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
