import re
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
        }
        for process, (filename, parameter) in expected.items():
            text = (ROOT / "modules" / filename).read_text(encoding="utf-8")
            with self.subTest(process=process):
                self.assertIn(f"process {process} {{", text)
                self.assertIn(f'container "${{params.{parameter}}}"', text)

    def test_container_defaults_are_digest_pinned(self):
        config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
        images = dict(re.findall(r"^\s*(container_\w+) = '([^']+)'", config, re.MULTILINE))
        self.assertGreater(len(images), 10)
        for parameter, image in images.items():
            with self.subTest(parameter=parameter):
                self.assertIn("@sha256:", image)


def environment_dependencies(path):
    dependencies, inside = [], False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("dependencies:"):
            inside = True
        elif inside and line.strip().startswith("- "):
            dependencies.append(line.strip()[2:])
    return dependencies


class CondaEnvironmentTests(unittest.TestCase):
    def test_every_environment_file_is_used(self):
        config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
        self.assertEqual(
            set(re.findall(r"envs/([\w.]+\.yml)", config)),
            {path.name for path in (ROOT / "envs").glob("*.yml")},
        )

    def test_every_dependency_is_pinned(self):
        for path in sorted((ROOT / "envs").glob("*.yml")):
            for dependency in environment_dependencies(path):
                with self.subTest(environment=path.name, dependency=dependency):
                    self.assertRegex(dependency, r"^[\w.+-]+=[\w.*+-]+$")

    def test_r_methods_share_one_environment(self):
        config = (ROOT / "nextflow.config").read_text(encoding="utf-8")
        self.assertRegex(
            config,
            r"withName: 'R_METHODS_INSTALL\|FEAST\|FASTST' \{\s*conda = \"\$\{projectDir\}/envs/attribution_r\.yml\"",
        )


@unittest.skipUnless(NEXTFLOW, "nextflow is required for profile tests")
class ExecutionProfileTests(unittest.TestCase):
    def test_only_conda_builds_the_r_methods_from_source(self):
        for profile in ("conda", "docker", "apptainer", "singularity"):
            with self.subTest(profile=profile):
                values = flat_config(profile)
                self.assertEqual(values["params.r_methods_from_source"], str(profile == "conda").lower())

    def test_profiles_enable_exactly_one_dependency_engine(self):
        engines = ("conda", "docker", "apptainer", "singularity")
        for profile in engines:
            with self.subTest(profile=profile):
                values = flat_config(profile)
                for engine in engines:
                    self.assertEqual(values.get(f"{engine}.enabled", "false"), str(engine == profile).lower())


if __name__ == "__main__":
    unittest.main()
