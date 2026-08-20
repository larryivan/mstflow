#!/usr/bin/env python3
"""Write docs/parameters.md from nextflow_schema.json. Re-run after editing the schema."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
schema = json.loads((ROOT / "nextflow_schema.json").read_text())


def cell(text):
    return str(text).replace("|", "\\|").replace("\n", " ")


def row(name, spec, required):
    kind = spec.get("type", "")
    kind = " or ".join(kind) if isinstance(kind, list) else kind
    if "enum" in spec:
        kind = " / ".join(f"`{v}`" for v in spec["enum"])
    default = spec.get("default", "")
    shown = json.dumps(default) if isinstance(default, bool) else default
    default = "required" if name in required else (f"`{shown}`" if default != "" else "")
    text = spec.get("description", "")
    if spec.get("help_text"):
        text += " " + spec["help_text"]
    return f"| `--{name}` | {cell(text)} | {kind} | {default} |"


lines = ["# Parameters", "",
         "Every option, generated from `nextflow_schema.json`. `nextflow run larryivan/mstflow --help`",
         "prints the same list; `--help <name>` explains one option. Options marked *advanced* are",
         "hidden from `--help` unless `--showHidden` is given.", ""]
for ref in schema["allOf"]:
    group = schema["$defs"][ref["$ref"].split("/")[-1]]
    required = set(group.get("required", []))
    props = group["properties"]
    shown = [(n, s) for n, s in props.items() if not s.get("hidden")]
    hidden = [(n, s) for n, s in props.items() if s.get("hidden")]
    lines += [f"## {group['title']}", ""]
    for title, items in (("", shown), ("Advanced", hidden)):
        if not items:
            continue
        if title:
            lines += [f"<details><summary>{title}</summary>", ""]
        lines += ["| Option | Description | Type | Default |", "|---|---|---|---|"]
        lines += [row(n, s, required) for n, s in items]
        lines += [""] + (["</details>", ""] if title else [])
(ROOT / "docs" / "parameters.md").write_text("\n".join(lines))
print("wrote docs/parameters.md")
