#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


PROKARYOTE_ROOTS = {"2", "2157"}
ROOT_TAXIDS = {"0", "1"}


def read_metadata(path):
    with open(path, newline="") as handle:
        return [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")]


def sample_id(path):
    name = Path(path).name
    for suffix in (".bracken.tsv", ".bracken.txt", ".tsv", ".txt"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path).stem


def feature_id(taxid, name):
    return f"taxid_{taxid}__{name.strip().replace(chr(9), ' ')}"


def read_bracken(path):
    sample = sample_id(path)
    counts = {}
    features = {}
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"name", "taxonomy_id", "taxonomy_lvl", "new_est_reads"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{path} is missing Bracken column(s): {', '.join(sorted(missing))}")
        for row in reader:
            taxid = row["taxonomy_id"].strip()
            name = row["name"].strip()
            fid = feature_id(taxid, name)
            counts[fid] = counts.get(fid, 0.0) + float(row["new_est_reads"] or 0)
            features[fid] = {"feature_id": fid, "taxid": taxid, "name": name, "rank": row["taxonomy_lvl"].strip()}
    return sample, counts, features


def read_nodes(taxdump):
    if not taxdump:
        raise SystemExit("--taxdump is required when --taxon-filter prokaryote")
    nodes = Path(taxdump)
    nodes = nodes / "nodes.dmp" if nodes.is_dir() else nodes
    if not nodes.exists():
        raise SystemExit(f"nodes.dmp not found: {nodes}")
    parents = {}
    with open(nodes) as handle:
        for line in handle:
            parts = [part.strip() for part in line.split("|")]
            if len(parts) >= 2:
                parents[parts[0]] = parts[1]
    return parents


def lineage_has(taxid, targets, parents, cache):
    if taxid in cache:
        return cache[taxid]
    seen = set()
    current = taxid
    while current and current not in seen:
        if current in targets:
            cache[taxid] = True
            return True
        if current in ROOT_TAXIDS or current not in parents:
            break
        seen.add(current)
        current = parents[current]
    cache[taxid] = False
    return False


def select_features(features, mode, taxdump, drop_taxids):
    parents = read_nodes(taxdump) if mode == "prokaryote" else {}
    cache = {}
    reasons = {}
    for fid, feature in features.items():
        taxid = feature["taxid"]
        if taxid in drop_taxids:
            reasons[fid] = "drop_taxid"
        elif mode == "prokaryote" and not lineage_has(taxid, PROKARYOTE_ROOTS, parents, cache):
            reasons[fid] = "not_prokaryote"
        else:
            reasons[fid] = "kept"
    kept = [fid for fid, reason in reasons.items() if reason == "kept"]
    return sorted(kept), reasons


def write_matrix(path, samples, features, counts, relative=False):
    totals = {sample: sum(counts.get(fid, {}).get(sample, 0.0) for fid in features) for sample in samples}
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["#OTU ID", *samples])
        for fid in features:
            row = [fid]
            for sample in samples:
                value = counts.get(fid, {}).get(sample, 0.0)
                if relative:
                    value = value / totals[sample] if totals[sample] else 0.0
                row.append(f"{value:.12g}")
            writer.writerow(row)


def write_features(path, feature_ids, features):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["feature_id", "taxid", "name", "rank"], lineterminator="\n")
        writer.writeheader()
        for fid in feature_ids:
            writer.writerow(features[fid])


def write_qc(path, raw_features, filtered_features, reasons, mode):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["taxon_filter", mode])
        writer.writerow(["raw_features", len(raw_features)])
        writer.writerow(["filtered_features", len(filtered_features)])
        writer.writerow(["dropped_features", len(raw_features) - len(filtered_features)])
        for reason in sorted(set(reasons.values())):
            writer.writerow([f"features_{reason}", sum(value == reason for value in reasons.values())])


def main():
    parser = argparse.ArgumentParser(description="Build Bracken feature tables with optional lineage filtering.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--taxon-filter", choices=["none", "prokaryote"], default="none")
    parser.add_argument("--taxdump")
    parser.add_argument("--drop-taxids", default="0,1")
    parser.add_argument("bracken_files", nargs="+")
    args = parser.parse_args()

    metadata_samples = read_metadata(args.metadata)
    counts = {}
    features = {}
    observed = []
    for path in args.bracken_files:
        sample, sample_counts, sample_features = read_bracken(path)
        observed.append(sample)
        features.update(sample_features)
        for fid, count in sample_counts.items():
            feature_counts = counts.setdefault(fid, {})
            feature_counts[sample] = feature_counts.get(sample, 0.0) + count

    observed_set = set(observed)
    samples = [sample for sample in metadata_samples if sample in observed_set]
    samples.extend(sorted(observed_set.difference(samples)))
    raw_features = sorted(features)
    drop_taxids = {taxid.strip() for taxid in args.drop_taxids.split(",") if taxid.strip()}
    filtered_features, reasons = select_features(features, args.taxon_filter, args.taxdump, drop_taxids)

    write_matrix(f"{args.prefix}_counts.tsv", samples, raw_features, counts, relative=False)
    write_matrix(f"{args.prefix}_relative_abundance.tsv", samples, raw_features, counts, relative=True)
    write_features(f"{args.prefix}_features.tsv", raw_features, features)

    filtered_prefix = f"{args.prefix}_{args.taxon_filter}"
    write_matrix(f"{filtered_prefix}_counts.tsv", samples, filtered_features, counts, relative=False)
    write_matrix(f"{filtered_prefix}_relative_abundance.tsv", samples, filtered_features, counts, relative=True)
    write_features(f"{filtered_prefix}_features.tsv", filtered_features, features)
    write_qc(f"{filtered_prefix}_filter_qc.tsv", raw_features, filtered_features, reasons, args.taxon_filter)


if __name__ == "__main__":
    main()
