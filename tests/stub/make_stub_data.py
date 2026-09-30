#!/usr/bin/env python3
import gzip
import os
from pathlib import Path


# Wiring-only data for tests/stub/run.sh. Generated on demand and git-ignored.
HERE = Path(__file__).resolve().parent
TEST_DIR = HERE / "data"
READ_DIR = TEST_DIR / "reads"
REF_DIR = TEST_DIR / "refs"
KRAKEN_DIR = TEST_DIR / "kraken2_stub_db"


HOST = (
    "TACGATGCTAGCTACCGATCGTACGTTAGCCTAGGCTAACGATCGGATCCGTAACG"
    "TTGACCTAGCTAGGATCGTACCTGATCGTAGCTAGGCTAATCGGATCGTAGCTAACG"
)
NONHOST_A = "TTGACCGTACGATCGTTAACGATCGATTCGATGCAATTCGATCGTACCGATTAACCGTA"
NONHOST_B = "GCTTAGGATCCGATTAACCGGTTAGCATCGATCGTAGGCTAACCGATCGTTAGCAATCG"
NONHOST_C = "CGATTCGGAATTCGATCGGATTAACCGATCGTAGCTAGCTTACGATCGGTTACCGATC"


def revcomp(seq):
    table = str.maketrans("ACGT", "TGCA")
    return seq.translate(table)[::-1]


def write_fastq(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for name, seq in records:
            qual = "I" * len(seq)
            handle.write(f"@{name}\n{seq}\n+\n{qual}\n")


def write_metadata(path):
    rows = [
        ("modern_pe", READ_DIR / "modern_pe_R1.fastq.gz", READ_DIR / "modern_pe_R2.fastq.gz", "PE", "modern", "source", "human_fecal_gastrointestinal", "remove"),
        ("modern_se", READ_DIR / "modern_se.fastq.gz", "NA", "SE", "modern", "source", "extraction_library_blank", "skip"),
        ("ancient_pe", READ_DIR / "ancient_pe_R1.fastq.gz", READ_DIR / "ancient_pe_R2.fastq.gz", "PE", "ancient", "sink", "human_fecal_gastrointestinal", "remove"),
        ("ancient_se", READ_DIR / "ancient_se.fastq.gz", "NA", "SE", "ancient", "sink", "extraction_library_blank", "skip"),
    ]
    with open(path, "w") as handle:
        handle.write(
            "sample_id\tfastq_1\tfastq_2\tlayout\tage_group\tpreprocess_mode\thost_action\thost_group\thost_ref\t"
            "source_sink\tmodel_source_env\tsample_label\tproject\trun\n"
        )
        for sample_id, fastq_1, fastq_2, layout, age_group, source_sink, env, host_action in rows:
            host_group = "human" if host_action == "remove" else "none"
            host_ref = str(REF_DIR / "chm13_tiny.fa") if host_action == "remove" else "NA"
            handle.write(
                f"{sample_id}\t{fastq_1}\t{fastq_2}\t{layout}\t{age_group}\tauto\t{host_action}\t{host_group}\t{host_ref}\t"
                f"{source_sink}\t{env}\t{sample_id}\ttiny\t{sample_id}\n"
            )


def main():
    READ_DIR.mkdir(parents=True, exist_ok=True)
    REF_DIR.mkdir(parents=True, exist_ok=True)

    with open(REF_DIR / "chm13_tiny.fa", "w") as handle:
        handle.write(">chm13_tiny\n")
        handle.write(HOST + "\n")
    # Stub mode never reads these indexes, but production preflight deliberately
    # requires a complete BWA bundle before any sample task is scheduled.
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"):
        (REF_DIR / f"chm13_tiny.fa{suffix}").write_text("mstflow stub index\n")

    KRAKEN_DIR.mkdir(parents=True, exist_ok=True)
    for filename in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (KRAKEN_DIR / filename).write_text("mstflow stub database\n")
    (KRAKEN_DIR / "database75mers.kmer_distrib").write_text("mstflow stub distribution\n")
    krakenuniq = TEST_DIR / "krakenuniq_stub_db"
    krakenuniq.mkdir(parents=True, exist_ok=True)
    for filename in ("database.kdb", "database.idx", "taxDB"):
        (krakenuniq / filename).write_text("mstflow stub database\n")

    write_fastq(
        READ_DIR / "modern_pe_R1.fastq.gz",
        [("modern_pe_host/1", HOST[:80]), ("modern_pe_nonhost/1", NONHOST_A[:60])],
    )
    write_fastq(
        READ_DIR / "modern_pe_R2.fastq.gz",
        [("modern_pe_host/2", revcomp(HOST[20:100])), ("modern_pe_nonhost/2", revcomp(NONHOST_B[:60]))],
    )
    write_fastq(
        READ_DIR / "modern_se.fastq.gz",
        [("modern_se_host", HOST[:80]), ("modern_se_nonhost", NONHOST_B[:60])],
    )

    ancient_insert = NONHOST_C[:54]
    write_fastq(
        READ_DIR / "ancient_pe_R1.fastq.gz",
        [("ancient_pe_host/1", HOST[:54]), ("ancient_pe_nonhost/1", ancient_insert)],
    )
    write_fastq(
        READ_DIR / "ancient_pe_R2.fastq.gz",
        [("ancient_pe_host/2", revcomp(HOST[:54])), ("ancient_pe_nonhost/2", revcomp(ancient_insert))],
    )
    write_fastq(
        READ_DIR / "ancient_se.fastq.gz",
        [("ancient_se_host", HOST[:54]), ("ancient_se_nonhost", NONHOST_C[:54])],
    )

    write_metadata(TEST_DIR / "metadata.tsv")

    # SourcePredict reads the NCBI taxonomy; SourceID-NMF's code is staged from a
    # directory. Stub mode needs only the files to exist, and no network.
    taxonomy = TEST_DIR / "taxonomy"
    taxonomy.mkdir(parents=True, exist_ok=True)
    (taxonomy / "nodes.dmp").write_text("1\t|\t1\t|\tno rank\t|\n")
    (taxonomy / "names.dmp").write_text("1\t|\troot\t|\t\t|\tscientific name\t|\n")
    nmf = TEST_DIR / "sourceid_nmf"
    nmf.mkdir(parents=True, exist_ok=True)
    for name in ("SourceID-NMF.py", "admm.py", "data_clustering.py", "data_estimation.py", "data_initialization.py"):
        (nmf / name).write_text("# stub\n")

    print(f"Wrote tiny test data under {TEST_DIR}")
    print(f"Metadata: {TEST_DIR / 'metadata.tsv'}")
    print(f"Host FASTA: {REF_DIR / 'chm13_tiny.fa'}")


if __name__ == "__main__":
    main()
