#!/usr/bin/env python3
import gzip
import os
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TEST_DIR = ROOT / "examples" / "test_data"
READ_DIR = TEST_DIR / "reads"
REF_DIR = TEST_DIR / "refs"


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
        ("modern_pe", READ_DIR / "modern_pe_R1.fastq.gz", READ_DIR / "modern_pe_R2.fastq.gz", "PE", "modern", "source", "human_fecal_gastrointestinal"),
        ("modern_se", READ_DIR / "modern_se.fastq.gz", "NA", "SE", "modern", "source", "extraction_library_blank"),
        ("ancient_pe", READ_DIR / "ancient_pe_R1.fastq.gz", READ_DIR / "ancient_pe_R2.fastq.gz", "PE", "ancient", "sink", "human_fecal_gastrointestinal"),
        ("ancient_se", READ_DIR / "ancient_se.fastq.gz", "NA", "SE", "ancient", "sink", "extraction_library_blank"),
    ]
    with open(path, "w") as handle:
        handle.write(
            "sample_id\tfastq_1\tfastq_2\tlayout\tage_group\tpreprocess_mode\thost_ref\t"
            "source_sink\tmodel_source_env\tsample_label\tproject\trun\n"
        )
        for sample_id, fastq_1, fastq_2, layout, age_group, source_sink, env in rows:
            handle.write(
                f"{sample_id}\t{fastq_1}\t{fastq_2}\t{layout}\t{age_group}\tauto\t{REF_DIR / 'chm13_tiny.fa'}\t"
                f"{source_sink}\t{env}\t{sample_id}\ttiny\t{sample_id}\n"
            )


def main():
    READ_DIR.mkdir(parents=True, exist_ok=True)
    REF_DIR.mkdir(parents=True, exist_ok=True)

    with open(REF_DIR / "chm13_tiny.fa", "w") as handle:
        handle.write(">chm13_tiny\n")
        handle.write(HOST + "\n")

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

    write_metadata(ROOT / "examples" / "metadata.tiny.tsv")

    print(f"Wrote tiny test data under {TEST_DIR}")
    print(f"Metadata: {ROOT / 'examples' / 'metadata.tiny.tsv'}")
    print(f"Host FASTA: {REF_DIR / 'chm13_tiny.fa'}")


if __name__ == "__main__":
    main()
