# MetaSourceFlow Workflow Notes

## 1. Read QC and Preprocessing

### 1.1 Modern Paired-End Samples

```bash
# 1. QC
fastqc data/raw/${SAMPLE}_R1.fastq.gz data/raw/${SAMPLE}_R2.fastq.gz -o logs/qc
multiqc logs/qc -o logs/qc

# 2. fastp trimming
fastp \
  -i data/raw/${SAMPLE}_R1.fastq.gz \
  -I data/raw/${SAMPLE}_R2.fastq.gz \
  -o results/trimmed/${SAMPLE}.trimmed_R1.fq.gz \
  -O results/trimmed/${SAMPLE}.trimmed_R2.fq.gz \
  --detect_adapter_for_pe \
  --cut_front --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --length_required 50 \
  --thread 16 \
  --html logs/qc/${SAMPLE}.fastp.html \
  --json logs/qc/${SAMPLE}.fastp.json
```

### 1.2 Ancient Paired-End Samples

```bash
# 1. QC
fastqc data/raw/${SAMPLE}_R1.fastq.gz data/raw/${SAMPLE}_R2.fastq.gz -o logs/qc
multiqc logs/qc -o logs/qc

# 2. AdapterRemoval trimming and collapse
AdapterRemoval \
  --file1 data/raw/${SAMPLE}_R1.fastq.gz \
  --file2 data/raw/${SAMPLE}_R2.fastq.gz \
  --basename results/trimmed/${SAMPLE} \
  --trimns \
  --trimqualities \
  --minquality 20 \
  --minlength 35 \
  --collapse \
  --gzip \
  --threads 16

# 3. Merge retained collapsed outputs
cat \
  results/trimmed/${SAMPLE}.collapsed.gz \
  results/trimmed/${SAMPLE}.collapsed.truncated.gz \
  > results/trimmed/${SAMPLE}.merged.fq.gz

# 4. CHM13 dehost with BWA aln
# 5. Kraken2 single-end profiling
```

## 2. Reference Data

### 2.1 Kraken2 PlusPF Database

```bash
aria2c \
  -x 16 -s 16 -k 1M \
  --continue=true \
  --file-allocation=none \
  --max-tries=0 \
  --retry-wait=10 \
  --timeout=60 \
  --connect-timeout=30 \
  -o k2_pluspf_20260226.tar.gz \
  https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20260226.tar.gz
```

### 2.2 Human Reference: T2T-CHM13

Download the CHM13 reference and index sidecar files:

```bash
mkdir -p ~/db/chm13
cd ~/db/chm13

aria2c -x 16 -s 16 -k 1M \
  --continue=true \
  --file-allocation=none \
  --max-tries=0 \
  --retry-wait=10 \
  -o chm13v2.0_maskedY.rCRS.fa.gz \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY.rCRS.fa.gz

aria2c -x 8 -s 8 -k 1M \
  --continue=true \
  --file-allocation=none \
  -o chm13v2.0_maskedY.rCRS.fa.gz.fai \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY.rCRS.fa.gz.fai

aria2c -x 8 -s 8 -k 1M \
  --continue=true \
  --file-allocation=none \
  -o chm13v2.0_maskedY.rCRS.fa.gz.gzi \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY.rCRS.fa.gz.gzi
```

Index the uncompressed FASTA:

```bash
gunzip -k chm13v2.0_maskedY.rCRS.fa.gz
bwa index chm13v2.0_maskedY.rCRS.fa
samtools faidx chm13v2.0_maskedY.rCRS.fa
```

## 3. Human Host Removal

Human-read removal was performed by aligning reads to T2T-CHM13. Ancient samples use BWA-backtrack (`bwa aln/samse`) with ancient-DNA-oriented parameters: seed disabled (`-l 1024`), relaxed edit-distance control (`-n 0.01`), and up to two gap opens (`-o 2`). These settings follow common aDNA mapping practice, where short and damaged molecules may fail seed-based alignment.

### 3.1 Ancient Single-End or Collapsed Reads

```bash
REF=~/db/chm13/chm13v2.0_maskedY.rCRS.fa
SAMPLE=sample_ancient
FASTQ=data/${SAMPLE}.fq.gz

mkdir -p results/dehost

bwa aln \
  -l 1024 \
  -n 0.01 \
  -o 2 \
  -t 32 \
  ${REF} \
  ${FASTQ} \
  > results/dehost/${SAMPLE}.chm13.sai

bwa samse \
  ${REF} \
  results/dehost/${SAMPLE}.chm13.sai \
  ${FASTQ} \
  | samtools view -@ 8 -b -f 4 -o results/dehost/${SAMPLE}.unmapped.bam -

samtools fastq \
  -@ 8 \
  results/dehost/${SAMPLE}.unmapped.bam \
  | gzip > results/dehost/${SAMPLE}.dehosted.fq.gz
```

### 3.2 Modern Single-End Samples

```bash
REF=~/db/chm13/chm13v2.0_maskedY.rCRS.fa
SAMPLE=sample_modern
FASTQ=data/${SAMPLE}.fq.gz

mkdir -p results/dehost

bwa mem \
  -t 32 \
  ${REF} \
  ${FASTQ} \
  | samtools view -@ 8 -b -f 4 -o results/dehost/${SAMPLE}.unmapped.bam -

samtools fastq \
  -@ 8 \
  results/dehost/${SAMPLE}.unmapped.bam \
  | gzip > results/dehost/${SAMPLE}.dehosted.fq.gz
```

### 3.3 Modern Paired-End Samples

```bash
REF=~/db/chm13/chm13v2.0_maskedY.rCRS.fa
SAMPLE=sample_modern
R1=data/${SAMPLE}_R1.fq.gz
R2=data/${SAMPLE}_R2.fq.gz

mkdir -p results/dehost

bwa mem \
  -t 32 \
  ${REF} \
  ${R1} \
  ${R2} \
  | samtools view -@ 8 -b -f 12 -F 256 -o results/dehost/${SAMPLE}.both_unmapped.bam -

samtools fastq \
  -@ 8 \
  -1 results/dehost/${SAMPLE}.dehosted_R1.fq.gz \
  -2 results/dehost/${SAMPLE}.dehosted_R2.fq.gz \
  -0 /dev/null \
  -s /dev/null \
  -n \
  results/dehost/${SAMPLE}.both_unmapped.bam
```

## 4. Taxonomic Profiling

### 4.1 Kraken2

```bash
DB=~/db/kraken2/k2_pluspf_20260226
SAMPLE=sample

kraken2 \
  --db ${DB} \
  --threads 32 \
  --gzip-compressed \
  --use-names \
  --report-minimizer-data \
  --confidence 0.05 \
  --minimum-hit-groups 2 \
  --report results/kraken2/${SAMPLE}.k2.report \
  --output results/kraken2/${SAMPLE}.k2.out \
  results/dehost/${SAMPLE}.dehosted.fq.gz
```
