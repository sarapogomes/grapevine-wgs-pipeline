# Grapevine WGS Pipeline — SNP-Based Cultivar Identification

Bioinformatics pipeline for whole genome sequencing (WGS) of *Vitis vinifera* samples, including read trimming, alignment, variant calling, and cultivar identification through SNP profile comparison against the VIVC (Vitis International Variety Catalogue) database.

> **Note:** The files `vivc112.bed`, `vivc110snpsv5v12x.txt` and `vivc_snp_db.txt` are not included in this repository as they are third-party resources. Please contact the authors or the VIVC directly to obtain them.

---

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Project Structure](#project-structure)
- [Pipeline](#pipeline)
  - [00. Conda Environment](#00-conda-environment)
  - [01. Quality Trimming](#01-quality-trimming)
  - [02. Quality Control](#02-quality-control)
  - [03. Reference Mapping](#03-reference-mapping)
  - [04. Variant Calling](#04-variant-calling)
  - [05. Genotype Processing in Excel](#05-genotype-processing-in-excel)
  - [06. VIVC Cultivar Identification](#06-vivc-cultivar-identification)
- [Reference Genome](#reference-genome)
- [Notes](#notes)

---

## Overview

```
Raw FASTQ
    ↓  Trimmomatic (adapter removal + quality filtering)
Trimmed FASTQ
    ↓  FastQC (quality report)
QC Reports
    ↓  BWA-MEM2 (alignment to T2T reference)
Sorted BAM
    ↓  GATK HaplotypeCaller (targeted variant calling — 110 SNPs)
GVCF
    ↓  bcftools query (genotype extraction)
Genotype table (.txt)
    ↓  Excel processing (letter conversion + strand correction)
SNP profile (.xlsx)
    ↓  R script comparison vs VIVC database
Cultivar match results (.txt)
```

**Number of samples processed:** 15 *Vitis vinifera* cultivar accessions

---

## Requirements

### System
- Windows 10/11 with WSL2 (Ubuntu)
- Miniconda or Anaconda

### Bioinformatics tools (via conda)

```bash
conda install -y \
    gatk4 \
    bwa \
    samtools \
    picard \
    trim-galore \
    fastqc \
    bcftools \
    snpeff \
    bedtools \
    vcftools \
    tabix \
    tree
```

> **Note:** BWA-MEM2 must be installed separately (not via conda in this setup):
> ```bash
> conda install -c bioconda bwa-mem2
> ```

### R packages (install once in RStudio)

```r
install.packages("readxl")
install.packages("openxlsx")
```

### Software versions used

| Tool | Version |
|---|---|
| Trimmomatic | 0.39 |
| FastQC | 0.12 |
| BWA-MEM2 | 2.2.1+ |
| SAMtools | 1.17+ |
| GATK | 4.6.2.0 |
| BCFtools | 1.17+ |
| R | 4.5.1 |

---

## Project Structure

```
wgs_videira/
│
├── raw_data/                        # Raw paired-end FASTQ files
│   ├── *_1.fastq
│   └── *_2.fastq
│
├── trimmed/                         # Trimmed reads (global output)
│   ├── *_1_paired.fastq
│   └── *_2_paired.fastq
│
├── qc/                              # FastQC reports on trimmed reads
│   ├── *_fastqc.html
│   └── *_fastqc.zip
│
└── vitis_project/
    │
    ├── scripts/                     # All pipeline scripts (this repository)
    │   ├── 00_ambiente_conda.txt
    │   ├── 01_trimming.txt
    │   ├── 02_qc.txt
    │   ├── 03_mapping_bwa-mem2.txt
    │   ├── 04_variant_calling.txt
    │   ├── 05_excel_.txt
    │   └── comparar_amostra_vivc_v4.R
    │
    ├── reference/
    │   ├── T2T_ref/                 # T2T reference genome (used for mapping)
    │   │   ├── T2T_ref.fasta
    │   │   ├── T2T_ref.fasta.fai
    │   │   ├── T2T_ref.dict
    │   │   └── [bwa-mem2 index files]
    │   ├── PN40024.v4/              # PN40024 v4 reference (alternative)
    │   ├── PN40024.12x_v2/          # PN40024 12X v2 (SNP panel origin)
    │   ├── vivc112.bed              # [NOT INCLUDED — third-party]
    │   ├── vivc110snpsv5v12x.txt    # [NOT INCLUDED — third-party]
    │   └── vivc_snp_db.txt          # [NOT INCLUDED — third-party]
    │
    ├── trimmed/                     # Trimmed reads (project copy)
    ├── mapping_V5/                  # Sorted BAM files + indexes
    ├── qc_bam_V5/                   # BAM statistics
    ├── vcf_files_V5/                # GVCF files + genotype tables (.txt)
    ├── genotypes/                   # Per-sample Excel genotype profiles (.xlsx)
    ├── results_vivc/                # VIVC comparison results
    └── resources/
        └── 72_und_40_SNPs_12Jan2021.xlsx  # Cabezas & Laucou SNP panels
```

---

## Pipeline

### 00. Conda Environment

```bash
# Create environment with Python 3.9
conda create -n wgs_analysis python=3.9

# Activate
conda activate wgs_analysis

# Add channels
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Install tools
conda install -y \
    gatk4 bwa samtools picard trim-galore \
    fastqc bcftools snpeff bedtools vcftools tabix tree

# Verify GATK
gatk --version
```

---

### 01. Quality Trimming

Trimmomatic 0.39 in paired-end mode. Parameters: adapter removal (TruSeq3-PE), quality sliding window (4:15), minimum length 36 bp.

```bash
java -jar /home/ser/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 \
    /mnt/c/Users/User/Desktop/wgs_videira/raw_data/SAMPLE_1.fastq \
    /mnt/c/Users/User/Desktop/wgs_videira/raw_data/SAMPLE_2.fastq \
    /mnt/c/Users/User/Desktop/wgs_videira/trimmed/SAMPLE_1_paired.fastq \
    /mnt/c/Users/User/Desktop/wgs_videira/trimmed/SAMPLE_1_unpaired.fastq \
    /mnt/c/Users/User/Desktop/wgs_videira/trimmed/SAMPLE_2_paired.fastq \
    /mnt/c/Users/User/Desktop/wgs_videira/trimmed/SAMPLE_2_unpaired.fastq \
    ILLUMINACLIP:/home/ser/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

Replace `SAMPLE` with the sample name for each accession.

---

### 02. Quality Control

FastQC on trimmed paired reads.

```bash
fastqc -o /mnt/c/Users/User/Desktop/wgs_videira/qc \
    "/mnt/c/Users/User/Desktop/wgs_videira/trimmed/SAMPLE_1_paired.fastq" \
    "/mnt/c/Users/User/Desktop/wgs_videira/trimmed/SAMPLE_2_paired.fastq"
```

---

### 03. Reference Mapping

BWA-MEM2 alignment to the T2T reference genome.

**Index the reference once:**

```bash
REF="/mnt/c/Users/User/Desktop/wgs_videira/vitis_project/reference/T2T_ref/T2T_ref.fasta"

bwa-mem2 index "$REF"
samtools faidx "$REF"
echo "Reference indexed!"
```

**Map each sample:**

```bash
cd /mnt/c/Users/User/Desktop/wgs_videira/vitis_project

REF="/mnt/c/Users/User/Desktop/wgs_videira/vitis_project/reference/T2T_ref/T2T_ref.fasta"
R1="/mnt/c/Users/User/Desktop/wgs_videira/vitis_project/trimmed/SAMPLE_1_paired.fastq"
R2="/mnt/c/Users/User/Desktop/wgs_videira/vitis_project/trimmed/SAMPLE_2_paired.fastq"

# Align and sort to BAM (pipe — no intermediate SAM file)
bwa-mem2 mem -t $(nproc) \
    -R "@RG\tID:SAMPLE\tSM:SAMPLE\tPL:ILLUMINA" \
    "$REF" "$R1" "$R2" | \
samtools view -@ $(nproc) -bS | \
samtools sort -@ $(nproc) -o SAMPLE_sorted.bam -

# Index BAM
samtools index SAMPLE_sorted.bam

echo "Mapping complete → SAMPLE_sorted.bam"
```

**BAM quality control:**

```bash
BAM="/mnt/c/Users/User/Desktop/wgs_videira/vitis_project/mapping_V5/SAMPLE_sorted.bam"

samtools flagstat "$BAM"
samtools stats "$BAM" > "${BAM%.bam}.stats.txt"
samtools depth -a "$BAM" | awk '{sum+=$3; n++} END {if(n>0) print "Mean coverage: " sum/n "x"}'
```

**Expected QC thresholds:**
- Mapped reads: > 95%
- Properly paired: > 80%
- Mean coverage: > 15x (30x recommended)

---

### 04. Variant Calling

GATK HaplotypeCaller in GVCF mode, restricted to the 110 SNP positions defined in `vivc112.bed`.

> **Requires:** `vivc112.bed` (not included — third-party file)

```bash
cd /mnt/c/Users/User/Downloads/gatk-4.6.2.0/gatk-4.6.2.0

./gatk HaplotypeCaller \
    -R /mnt/c/Users/User/Desktop/wgs_videira/vitis_project/reference/T2T_ref/T2T_ref.fasta \
    -I /mnt/c/Users/User/Desktop/wgs_videira/vitis_project/mapping_V5/SAMPLE_sorted.bam \
    -L /mnt/c/Users/User/Desktop/wgs_videira/vitis_project/vivc112.bed \
    -O /mnt/c/Users/User/Desktop/wgs_videira/vitis_project/vcf_files_V5/SAMPLE.g.vcf.gz \
    -ERC GVCF
```

**Extract genotypes:**

```bash
bcftools query \
    -f '%CHROM %POS %REF %ALT [%GT]\n' \
    /mnt/c/Users/User/Desktop/wgs_videira/vitis_project/vcf_files_V5/SAMPLE.g.vcf.gz \
    > /mnt/c/Users/User/Desktop/wgs_videira/vitis_project/vcf_files_V5/SAMPLE.txt
```

The output `.txt` file has columns: `CHROM POS REF ALT GT`

---

### 05. Genotype Processing in Excel

Import each sample's `.txt` file into Excel and build a genotype profile. Each sample's `.xlsx` contains two sheets:

**Sheet 1 — `SAMPLE_NAME`**

| Column | Description |
|---|---|
| CHROM | Chromosome |
| POS | Position (T2T coordinates) |
| REF | Reference allele |
| ALT | Alternative allele |
| GT | Numeric genotype (0/0, 0/1, 1/1) |
| GT_L | Genotype in letters with slash (e.g. `A/G`) |
| GT_CORRECT_2 | Strand-corrected genotype, no slash (e.g. `AG`) |
| SNP_ID | SNP name from VIVC panel |
| SUBJECT_STRAND | Strand orientation (`plus` or `minus`) |

**Sheet 2 — `vivc110snps`**

Reference panel with SNP positions in 12x and T2T coordinates.

**Excel formulas used:**

```excel
# Convert numeric GT to letter genotype (column GT_L)
=SE(F2="0/0";C2&"/"&C2;SE(F2="0/1";C2&"/"&D2;SE(F2="1/1";D2&"/"&D2;"")))

# Look up SNP name from reference panel
=ÍNDICE(vivc110snps!$A:$A;CORRESP(B2;vivc110snps!$C:$C;0))

# Look up strand
=ÍNDICE(vivc110snps!$D:$D;CORRESP(B2;vivc110snps!$C:$C;0))
```

**Strand correction for `minus` strand SNPs (column GT_CORRECT_2):**

SNPs with `SUBJECT_STRAND = minus` require complementing both alleles before comparison with the VIVC database:

```
A ↔ T
C ↔ G
```

Example: `GT_L = T/C` on minus strand → `GT_CORRECT_2 = AG`

For `plus` strand SNPs, `GT_CORRECT_2` is simply `GT_L` without the `/`.

---

### 06. VIVC Cultivar Identification

The R script `comparar_amostra_vivc_v4.R` compares each sample's corrected SNP profile against the VIVC database (~1874 varieties) using the 110 diagnostic SNPs.

> **Requires:** `vivc_snp_db.txt` (not included — third-party file)

**Usage in RStudio:**

1. Open `scripts/comparar_amostra_vivc_v4.R`
2. Edit only the configuration section at the top:

```r
rm(list = ls())  # clear environment before each run
library(readxl)

FICHEIRO_AMOSTRA <- "C:/path/to/genotypes/SAMPLE_NAME.xlsx"
NOME_AMOSTRA     <- "SAMPLE_NAME"   # must match Excel sheet name
FICHEIRO_DB      <- "C:/path/to/reference/vivc_snp_db.txt"
PASTA_OUTPUT     <- "C:/path/to/results_vivc/"
```

3. Run with `Ctrl+Shift+Enter`

**Output:** `SAMPLE_NAME_results.txt` — tab-separated file sorted by `score` descending:

| Column | Description |
|---|---|
| `score` | Number of SNPs matching the database variety |
| `md_count` | Number of SNPs with missing data |
| `pct_md` | Proportion of missing data (md_count / 110) |
| `accsession` | Variety name in VIVC |
| `SNP*` | TRUE/FALSE per SNP — TRUE = match |

**Interpreting results:** A variety with a high `score` and low `md_count` is a strong candidate match. Exact cultivar identity (score = 110, md_count = 0) confirms identification.

---

## Reference Genome

The **T2T (Telomere-to-Telomere) PN40024** reference genome was used for all mapping and variant calling steps. It is available for download at [Grapedia](https://grapedia.org/files-download/).

> Shi X. et al. (2023). The complete reference genome for grapevine (*Vitis vinifera* L.) genetics and breeding. *Horticulture Research*, 10(5), uhad061. https://doi.org/10.1093/hr/uhad061

---

## Notes

- The pipeline runs in **WSL2** (Windows Subsystem for Linux 2) with conda
- Large files (BAMs, FASTQs, GVCFs) are not included due to size — add them to `.gitignore`
- For each new sample, only the `SAMPLE` name needs to be updated in the scripts
- The `vivc112.bed` BED file targets 110 SNP positions in T2T coordinates and is required for GATK `--intervals`
- SNP strand correction is critical — incorrect strand orientation leads to systematically wrong genotype calls during VIVC comparison
- The script must be run with `rm(list = ls())` at the start to avoid carrying over variables between samples in RStudio
