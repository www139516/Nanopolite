# NanoPolite

**NanoPolite** is an open-source Bash pipeline for automated processing of Oxford Nanopore Technology (ONT) multiplex amplicon sequencing data. It takes per-sample demultiplexed FASTQ files as input and produces polished, orientation-corrected, locus-assigned consensus sequences ready for downstream genotyping or variant analysis — without requiring a per-locus reference configuration.

NanoPolite was developed for cost-effective genotyping of maize inbred lines using multiplex PCR amplicon nanopore sequencing, and is described in:

> Wang Y. et al. (2026) *Cost-effective genotyping of maize inbred lines using multiplex amplicon nanopore sequencing.* Plant Methods. [DOI pending]

---

## Table of Contents

- [Overview](#overview)
- [Pipeline Stages](#pipeline-stages)
- [Output Files](#output-files)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Parameter Reference](#parameter-reference)
- [Input File Formats](#input-file-formats)
- [Example Commands](#example-commands)
- [Citation](#citation)
- [License](#license)

---

## Overview

NanoPolite implements four sequential processing stages:

```
FASTQ input(s)
    │
    ▼
[QC]  NanoFilt quality filtering
    │
    ▼
[Cluster]  vsearch read clustering
    │
    ▼
[Polish]  minimap2 + Racon (or Medaka) consensus polishing
    │
    ├──► 1_Total_Clusters.fasta
    └──► 2_Polished_Clusters.fasta
         │
         ▼
[Stage A]  BLASTn-based strand orientation (optional, requires -r)
         │
         ▼
[Stage B]  Primer-based amplicon trimming & locus assignment (optional, requires -p)
         │
         ├──► 3_Total_Clusters_Trimmed.fasta
         └──► 4_Polished_Clusters_Trimmed.fasta
              │
              ▼
[Stage C]  SNP/InDel variant calling HTML report (optional, requires -r and -V)

[Stage D]  Best-per-locus sequence extraction (see NanoPolite_alpha_1.0.sh)
              │
              └──► 5_Best_Per_Locus.fasta
```

The pipeline is designed around a key principle: **locus assignment is primer-based, not reference-based**. Stage B uses the primer CSV you already have to assign each polished consensus to the correct locus, so no per-locus reference FASTA is needed for basic genotyping workflows.

---

## Pipeline Stages

### Read QC
Per-sample FASTQ files are quality-filtered using **NanoFilt**. Reads below the minimum Phred quality score or length threshold are discarded. Both plain and gzip-compressed FASTQ files are supported transparently.

### Read Clustering
Quality-filtered reads are clustered using **vsearch** (`--cluster_fast`) at a configurable sequence identity threshold (default 85%), with both strand orientations considered. Each cluster represents reads from a single amplicon locus. Clusters below the minimum depth threshold are retained as unpolished representatives; clusters meeting the threshold proceed to polishing.

### Consensus Polishing
For each cluster meeting the depth threshold, the longest read is used as a draft consensus. Reads are aligned to the draft using **minimap2** (ONT preset), and the draft is polished with **Racon** for one or more rounds. An experimental **Medaka** polisher is also supported via `--polisher medaka`. Polished sequences are labelled with sample name, cluster ID, and supporting read count in the FASTA header.

### Stage A — Reference-Based Orientation (optional)
Polished consensus sequences are aligned to a user-supplied reference FASTA using **BLASTn**. Sequences aligning to the minus strand are reverse-complemented, standardising all output to the forward orientation. Requires `-r`.

### Stage B — Primer Trimming and Locus Assignment (optional)
For each primer pair in the panel CSV, **seqkit amplicon** extracts and trims the matching amplicon from the consensus pool. Extracted sequences are labelled with the primer pair ID as a prefix, enabling locus assignment without a per-locus reference. A reverse-complement fallback recovers sequences that escaped Stage A orientation correction. Requires `-p`.

### Stage C — Variant Calling (optional)
Given a mapping table linking reference names to sample identifiers, NanoPolite calls SNPs and InDels by BLASTn alignment of polished trimmed sequences against reference sequences and parses the BTOP alignment string. Results are written as an interactive HTML report. Requires `-r` and `-V`.

---

## Output Files

| File | Description |
|---|---|
| `1_Total_Clusters.fasta` | All clusters (polished + unpolished) across all samples |
| `2_Polished_Clusters.fasta` | Polished clusters only |
| `3_Total_Clusters_Trimmed.fasta` | All clusters after primer trimming and locus assignment (Stage B) |
| `4_Polished_Clusters_Trimmed.fasta` | Polished clusters after primer trimming (Stage B) |
| `5_Best_Per_Locus.fasta` | One best-supported consensus per locus per sample (Stage D, alpha v1.0) |
| `nanopolite_stats_*.tsv` | Per-sample statistics: raw reads, clean reads, clusters, polished/unpolished counts |
| `Nanopolite_Report_*.html` | Summary HTML report with run statistics table |
| `Nanopolite_Variant_Report_*.html` | SNP/InDel variant HTML report (Stage C, if `-V` supplied) |
| `nanopolite_run_*.log` | Full timestamped run log including dependency versions and per-sample timing |

FASTA headers follow the format:
```
>SampleName_cluster_N_Reads=X_Polished
>SampleName_cluster_N_Reads=X_Unpolished
```
After Stage B, locus-assigned headers follow:
```
>LocusID_SampleName_cluster_N_Reads=X_Polished
```

---

## Dependencies

NanoPolite requires the following tools to be available in your `$PATH`:

| Tool | Version tested | Required for |
|---|---|---|
| [NanoFilt](https://github.com/wdecoster/nanofilt) | v2.8.0 | QC filtering |
| [vsearch](https://github.com/torognes/vsearch) | v2.30.4 | Read clustering |
| [minimap2](https://github.com/lh3/minimap2) | v2.30 | Draft alignment |
| [Racon](https://github.com/lbcb-sci/racon) | v1.5.0 | Consensus polishing (default) |
| [seqkit](https://bioinf.shenwei.me/seqkit/) | v2.13.0 | Sequence manipulation and amplicon trimming |
| [BLASTn](https://blast.ncbi.nlm.nih.gov/) | v2.9.0+ | Stage A orientation and Stage C variant calling |
| python3 | ≥3.7 | Reporting and Stage C variant calling |
| pandas | any recent | HTML report generation |
| [Medaka](https://github.com/nanoporetech/medaka) | optional | Alternative polisher (`--polisher medaka`) |

---

## Installation

NanoPolite is a single Bash script with no compilation required.

```bash
# Clone the repository
git clone https://github.com/www139516/Nanopolite.git
cd Nanopolite

# Make the script executable
chmod +x Nanopolite_alpha_1.0.sh

# Optional: add to PATH
export PATH=$PATH:$(pwd)
```

Install dependencies using conda (recommended):

```bash
conda create -n nanopolite -c bioconda -c conda-forge \
    nanofilt vsearch minimap2 racon seqkit blast python=3.10 pandas
conda activate nanopolite
```

---

## Usage

```
bash Nanopolite_alpha_1.0.sh [options]
```

At minimum, supply one or more input FASTQ files via `-i` or a directory via `-d`:

```bash
# Process a single sample
bash Nanopolite_alpha_1.0.sh -i sample01.fastq.gz -o results/

# Process all FASTQ files in a directory
bash Nanopolite_alpha_1.0.sh -d /path/to/fastq_dir/ -o results/

# Full workflow with orientation correction and primer trimming
bash Nanopolite_alpha_1.0.sh \
    -d /path/to/fastq_dir/ \
    -r reference.fasta \
    -p primers.csv \
    -o results/ \
    -T 16
```

---

## Parameter Reference

### Input / Output

| Flag | Description | Default |
|---|---|---|
| `-i`, `--input` | Input FASTQ file(s). Accepts `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz`, and wildcards | — |
| `-d`, `--dir` | Directory containing FASTQ files (compressed or plain) | — |
| `-o`, `--out` | Output directory path | `Nanopolite_results_TIMESTAMP` |

### QC Parameters

| Flag | Description | Default |
|---|---|---|
| `-q`, `--quality` | Minimum mean Phred quality score | `10` |
| `-m`, `--min-len` | Minimum read length (bp) | `20` |
| `--headcrop` | Trim N bases from read start | `0` |
| `--tailcrop` | Trim N bases from read end | `0` |
| `--chimera-check` | Enable vsearch de-novo chimera filtering after QC | off |

### Clustering Parameters

| Flag | Description | Default |
|---|---|---|
| `--id` | vsearch clustering identity threshold | `0.9` |
| `--min-cluster` | Minimum reads per cluster to trigger polishing | `3` |

### Polishing Parameters

| Flag | Description | Default |
|---|---|---|
| `--polisher` | Polisher to use: `racon` or `medaka` | `racon` |
| `--rounds` | Number of Racon polishing rounds | `1` |
| `--medaka-model` | Medaka model name (if `--polisher medaka`) | `r1041_e82_400bps_sup_v5.0.0` |

### Orientation and Trimming Parameters

| Flag | Description | Default |
|---|---|---|
| `-r`, `--ref-seq` | Reference FASTA for strand orientation (Stage A) | — |
| `-p`, `--primers` | Primer CSV for amplicon trimming and locus assignment (Stage B) | — |
| `--amp-min` | Minimum retained amplicon length (bp) after trimming | `500` |
| `--amp-max` | Maximum retained amplicon length (bp) after trimming | `1200` |
| `--mismatch` | Allowed mismatches for primer binding | `5` |
| `--search-range` | seqkit amplicon extraction region: `1:-1` (include primers) or `2:-2` (exclude primers) | `1:-1` |

### Variant Calling Parameters

| Flag | Description | Default |
|---|---|---|
| `-V`, `--table` | Mapping table for Stage C variant calling (requires `-r`) | — |

### Performance

| Flag | Description | Default |
|---|---|---|
| `-T`, `--threads` | CPU threads | `8` |

---

## Input File Formats

### Primer CSV (`-p`)

A comma-separated file with a header row. Each data row provides one primer pair:

```
ID,Forward,Reverse
M0143,ATGCATGCATGCATGC,GCTAGCTAGCTAGCTA
M0378,TTGGCCAATGCTAGCT,AGCTTAGCTAGCTAAG
```

- `ID`: locus identifier used as prefix in output FASTA headers
- `Forward`: forward primer sequence (5′→3′)
- `Reverse`: reverse primer sequence (5′→3′)

### Variant Mapping Table (`-V`)

A space- or tab-separated file (no header). Each row maps one reference sequence name to one or more partial sample name strings:

```
Locus01    sample01,sample02,sample03
Locus02    sample01,sample02
```

---

## Example Commands

### Minimal run (QC + clustering + polishing only)
```bash
bash Nanopolite_alpha_1.0.sh \
    -d ./fastq/ \
    -o results_basic/
```

### Standard multiplex amplicon genotyping workflow
```bash
bash Nanopolite_alpha_1.0.sh \
    -d ./fastq/ \
    -r ZmB73-V5.fa \
    -p ONTprimers.csv \
    --amp-min 500 \
    --amp-max 1200 \
    --mismatch 5 \
    -q 10 \
    -m 100 \
    --id 0.85 \
    --min-cluster 3 \
    --rounds 1 \
    -T 32 \
    -o results/
```
This is the command used in the associated Plant Methods manuscript for processing 142 maize inbred lines across 11 primer pairs.

### With chimera filtering
```bash
bash Nanopolite_alpha_1.0.sh \
    -d ./fastq/ \
    -r reference.fasta \
    -p primers.csv \
    --chimera-check \
    -T 16 \
    -o results_chimera_filtered/
```

### With Medaka polishing (experimental)
```bash
bash Nanopolite_alpha_1.0.sh \
    -d ./fastq/ \
    --polisher medaka \
    --medaka-model r1041_e82_400bps_sup_v5.0.0 \
    -T 8 \
    -o results_medaka/
```

### Haplotype or MSA analysis (exclude primer sequences from output)
```bash
bash Nanopolite_alpha_1.0.sh \
    -d ./fastq/ \
    -p primers.csv \
    --search-range 2:-2 \
    -o results_haplotype/
```

---

## Citation

If you use NanoPolite in your research, please cite this websit.

---

## License

NanoPolite is released under the [MIT License](LICENSE).

Copyright © 2026 Yuancong Wang, Institute of Germplasm Resources and Biotechnology, Jiangsu Academy of Agricultural Sciences (JAAS).
