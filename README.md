# igvf-data-collect

Here are some scripts to collect analysis-sets and per-sample information from the IGVF portal.

## Overview

This repository contains tools to extract sequencing data information from the IGVF portal and download the associated FASTQ files.

## Scripts

### 1. analysis_to_sample.py

Extracts data from the IGVF portal for a given analysis set ID and generates two TSV files:
- An analysis sets file mapping IDs to measurement/auxiliary sets (tsv)
- A per-sample file containing sequencing details (paths, checksums, modalities, etc.) (tsv)

#### Usage

```bash
python analysis_to_sample.py -i ANALYSIS_SET_ID --access-key ACCESS_KEY --secret-key SECRET_KEY -per-sample-output per_sample_file.tsv -analysis-set-output analysis_sets.tsv
```

### 2. download_fastq.py

Downloads all FASTQ files listed in the per-sample file generated from the previous step.

#### Usage

```bash
python download_fastq.py --sample per-sample_file.tsv --access-key ACCESS_KEY --secret-key SECRET_KEY
```

This will create a directory `fastq_files` and download all FASTQ files into this directory.

