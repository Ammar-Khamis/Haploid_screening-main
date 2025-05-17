# Haploid Screening Workflow



# Haploid Screening Pipeline

# README.md


This file provides guidance for working with code in this repository.

## Repository Purpose

This repository contains scripts for a haploid genetic screening workflow using Oxford Nanopore Technology (ONT) sequencing data. The workflow analyzes eHAP1 cells with PiggyBac (PB) transposon insertions across different experimental conditions, particularly comparing DMSO (control) versus Cordycepin treatment.

## Core Workflow

The main workflow is defined in `enhanced_complete_Haploid_Screening_workflow.sh` and consists of these key steps:

1. Genome preparation and indexing
2. Data processing for each barcode/condition:
   - Adapter trimming (`Improved_trim_fastq.py`)
   - TTAA site filtering (`Filter_TTAA.py`)
   - Mapping to reference genome (minimap2)
   - Insertion site extraction (`Extract_Insertion_Sites.py`)
3. Gene annotation and frequency analysis
4. Visualization and statistics generation

## Key Scripts and Usage

### Data Processing

```bash
# Process raw ONT data through the complete pipeline
./enhanced_complete_Haploid_Screening_workflow.sh
```

```bash
# Extract insertion sites from a BAM file
samtools view mapped.bam | python Extract_Insertion_Sites.py --output insertion_sites.bed
```

```bash
# Filter FASTQ reads for TTAA insertion sites
python Filter_TTAA.py --input reads.fastq --output filtered.fastq
```

### Visualization

Run the visualization scripts on processed matrices:

```bash
# Generate visualizations for matrix data (v1)
./Run_visualization_script.sh

# Generate enhanced visualizations with Excel output (v2)
./Run_visualization_script_v2.sh
```

### Pathway Analysis

Replot IPA (Ingenuity Pathway Analysis) data:

```bash
# Run the latest version of the IPA visualization tool
python replot_ipa_bars_v6.py
```

## Data Organization

- Primary experiment data is processed into gene and chromosome matrices
- Processed data is stored in matrices:
  - `normalized_matrix_full.txt`: Complete matrix with all genes
  - `normalized_matrix_clean.txt`: Filtered matrix
  - `normalized_matrix_split.txt`: Matrix split by gene

- Visualization outputs include:
  - Heatmaps showing insertion patterns
  - Chromosome-level visualizations
  - Condition comparison plots
  - Pathway visualization plots

## Analysis Parameters

- Treatment conditions: DMSO (control), Cordycepin (150µM, 200µM)
- Vector systems: hyPB (hyperactive PiggyBac), suPB (super PiggyBac)
- TTAA sites are used for filtering PiggyBac insertions
- Insertion site merging window: 5bp (configurable in Extract_Insertion_Sites.py)
- Minimum mapping quality: 10 (configurable in Extract_Insertion_Sites.py)

## Environment Setup

This workflow requires several Conda environments:

```bash
# For running genomic tools (samtools, bedtools, minimap2)
conda activate biotools

# For running LAST and perl scripts
conda activate last_env
```

## Common Analysis Tasks

```bash
# To analyze specific barcodes for a new dataset
# Modify the ONT_BARCODE parameter in the script
ONT_BARCODE="barcode07" ./enhanced_complete_Haploid_Screening_workflow.sh

# To regenerate visualizations for existing data
# Modify the ANALYSIS_DIR, GENE_MATRIX, and CHROM_MATRIX paths in Run_visualization_script_v2.sh
./Run_visualization_script_v2.sh

# To replot pathway analysis with different thresholds
# Edit the GENE_XLSX and IPA_FILE variables at the top of replot_ipa_bars_v6.py
python replot_ipa_bars_v6.py
```

## Dependencies

- Python packages: pandas, numpy, matplotlib, seaborn, xlsxwriter, scipy, statsmodels
- Bioinformatics tools: samtools, bedtools, minimap2, LAST

- For pathway analysis: Ingenuity Pathway Analysis (IPA) export files

- For pathway analysis: Ingenuity Pathway Analysis (IPA) exports

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.



## Installation

Install the required Python packages:

```bash
pip install -r requirements.txt
```

Install the external tools (samtools, bedtools, minimap2, LAST). The easiest
approach is via Bioconda:

```bash
conda install -c bioconda samtools bedtools minimap2 last
```



