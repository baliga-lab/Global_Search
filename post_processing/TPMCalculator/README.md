
# mRNA Quantification using TPMCalculator

## Introduction

This pipeline quantifies mRNA abundance directly from BAM file alignments using **TPMCalculator** and generates visualizations to assess the results. It is designed to process multiple samples and combine results into comprehensive matrices for downstream analysis. 

The pipeline is organized into four scripts, each responsible for a different stage of processing and analysis. These scripts require access to a **SLURM** job manager for submitting and running the workflow.

## Installation
This repository can be cloned locally by running the following `git` command:
```bash
git clone https://github.com/baliga-lab/Global_Search.git
```
Please note that Git is required to run the above command. For instructions on downloading Git, please see [the Git Guide](https://github.com/git-guides/install-git).

<a id="environment"></a>
### Environment

#### Conda
This application is built on top of multiple Python packages with specific version requirements. We recommend using `conda` to create an isolated Python environment with all necessary packages. If you do not have it installed, follow the [Miniconda Instillation Instructions](https://www.anaconda.com/docs/getting-started/miniconda/install). The list of necessary packages can be found at in the [`environment.yml`](./environment.yml) file.

To create the specified `gs-coral-env` Conda environment, run the following command:
```bash
conda env create -f environment.yml
```

Once the Conda environment is created, it can be activated by:
```bash
conda activate gs-coral-env
```
After coding inside the environment, it can be deactivated with the command:
```bash
conda deactivate
```

**Please note! This environment file does not contain the packages required to run the whole Global_Search pipeline. The packages included are only enough to run the scripts in the TPMCalculator section**. Installation instructions for the Gloabl_Search pipleine can be found in the main [README](../../README.md)

## TPMCalculator

TPMCalculator requires two inputs:  
- A **GTF annotation file** describing gene models.  
- **BAM alignment files** for each sample.  

It produces four output files per sample containing TPM values and raw read counts for:  
- **Genes**  
- **Transcripts**  
- **Exons**  
- **Introns**

## GFF to GTF

While there are several ways to convert a GFF file to a GTF file, a custom script was used to facilitate this conversion. The [`gff3_to_gtf.py`](./scripts/gff3_to_gtf.py) script...
- Reads a GFF3 file
- Converts `gene`, `mRNA/transcript`, `exon`, and `CDS` features into GTF format
- Constructs valid `gene_id` and `transcript_id` fields
- Writes a new GTF file

Run the following command to convert a GFF3 file to a GTF file:
```bash
python3 gff3_to_gtf.py <input_gff3> <output_gtf>
```

## Pipeline Overview

The scripts should be run in the following order:

0. **Generate the GTF file (if applicable)**
   - See instructions above. 

1. **check_bam_sorted.py**  
   - **Purpose:** Verify that all BAM files are correctly sorted and log their sort order.  
   - **Inputs:**  
     - `base_dir`: Path to the folder containing BAM files. BAM files can be nested within sample directories.  
   - **Outputs:**  
     - TSV file logging the sort order of each BAM file (`bam_sort_status_<sample>.tsv`).  

2. **generate_tpmcalculator_results.py**  
   - **Purpose:** Generate TPMCalculator results for each sample and create SLURM job scripts to run the analyses.  
   - **Inputs:**  
     - `base_dir`: Path to the folder containing BAM files. Bam files may be nested.
     - `gtf_file`: Path to the reference GTF annotation file.
     - `tpmcalculator_dir`: Desired output directory for TPMCalculator results.  
   - **Outputs:**  
     - TPMCalculator output files for each sample (genes, transcripts, exons, introns).  
     - SLURM job scripts for each sample (`*_tpmcalculator_job.sh`).  

3. **combine_tpmcalculator_results.py**  
   - **Purpose:** Combine TPMCalculator results from all samples into unified TPM and reads matrices.  
   - **Inputs:**  
     - `base_dir`: Folder containing all TPMCalculator outputs for all samples. Output files may be nested.
     - `output_dir`: Desired output directory for combined matrices.  
   - **Outputs:**  
     - TPM matrix (`STAR_TPMcalculator_<project>_TPMs_Matrix_Merged.csv`).  
     - Reads matrix (`STAR_TPMcalculator_<project>_Reads_Matrix_Merged.csv`).  

4. **analyze_organism_totals.py**  
   - **Purpose:** Summarize total TPMs and reads by organism (e.g., Host, Symbiont species) and generate visualizations to verify results.  
   - **Inputs:**  
     - `base_dir`: Folder containing combined TPMCalculator results. This is the output directory.  
     - `tpm_file`: Path to the merged TPM matrix.  
     - `reads_file`: Path to the merged reads matrix. 
   - **Outputs:**  
     - Summary CSV files of TPM and reads totals by organism (`TPM_totals_by_organism.csv`, `Reads_totals_by_organism.csv`).  
     - Barplots and summary figures (`organism_totals_barplots.png/pdf`, `organism_summary_stats.png/pdf`).  

## Recommended Directory Structure

For clarity and organization, the pipeline outputs are recommended to be structured as follows:
```
├── sample_organism_name/
│   ├── combined_results/
│   │   ├── organism_analysis/
│   │   │   ├── <summary visualization files>
│   │   │   ├── ...
│   │   ├── STAR_TPMCalculator_<sample_organism_name>_Reads_Matrix_Merged.csv
│   │   ├── STAR_TPMCalculator_<sample_organism_name>_TPMs_Matrix_Merged.csv
│   ├── tpmcalculator/
│   │   ├── sample_1/
│   │   │   ├── <TPMCalculator output files>
│   │   │   ├── ...
│   │   ├── sample_2/
│   │   │   ├── <TPMCalculator output files>
│   │   │   ├── ...
│   │   ├── sample_.../
│   │   │   ├── <TPMCalculator output files>
│   │   │   ├── ...
│   ├── reference/
│   │   ├── <sample_organism_name>.gff3 or gff
│   │   ├── <sample_organism_name>.gtf
```
