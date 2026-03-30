#!/usr/bin/env python3
"""
generate_tpmcalculator_jobs.py

This module generates and submits SLURM job scripts for running TPMCalculator on
multiple RNA-seq samples. It scans a PerSample directory containing sample-specific
subdirectories, identifies the appropriate BAM file for each sample,
and generates a SLURM job script that runs TPMCalculator using a
provided GTF annotation file.

Functions:
    - create_slurm_job_script(sample_id, bam_file, gtf_file, tmp_output_dir)
    - main()

Usage:
    python3 generate_tpmcalculator_slurm_jobs.py
"""

import os
import glob
import subprocess

def create_slurm_job_script(sample_id, bam_file, gtf_file, tmp_output_dir):
    """
    Create a SLURM job script for TPMCalculator analysis.
    
    Args:
        - sample_id (string): Sample identifier
        - bam_file (string): Path to the BAM file
        - gtf_file (string): Path to the GTF file
        - tmp_output_dir (string): Directory to save TPMCalculator output

    Returns:
        - job_script_path (string): Path to SLURM job bash script
    """

    job_script_content = f"""#!/bin/bash
#SBATCH --job-name=tpmcalc_{sample_id}
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output={tmp_output_dir}/{sample_id}_tpmcalculator_%j.out
#SBATCH --error={tmp_output_dir}/{sample_id}_tpmcalculator_%j.err

# TPMCalculator is available via conda installation

# Define variables
BAM="{bam_file}"
GTF="{gtf_file}"
OUTPUT_DIR="{tmp_output_dir}"

# Run TPMCalculator
echo "Starting TPMCalculator analysis for {sample_id}..."
mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR
TPMCalculator -g $GTF -b $BAM -o $OUTPUT_DIR/{sample_id}_tpmcalculator.txt

echo "TPMCalculator analysis completed for {sample_id}"
"""
    
    # Create the bash file containing the SLURM job
    job_script_path = os.path.join(tmp_output_dir, f"{sample_id}_tpmcalculator_job.sh")
    
    with open(job_script_path, 'w') as f:
        f.write(job_script_content)
    
    # Make the script executable
    os.chmod(job_script_path, 0o755)
    
    return job_script_path

def main():
    """
    This function scans a project results directory, identifies sample
    BAM files produced by the STAR pipeline, generates a SLURM job script
    for each sample, and submits those jobs to the cluster.

    Expected Input Directory Structure
    ----------------------------------
    The script expects the following directory structure:

        base_dir/
        ├── PerSample/
        │   ├── SAMPLE_1/
        │   │   └── results_STAR_htseq/
        │   │       └── *_dedup_NoSingletonCollated.out.bam
        │   │
        │   ├── SAMPLE_2/
        │   │   └── results_STAR_htseq/
        │   │       └── *_dedup_NoSingletonCollated.out.bam
        │   │
        │   └── ...

    Where:
        - Each directory inside `PerSample/` represents one RNA-seq sample.
        - Each sample directory must contain a `results_STAR_htseq/` folder.
        - The BAM file used for TPMCalculator must match the pattern:
          `*dedup_NoSingletonCollated.out.bam`.

    Required Paths
    --------------
    req_paths is a dictionary containing required input and output paths:

        base_dir
            Root directory of the RNA-seq project results. This directory
            must contain the `PerSample/` directory with all sample folders.

        gtf_file
            Path to the reference GTF annotation file used by TPMCalculator
            to assign reads to genes and transcripts.

        output_dir
            Directory where TPMCalculator outputs and generated SLURM job
            scripts will be written. A subdirectory for each sample will
            be created automatically if it does not already exist.
    """

    #################
    # Define paths
    ##################
    req_paths ={
        "base_dir": "</path/to/base/dir>",
        "gtf_file": "</path/to/gtf/file>",
        "output_dir": "</path/to/output/dir>"
    }
    per_sample_dir = os.path.join(req_paths["base_dir"], "PerSample")
    
    # Check if GTF file exists
    if not os.path.exists(req_paths["gtf_file"]):
        print(f"Error: GTF file not found at {req_paths["gtf_file"]}")
        return
    
    # Check if PerSample directory exists
    if not os.path.exists(per_sample_dir):
        print(f"Error: PerSample directory not found at {per_sample_dir}")
        return
    
    # Find all sample directories
    sample_dirs = [d for d in os.listdir(per_sample_dir) 
                   if os.path.isdir(os.path.join(per_sample_dir, d))]

    print(f"Found {len(sample_dirs)} sample directories")

    job_scripts_created = []
    
    for sample_dir in sample_dirs:
        sample_path = os.path.join(per_sample_dir, sample_dir)
        results_star_htseq_path = os.path.join(sample_path, "results_STAR_htseq")
        
        # Check if results_STAR_htseq directory exists
        if not os.path.exists(results_star_htseq_path):
            print(f"Warning: results_STAR_htseq directory not found for {sample_dir}")
            continue
        
        # Find BAM file
        bam_pattern = os.path.join(results_star_htseq_path, "*dedup_NoSingletonCollated.out.bam")
        bam_files = glob.glob(bam_pattern)

        if not bam_files:
            print(f"Warning: No BAM file found for {sample_dir}")
            continue
        elif len(bam_files) > 1:
            print(f"Warning: Multiple BAM files found for {sample_dir}, using the first one")
        
        bam_file = bam_files[0]
        
        # Create tpmcalculator directory
        tpmcalculator_dir = req_paths["output_dir"] + sample_dir
        os.makedirs(tpmcalculator_dir, exist_ok=True)
        
        # Extract sample ID from directory name
        sample_id = sample_dir
        
        # Create SLURM job script
        job_script_path = create_slurm_job_script(sample_id, bam_file, req_paths["gtf_file"], tpmcalculator_dir)
        job_scripts_created.append(job_script_path)
        
        print(f"Created job script for {sample_id}: {job_script_path}")
    
    print(f"\nTotal job scripts created: {len(job_scripts_created)}")
    
    # Submit all jobs automatically
    print("\nSubmitting jobs to SLURM cluster...")
    submitted_jobs = []

    for job_script in job_scripts_created:
        try:
            result = subprocess.run(['sbatch', job_script], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            if result.returncode == 0:
                job_id = result.stdout.strip().split()[-1]
                submitted_jobs.append(job_id)
                print(f"Submitted job {job_id} for {job_script}")
            else:
                print(f"Failed to submit {job_script}: {result.stderr.strip()}")
        except Exception as e:
            print(f"Error submitting {job_script}: {e}")

    print(f"\nSuccessfully submitted {len(submitted_jobs)} jobs out of {len(job_scripts_created)} total jobs.")


if __name__ == "__main__":
    main()