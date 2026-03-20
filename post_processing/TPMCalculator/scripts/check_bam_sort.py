#!/usr/bin/env python3
"""
check_bam_sort.py

The script scans each sample directory, identifies the STAR-aligned BAM file,
determines its sort order using `samtools view -H`, and logs the results to
a TSV file for quick verification of alignment output consistency.

Functions:
    - get_sort_status(bam_path)
    - main()

Usage:
    python3 check_bam_sort.py
"""

import os
import glob
import subprocess
import csv

def get_sort_status(bam_path):
    """
    Returns samtools sort order: "coordinate", "queryname", "unsorted",
    or "unknown".

    Args:
        - bam_path (string): Path to bam file.

    Returns:
        - (string): samtools sort order.
    """
    try:

        # Run samtools command
        header = subprocess.run(
            ["samtools", "view", "-H", bam_path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        # Return sort order
        if header.returncode != 0:
            return "unknown"

        for line in header.stdout.splitlines():
            if line.startswith("@HD"):
                if "SO:coordinate" in line:
                    return "coordinate"
                elif "SO:queryname" in line:
                    return "queryname"
                elif "SO:unsorted" in line:
                    return "unsorted"
                else:
                    return "unsorted"

        return "unsorted"
    except Exception:
        return "unknown"


def main():
    """
    Finds bam files within a provided folder and logs the
    samtools sort order.

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
    """

    #################
    # Define paths
    ##################
    base_dir = "</path/to/base/dir>"
    per_sample_dir = os.path.join(base_dir, "PerSample")

    log_file_name = "bam_sort_status_" + os.path.basename(base_dir) + ".tsv"
    log_path = os.path.join(os.getcwd(), log_file_name)

    # Check if PerSample directory exists
    if not os.path.exists(per_sample_dir):
        print(f"Error: PerSample directory not found at {per_sample_dir}")
        return

    # Find all sample directories
    sample_dirs = [d for d in os.listdir(per_sample_dir) 
                   if os.path.isdir(os.path.join(per_sample_dir, d))]
    print(f"\nFound {len(sample_dirs)} sample directories")

    # Open TSV
    with open(log_path, "w", newline="") as log_file:
        writer = csv.writer(log_file, delimiter="\t")
        writer.writerow(["sample_dir", "bam_file", "sort_status"])

        # Loop over samples
        for sample_dir in sample_dirs:
            sample_path = os.path.join(per_sample_dir, sample_dir)
            results_star_htseq_path = os.path.join(sample_path, "results_STAR_htseq")

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

            # Check sort status and log result
            sort_status = get_sort_status(bam_file)
            writer.writerow([sample_dir, bam_file, sort_status])

    print(f"\nLogged BAM sort status to {log_path}\n")


if __name__ == "__main__":
    main()