#!/usr/bin/env python3
"""
combine_tpmcalculator_results.py

Combine TPMCalculator output files from multiple RNA-seq samples into merged matrices.
This script scans the TPMCalculator output directories and aligns TPM and Reads
values across all samples by gene. Two merged CSV are generated.

Functions:
    - extract_sample_name(filepath)
    - find_tpmcalculator_outputs(base_dir)
    - combine_tpm_matrices(output_files, output_dir, req_paths)
    - main()

Usage:
    python3 combine_tpmcalculator_results.py
"""

import os
import pandas as pd
import glob
import re
from pathlib import Path

def extract_sample_name(filepath):
    """
    Extract sample name from filepath by removing '_star_10_0.3_0.66_0_dedup' part.

    Args:
        - filepath (str): Path to the TPMCalculator output file

    Returns:
        - sample_name (str): Sample name
    """
    filename = os.path.basename(filepath)
    # Remove the file extension and TPMCalculator suffix
    base_name = filename.replace('_tpmcalculator.txt_genes.out', '')
    # Remove the STAR alignment suffix
    sample_name = re.sub(r'_star_10_0\.3_0\.66_0_dedup.*', '', base_name)
    return sample_name

def find_tpmcalculator_outputs(base_dir):
    """
    Find all TPMCalculator output files with various possible patterns.

    Args:
        - base_dir (str): Base directory containing PerSample folders

    Returns:
        - list: List of file paths
    """
    # per_sample_dir = os.path.join(base_dir, "PerSample")
    per_sample_dir = os.path.join(base_dir, "tpmcalculator_output")

    # Try multiple patterns for TPMCalculator output files
    patterns = [
        os.path.join(per_sample_dir, "*", "*Collated.out_genes"),
        os.path.join(per_sample_dir, "*", "*_genes.out"),
        os.path.join(per_sample_dir, "*", "*NoSingletonCollated.out_genes.out"),
        os.path.join(per_sample_dir, "*", "*dedup*.out_genes.out"),
    ]

    print(patterns)

    output_files = []
    for pattern in patterns:
        files = glob.glob(pattern)
        if files:
            output_files = files
            print(f"Found files with pattern: {pattern}")
            break

    return sorted(output_files)

def combine_tpm_matrices(output_files, output_dir, req_paths):
    """
    Combine TPMCalculator results into TPM and Reads matrices.

    Args:
        - output_files (list): List of TPMCalculator output file paths
        - output_dir (str): Directory to save the combined matrices
        - req_paths (dict): Dictionary containing names of TPM and Reads CSVs
    """
    tpm_data = {}
    reads_data = {}
    gene_ids = None

    print(f"Processing {len(output_files)} TPMCalculator output files...")

    for i, filepath in enumerate(output_files):
        try:
            # Extract sample name
            sample_name = extract_sample_name(filepath)
            print(f"Processing {i+1}/{len(output_files)}: {sample_name}")

            # Check file size first
            file_size = os.path.getsize(filepath)
            if file_size < 300:  # Less than 300 bytes is likely empty (just headers)
                print(f"Warning: File {filepath} is too small ({file_size} bytes), likely empty")
                continue

            # Read the TPMCalculator output file
            df = pd.read_csv(filepath, sep='\t')

            # Identify duplicated Gene_Id
            dup_mask = df['Gene_Id'].duplicated(keep=False)
            if dup_mask.any():
                skipped_file = os.path.join(output_dir, "Skipped_Duplicate_Gene_Ids.csv")

                # Save skipped rows
                df[dup_mask].assign(sample=sample_name).to_csv(
                    skipped_file,
                    sep=",",
                    index=False,
                    mode="a",
                    header=not os.path.exists(skipped_file)
                )
                df = df[~dup_mask]

            # Check if dataframe is empty
            if df.empty:
                print(f"Warning: File {filepath} contains no data")
                continue

            # Check if required columns exist
            required_cols = ['Gene_Id', 'TPM', 'Reads']
            missing_cols = [col for col in required_cols if col not in df.columns]

            if missing_cols:
                print(f"Warning: Missing columns {missing_cols} in {filepath}")
                print(f"Available columns: {list(df.columns)}")
                continue

            # Check if there are any genes
            if len(df) == 0:
                print(f"Warning: No genes found in {filepath}")
                continue

            # Set gene IDs from first file
            if gene_ids is None:
                gene_ids = df['Gene_Id'].tolist()
                print(f"Using gene list from {sample_name} with {len(gene_ids)} genes")

            # Check if gene lists match
            if gene_ids != df['Gene_Id'].tolist():
                print(f"Warning: Gene list mismatch in {sample_name}")
                # Try to align with existing gene_ids
                df_aligned = df.set_index('Gene_Id').reindex(gene_ids)
                tpm_data[sample_name] = df_aligned['TPM'].fillna(0).tolist()
                reads_data[sample_name] = df_aligned['Reads'].fillna(0).tolist()
            else:
                # Extract TPM and Reads columns
                tpm_data[sample_name] = df['TPM'].tolist()
                reads_data[sample_name] = df['Reads'].tolist()

        except Exception as e:
            print(f"Error processing {filepath}: {e}")
            continue

    if not tpm_data:
        print("No valid data found. Exiting.")
        return

    print(f"Successfully processed {len(tpm_data)} samples")
    print(f"Found {len(gene_ids)} genes")

    # Create TPM matrix
    tpm_df = pd.DataFrame(tpm_data)
    tpm_df.insert(0, 'Gene_Id', gene_ids)

    # Create Reads matrix
    reads_df = pd.DataFrame(reads_data)
    reads_df.insert(0, 'Gene_Id', gene_ids)

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Save matrices
    tpm_output = os.path.join(output_dir, req_paths["tpm_csv"])
    reads_output = os.path.join(output_dir, req_paths["reads_csv"])

    tpm_df.to_csv(tpm_output, index=False)
    reads_df.to_csv(reads_output, index=False)

    print(f"\nMatrices saved:")
    print(f"TPM matrix: {tpm_output}")
    print(f"Reads matrix: {reads_output}")
    print(f"Matrix dimensions: {tpm_df.shape}")
    print(f"Samples: {list(tpm_data.keys())}")

def main():
    """
    This function combines the TPMCalculator results across all samples,
    by gene, into merged TPM and Reads matrices. 

    Expected Input Directory Structure
    ----------------------------------
    The script expects the following structure:

        base_dir/
        └── tpmcalculator_output/
            ├── SAMPLE_1/
            │   └── SAMPLE_1...out_genes.out
            ├── SAMPLE_2/
            │   └── SAMPLE_2...out_genes.out
            └── ...

    Required Paths
    --------------------
    req_paths contains user-defined paths and output file names:

        base_dir
            Root directory containing 'tpmcalculator_output/'. This is the output directory
            from the previous script.

        tpm_csv
            File name for the merged TPM matrix.

        reads_csv
            File name for the merged Reads matrix.
    """

    #################
    # Define paths
    ##################
    req_paths = {
        "base_dir": "</path/to/base/dir>",
        "tpm_csv": "STAR_TPMcalculator_TPMs_Matrix_Merged_rerun.csv",
        "reads_csv": "STAR_TPMcalculator_Reads_Matrix_Merged_rerun.csv"
    }
    output_dir = req_paths["base_dir"] + "/combined_results"

    # Find TPMCalculator output files
    print("Searching for TPMCalculator output files...")
    output_files = find_tpmcalculator_outputs(req_paths["base_dir"])

    if not output_files:
        print("No TPMCalculator output files found!")
        print("Expected pattern: tpmcalculator/*/*Collated.out_sorted_genes")
        return

    print(f"Found {len(output_files)} TPMCalculator output files")

    # Show first few files for verification
    print("\nFirst few files found:")
    for f in output_files[:5]:
        sample_name = extract_sample_name(f)
        print(f"  {sample_name}: {f}")

    if len(output_files) > 5:
        print(f"  ... and {len(output_files) - 5} more")

    # Combine matrices
    combine_tpm_matrices(output_files, output_dir, req_paths)
    

if __name__ == "__main__":
    main()
