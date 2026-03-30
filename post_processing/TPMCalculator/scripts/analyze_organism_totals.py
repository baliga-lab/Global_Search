#!/usr/bin/env python3
"""
analyze_organism_totals.py

This script takes merged TPM and Reads matrices produced by the previous 
TPMCalculator combination step and verifes that expression values were
assinged correctly to each organism for the multi-organism dataset. 

It also creates summary dataframes and barplots to verify processing worked correctly,
and generates the following outputs:
    - TPM_totals_by_organism.csv
    - Reads_totals_by_organism.csv
    - organism_totals_barplots.png/.pdf
    - organism_summary_stats.png/.pdf

Functions:
    - analyze_organism_totals(tpm_file, reads_file, output_dir, organisms)
    - create_barplots(tpm_totals, reads_totals, output_dir, organisms)
    - create_summary_stats_plot(tpm_totals, reads_totals, output_dir, organisms)
    - main()

Usage:
    python3 analyze_organism_totals.py
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def analyze_organism_totals(tpm_file, reads_file, output_dir, organisms):
    """
    Analyze total TPMs and reads by organism prefix.

    Args:
        - tpm_file (str): Path to TPM matrix CSV file
        - reads_file (str): Path to reads matrix CSV file
        - output_dir (str): Directory to save outputs
        - organisms (dict): Dictionary containing the names and desired colors used
            for graph labeling

    Returns:
        - tpm_totals (pandas.DataFrame): DataFrame with TPM totals by organism
        - reads_totals (pandas.DataFrame): DataFrame with reads totals by organism
    """
    print("Loading TPM and reads matrices...")

    # Load matrices
    tpm_df = pd.read_csv(tpm_file)
    reads_df = pd.read_csv(reads_file)

    print(f"TPM matrix shape: {tpm_df.shape}")
    print(f"Reads matrix shape: {reads_df.shape}")

    # Get sample columns (all except Gene_Id)
    sample_columns = [col for col in tpm_df.columns if col != 'Gene_Id']
    print(f"Found {len(sample_columns)} samples")

    # Define organism prefixes
    organism_prefixes = organisms["names"]

    # Initialize result dataframes
    tpm_totals = pd.DataFrame(index=sample_columns, columns=organism_prefixes)
    reads_totals = pd.DataFrame(index=sample_columns, columns=organism_prefixes)

    print("Calculating totals by organism...")

    # Calculate totals for each organism and sample
    for sample in sample_columns:
        for prefix in organism_prefixes:
            # Filter genes by prefix
            gene_mask = tpm_df['Gene_Id'].str.startswith(prefix)

            # Calculate totals
            tpm_total = tpm_df.loc[gene_mask, sample].sum()
            reads_total = reads_df.loc[gene_mask, sample].sum()

            tpm_totals.loc[sample, prefix] = tpm_total
            reads_totals.loc[sample, prefix] = reads_total

    # Convert to numeric
    tpm_totals = tpm_totals.astype(float)
    reads_totals = reads_totals.astype(float)

    # Add total columns
    tpm_totals['Total'] = tpm_totals.sum(axis=1)
    reads_totals['Total'] = reads_totals.sum(axis=1)

    print("Summary statistics:")
    print("\nTPM totals:")
    print(tpm_totals.describe())
    print("\nReads totals:")
    print(reads_totals.describe())

    # Save summary dataframes
    tpm_output = os.path.join(output_dir, "TPM_totals_by_organism.csv")
    reads_output = os.path.join(output_dir, "Reads_totals_by_organism.csv")

    tpm_totals.to_csv(tpm_output)
    reads_totals.to_csv(reads_output)

    print(f"\nSaved summary files:")
    print(f"TPM totals: {tpm_output}")
    print(f"Reads totals: {reads_output}")

    # Create visualizations
    create_barplots(tpm_totals, reads_totals, output_dir, organisms)

    return tpm_totals, reads_totals

def create_barplots(tpm_totals, reads_totals, output_dir, organisms):
    """
    Create barplots for TPM and reads totals by organism.

    Args:
        - tpm_totals (pandas.DataFrame): DataFrame with TPM totals by organism
        - reads_totals (pandas.DataFrame): DataFrame with reads totals by organism
        - output_dir (str): Directory to save plots
        - organisms (dict): Dictionary containing the names and desired colors used
            for graph labeling

    Returns:
        - N/A
    """
    print("Creating visualizations...")

    # Set up the plotting style
    plt.style.use('default')
    colors = organisms["colors"]

    # Remove Total column for plotting
    organism_cols = organisms["names"]

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('TPM and Reads Distribution by Organism', fontsize=16, fontweight='bold')

    # Plot 1: TPM totals by organism (stacked bar)
    ax1 = axes[0, 0]
    tpm_plot_data = tpm_totals[organism_cols]
    tpm_plot_data.plot(kind='bar', stacked=True, ax=ax1, width=0.8, color=colors)
    ax1.set_title('TPM Totals by Organism (Stacked)', fontweight='bold')
    ax1.set_xlabel('Sample Index')
    ax1.set_ylabel('TPM')
    ax1.legend(title='Organism', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.set_xticks([])  # Remove x-axis labels for readability

    # Plot 2: Reads totals by organism (stacked bar)
    ax2 = axes[0, 1]
    reads_plot_data = reads_totals[organism_cols]
    reads_plot_data.plot(kind='bar', stacked=True, ax=ax2, width=0.8, color=colors)
    ax2.set_title('Reads Totals by Organism (Stacked)', fontweight='bold')
    ax2.set_xlabel('Sample Index')
    ax2.set_ylabel('Reads')
    ax2.legend(title='Organism', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.set_xticks([])  # Remove x-axis labels for readability

    # Plot 3: TPM percentage by organism
    ax3 = axes[1, 0]
    tpm_percentages = tpm_plot_data.div(tpm_plot_data.sum(axis=1), axis=0) * 100
    tpm_percentages.plot(kind='bar', stacked=True, ax=ax3, width=0.8, color=colors)
    ax3.set_title('TPM Distribution by Organism (%)', fontweight='bold')
    ax3.set_xlabel('Sample Index')
    ax3.set_ylabel('Percentage (%)')
    ax3.legend(title='Organism', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax3.set_xticks([])  # Remove x-axis labels for readability

    # Plot 4: Reads percentage by organism
    ax4 = axes[1, 1]
    reads_percentages = reads_plot_data.div(reads_plot_data.sum(axis=1), axis=0) * 100
    reads_percentages.plot(kind='bar', stacked=True, ax=ax4, width=0.8, color=colors)
    ax4.set_title('Reads Distribution by Organism (%)', fontweight='bold')
    ax4.set_xlabel('Sample Index')
    ax4.set_ylabel('Percentage (%)')
    ax4.legend(title='Organism', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax4.set_xticks([])  # Remove x-axis labels for readability

    # Adjust layout
    plt.tight_layout()

    # Save plot
    plot_output = os.path.join(output_dir, "organism_totals_barplots.png")
    plt.savefig(plot_output, dpi=300, bbox_inches='tight')
    plt.savefig(plot_output.replace('.png', '.pdf'), bbox_inches='tight')

    print(f"Saved plots: {plot_output}")
    plt.close()

    # Create summary statistics plot
    create_summary_stats_plot(tpm_totals, reads_totals, output_dir, organisms)

def create_summary_stats_plot(tpm_totals, reads_totals, output_dir, organisms):
    """
    Create summary statistics visualization.

    Args:
        - tpm_totals (pandas.DataFrame): DataFrame with TPM totals by organism
        - reads_totals (pandas.DataFrame): DataFrame with reads totals by organism
        - output_dir (str): Directory to save plots
        - organisms (dict): Dictionary containing the names and desired colors used
            for graph labeling

    Returns:
        - N/A
    """
    organism_cols = organisms["names"]
    color_cols = organisms["colors"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle('Summary Statistics by Organism', fontsize=16, fontweight='bold')

    # TPM summary
    ax1 = axes[0]
    tpm_summary = tpm_totals[organism_cols].mean()
    bars1 = ax1.bar(range(len(organism_cols)), tpm_summary.values,
                   color=color_cols)
    ax1.set_title('Mean TPM by Organism', fontweight='bold')
    ax1.set_xlabel('Organism')
    ax1.set_ylabel('Mean TPM')
    ax1.set_xticks(range(len(organism_cols)))
    ax1.set_xticklabels(organism_cols, rotation=45)

    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars1, tpm_summary.values)):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + val*0.01,
                f'{val:.0f}', ha='center', va='bottom')

    # Reads summary
    ax2 = axes[1]
    reads_summary = reads_totals[organism_cols].mean()
    bars2 = ax2.bar(range(len(organism_cols)), reads_summary.values,
                   color=color_cols)
    ax2.set_title('Mean Reads by Organism', fontweight='bold')
    ax2.set_xlabel('Organism')
    ax2.set_ylabel('Mean Reads')
    ax2.set_xticks(range(len(organism_cols)))
    ax2.set_xticklabels(organism_cols, rotation=45)

    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars2, reads_summary.values)):
        ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + val*0.01,
                f'{val:.0f}', ha='center', va='bottom')

    plt.tight_layout()

    # Save plot
    summary_output = os.path.join(output_dir, "organism_summary_stats.png")
    plt.savefig(summary_output, dpi=300, bbox_inches='tight')
    plt.savefig(summary_output.replace('.png', '.pdf'), bbox_inches='tight')

    print(f"Saved summary plot: {summary_output}")
    plt.close()

def main():
    """
    Genenerates analysis and visualization written to the
    base_dir/organism_analysis/ directory.

    Expected Input Directory Structure
    ----------------------------------
    The script expects the following structure:

        base_dir/
        └── combined_results/
        │   ├── STAR_TPMcalculator_*_TPMs_Matrix_Merged.csv
        │   ├── STAR_TPMcalculator_*_Reads_Matrix_Merged.csv
        └── tpmcalculator_output/
            └── ...

    Organism Dictionary
    -------------------
    organisms contains two fields used throughout the analysis:

        names
            List of gene ID prefixes used to identify which organism each
            gene belongs to. For example, a gene ID like "Sym_Cgor_gene123"
            should be labeled as "Sym_Cgor_".

        colors
            List of colors corresponding to each organism prefix.

    Required Paths
    --------------------
    req_paths contains user-defined paths and output file names:

        combined_results_dir
            Directory containing the 'combined_results/' folder from the
            previous script in the pipeline.

        tpm_file
            File name for the merged TPM matrix.

        reads_file
            File name for the merged Reads matrix.
    """

    ####################
    # Define Organisms
    ####################
    organisms = {
        "names": ['Host_', 'Sym_Cgor_', 'Sym_Dtre_', 'Sym_Smic_'],
        "colors": ['#1f77b4', '#ff7f0e', '#2ca02c' ,'#d62728'] 
    }

    ####################
    # Define paths
    ####################
    req_paths = {
        "combined_results_dir": "</path/to/combined_results/dir>",
        "tpm_file": "STAR_TPMcalculator_TPMs_Matrix_Merged_rerun.csv",
        "reads_file": "STAR_TPMcalculator_Reads_Matrix_Merged_rerun.csv"
    }

    tpm_file = os.path.join(req_paths["combined_results_dir"], req_paths["tpm_file"])
    reads_file = os.path.join(req_paths["combined_results_dir"], req_paths["reads_file"])

    if os.path.exists(reads_file):
        print("Using original reads matrix (may have different gene set)")
    else:
        reads_file = tpm_file
        print("No reads matrix found, using TPM matrix for both analyses")

    # Check if input files exist
    if not os.path.exists(tpm_file):
        print(f"TPM matrix file not found: {tpm_file}")
        print("Please run combine_tpmcalculator_results.py first")
        return

    # Create output directory for analysis results
    output_dir = os.path.join(req_paths["combined_results_dir"], "organism_analysis")
    os.makedirs(output_dir, exist_ok=True)
    
    # Analyze organism totals
    tpm_totals, reads_totals = analyze_organism_totals(tpm_file, reads_file, output_dir, organisms)

    print("\nAnalysis completed successfully!")
    print(f"Results saved in: {output_dir}")

if __name__ == "__main__":
    main()