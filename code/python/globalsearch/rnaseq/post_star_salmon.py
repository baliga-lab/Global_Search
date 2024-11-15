#!/usr/bin/env python3
from .variantcall import run_bcftools

"""
Post-processing step for STAR Salmon analysis.

This includes

  1. extracting TPM and read information from salmon quant
  2. running MultiQC
"""
import argparse
import json
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector
import subprocess
import os

DESCRIPTION = """post_star_salmon.py - Post-run step for STAR Salmon"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('configfile', help='configuration file')
    args = parser.parse_args()
    with open(args.configfile) as infile:
        config = json.load(infile)

    # Determine post run directory. If not specified, make it sub directory of
    # output directory
    output_dir = config['output_dir']
    try:
        postrun_outdir = config['postrun_output_dir']
    except:
        postrun_outdir = os.path.join(output_dir, "Post_Run_Results")
    if not os.path.exists(postrun_outdir):
        os.makedirs(postrun_outdir)

    genome_dir = config['genome_dir']
    organisms = config['organisms']
    #org1, org2 = os.path.basename(genome_dir).split('_')
    print('\033[33mExtracting salmon quant files...\033[0m')
    global_search = importr("GlobalSearch")
    global_search.extract_salmon_quants(StrVector(organisms), output_dir, postrun_outdir)
    # now run MultiQC
    print('\033[33mRunning MultiQC...\033[0m')
    multiqc_outdir = os.path.join(postrun_outdir, 'MultiQC')
    if not os.path.exists(multiqc_outdir):
        os.makedirs(multiqc_outdir)
    command = ['multiqc', '--outdir', multiqc_outdir, output_dir]
    compl_proc = subprocess.run(command, check=True, capture_output=False, cwd=output_dir)

    # run bcftools
    run_bcftools(config, postrun_outdir)

