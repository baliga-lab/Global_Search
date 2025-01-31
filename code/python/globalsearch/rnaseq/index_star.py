#!/usr/bin/env python3

"""
Module to encapsulate STAR indexing.

It takes the genome directory and one or more FASTA files
and passes them to STAR to generate a genome index within
the genome directory.
If the index exists, it will skip the generation, to avoid
wasting time as this is a very costly step.
"""
import argparse
import os
import subprocess
import glob
import json


DESCRIPTION = """index_star_salmon.py - Create genome index using STAR"""

class STARIndexArgs:
    def __init__(self, runThreadN, genomeChrBinNbits, genomeSAindexNbases,
                 sjdbGTFfeatureExon, sjdbGTFtagExonParentTranscript,
                 sjdbGTFtagExonParentGene):
        # mandatory
        self.runThreadN = runThreadN
        self.genomeChrBinNbits = genomeChrBinNbits
        self.genomeSAindexNbases = genomeSAindexNbases
        # optional
        self.sjdbGTFfeatureExon = None
        self.sjdbGTFtagExonParentTranscript = None
        self.sjdbGTFtagExonParentGene = None


####################### Create STAR index ###############################
### This should be specific for the organism
### Use the equation file maybe another script to create references
def create_genome_index(genome_dir, genome_fasta, args):
    index_command = ['STAR', '--runMode', 'genomeGenerate',
                     '--runThreadN', str(args.runThreadN),
                     '--genomeDir', genome_dir,
                     '--genomeFastaFiles', genome_fasta,
                     '--genomeChrBinNbits', str(args.genomeChrBinNbits),
                     '--genomeSAindexNbases', str(args.genomeSAindexNbases)]
    # optional commands
    if args.sjdbGTFfeatureExon is not None:
        index_command += ["sjdbGTFfeatureExon", args.sjdbGTFfeatureExon]
    if args.sjdbGTFtagExonParentTranscript is not None:
        index_command += ["sjdbGTFtagExonParentTranscript", args.sjdbGTFtagExonParentTranscript]
    if (args.sjdbGTFtagExonParentTranscript is not None and
        args.sjdbGTFtagExonParentGene is not None):
        index_command += ["sjdbGTFtagExonParentGene", args.sjdbGTFtagExonParentGene]

    index_cmd = ' '.join(index_command)
    print("RUNNING STAR in index MODE: '%s'" % index_cmd, flush=True)

    print ("\033[34m %s Indexing genome... \033[0m", flush=True)
    if os.path.exists('%s/SAindex' % (genome_dir)):
        print ('Genome indexes exist. Not creating!', flush=True)
    else:
        print('Creating genome indexes', flush=True)
        compl_proc = subprocess.run(index_command, check=True, capture_output=False, cwd=genome_dir)
        print('finished indexing with STAR', flush=True)


def run_config(configfile):
    with open(configfile) as infile:
        config = json.load(infile)

    genome_fasta = config['genome_fasta']
    genome = os.path.basename(os.path.normpath(config['genome_dir']))
    print("running index on: ", genome, " genome fasta: ", genome_fasta)
    star_index_options = config['star_index_options']
    try:
        sjdbGTFfeatureExon = star_index_options['sjdbGTFfeatureExon']
    except KeyError:
        sjdbGTFfeatureExon = None

    try:
        sjdbGTFtagExonParentTranscript = star_index_options['sjdbGTFtagExonParentTranscript']
    except KeyError:
        sjdbGTFtagExonParentTranscript = None

    try:
        sjdbGTFtagExonParentGene = star_index_options['sjdbGTFtagExonParentGene']
    except KeyError:
        sjdbGTFtagExonParentGene = None

    index_args = STARIndexArgs(star_index_options['runThreadN'],
                               star_index_options['genomeChrBinNbits'],
                               star_index_options['genomeSAindexNbases'],
                               sjdbGTFfeatureExon, sjdbGTFtagExonParentTranscript,
                               sjdbGTFtagExonParentGene)
    genome_dir = config['genome_dir']
    create_genome_index(genome_dir, genome_fasta, index_args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('genomedir', help='genome directory')
    parser.add_argument('--config', help='config file, overrides everything else')
    parser.add_argument('--genome_fasta', help='genome FASTA file')
    parser.add_argument('--runThreadN', type=int, default=32)
    parser.add_argument('--genomeChrBinNbits', type=int, default=16)
    parser.add_argument('--genomeSAindexNbases', type=int, default=12)
    parser.add_argument("--sjdbGTFfeatureExon")
    parser.add_argument("--sjdbGTFtagExonParentTranscript")
    parser.add_argument("--sjdbGTFtagExonParentGene")

    args = parser.parse_args()
    if args.config is not None:
        # override: using config file
        run_config(args.config)
    else:
        if args.genome_fasta is not None and os.path.exists(args.genome_fasta):
            genome_fasta = args.genome_fasta
        else:
            genome_fasta = glob.glob('%s/*.fasta' % (args.genomedir))[0]
        create_genome_index(args.genomedir, genome_fasta, args)
