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


DESCRIPTION = """index_star_salmon.py - Create genome index using STAR"""


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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('genomedir', help='genome directory')
    parser.add_argument('--genome_fasta', help='genome FASTA file')
    parser.add_argument('--runThreadN', type=int, default=32)
    parser.add_argument('--genomeChrBinNbits', type=int, default=16)
    parser.add_argument('--genomeSAindexNbases', type=int, default=12)
    parser.add_argument("--sjdbGTFfeatureExon")
    parser.add_argument("--sjdbGTFtagExonParentTranscript")
    parser.add_argument("--sjdbGTFtagExonParentGene")

    args = parser.parse_args()
    if args.genome_fasta is not None and os.path.exists(args.genome_fasta):
        genome_fasta = args.genome_fasta
    else:
        genome_fasta = glob.glob('%s/*.fasta' % (args.genomedir))[0]
    create_genome_index(args.genomedir, genome_fasta, args)
