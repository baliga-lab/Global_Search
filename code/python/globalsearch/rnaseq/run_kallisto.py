#!/usr/bin/env python3

"""
RNASeq Analysis pipeline using Kallisto
"""
import glob, sys, os, re, subprocess
from .find_files import find_fastq_files
from .trim_galore import trim_galore, collect_trimmed_data, create_result_dirs
from globalsearch.rnaseq.find_files import rnaseq_data_folder_list
import argparse
import json

# data and results directories
# These were the original examples
#RUN_DIR = "/proj/omics4tb2/Global_Search"
#DATA_DIR = "%s/Pilot_Pass/X204SC21081158-Z01-F003/raw_data" % RUN_DIR
#GENOME_DIR = "%s/reference_genomes/acerv_smic-reefGenomics" % RUN_DIR
#TRANSCRIPTOME_FILE = "%s/Acerv_Smic-ReefGenomics_merged.fasta" % GENOME_DIR

############# Functions ##############
####################### Run Kalisto ###############################
def run_kallisto(index_path, results_dir, pair_files):
    print('\033[33mRunning kallisto! \033[0m')
    # flatten the pair_files list into an input file list
    input_files = []
    for pair_file in pair_files:
        input_files.extend(pair_file)

    command = ['kallisto',
              'quant', '-i', index_path]
    command.extend(input_files)
    command.extend(['-o', results_dir,
                    '-b', '100', '--bias', '-t', '4', '--rf-stranded'])
    kallisto_cmd = ' '.join(command)
    print('Kallisto run command: "%s"' % kallisto_cmd)
    compl_proc = subprocess.run(command, check=True, capture_output=False)

 ####################### Create Kallisto index ###############################
def kallisto_index(index_path, transcriptome_path):
    print('\033[33mRunning kallisto index! \033[0m')
    if os.path.exists(index_path):
        print("index path '%s' exists, skipping the indexing" % index_path)
        return

    command = ['kallisto',
               'index',
               '-i', index_path,
               transcriptome_path]
    kallistoindex_cmd = ' '.join(command)
    print('kallisto index command: %s' % kallistoindex_cmd)
    compl_proc = subprocess.run(command, check=True, capture_output=False)

####################### Running the Pipeline ###############################
def run_pipeline(data_folder, results_folder, genome_dir, transcriptome_file, args):
    folder_count = 1

    # Loop through each data folder
    folder_name = data_folder.split('/')[-1]
    print('\033[33mProcessing Folder: %s\033[0m' %(folder_name))

    # Get the list of first file names in paired end sequences
    #first_pair_files = glob.glob('%s/*_1.fq*' %(data_folder))
    #print(first_pair_files)
    patterns = args.fastq_patterns.split(',')
    pair_files = find_fastq_files(data_folder, args.fastq_patterns.split(','))
    print("PATTERNS: ", patterns)
    print("PAIR FILES: ", pair_files)

    # Program specific results directories
    organism = os.path.basename(genome_dir)
    data_trimmed_dir = os.path.join(results_folder, folder_name, "trimmed")
    fastqc_dir = os.path.join(results_folder, folder_name, "fastqc_results")
    results_dir = os.path.join(results_folder, folder_name, organism)
    htseq_dir = os.path.join(results_folder, "htseq_counts")

    # Run create directories function to create directory structure
    create_result_dirs(data_trimmed_dir,fastqc_dir,results_dir, htseq_dir)
    print("PAIR_FILES: ", pair_files, flush=True)

    index_path = os.path.join(genome_dir, "%s_kallistoindex" % organism)
    kallisto_index(index_path, transcriptome_file)

    # Loop through each file and create filenames
    file_count = 1
    is_gzip = True
    is_paired_end = True
    for pair_file in pair_files:
        first_pair_file, second_pair_file = pair_file
        if second_pair_file is None:
            is_paired_end = False

        #first_file_name_full = first_pair_file.split('/')[-1]
        #second_pair_file = first_pair_file.replace('_1.fq', '_2.fq')
        #second_file_name_full = second_pair_file.split('/')[-1]
        #file_ext = first_pair_file.split('.')[-1]
        fastq_fname = os.path.basename(first_pair_file)
        is_gzip = fastq_fname.endswith("gz")
        if is_gzip:
            file_ext = '.'.join(fastq_fname.split('.')[-2:])
        else:
            file_ext = fastq_fname.split('.')[-1]

        #print ('\033[32m Processing File: %s of %s (%s)\033[0m' %(file_count, len(first_pair_files), first_file_name_full ))
        print('\033[32m Processing file set: %s of %s (first is "%s")\033[0m' % (file_count, len(pair_files),
                                                                                 first_pair_file),
              flush=True)

        #first_file_name = re.split('.fq|.fq.gz',first_file_name_full)[0]
        #second_file_name = re.split('.fq|.fq.gz',second_file_name_full)[0]
        #print('first_file_name:%s, second_file_name:%s' %(first_file_name,second_file_name))

        # Collect Sample attributes
        #exp_name = folder_name
        #print("exp_name: %s" %(exp_name))
        #lane = first_file_name.split("_")[-1]
        #print("Lane: %s" %(lane))
        #sample_id = re.split('.fq|.fq.gz', first_file_name)[0]
        #print("sample_id: %s" %(sample_id))
        sample_id = fastq_fname.replace(file_ext, "")
        print("sample_id: %s" % sample_id, flush=True)

        # 01. Run TrimGalore
        # WW: Activate/deactivate for debugging
        trim_galore(first_pair_file, second_pair_file, folder_name, sample_id,
                    data_trimmed_dir, fastqc_dir)
        file_count += 1

        # Run folder level salmon analysis
        first_pair_group, second_pair_group = collect_trimmed_data(data_trimmed_dir, is_gzip, is_paired_end)
        run_kallisto(index_path, results_dir, [first_pair_group.split(' '), second_pair_group.split(' ')])

        folder_count += 1

    return data_trimmed_dir,fastqc_dir,results_dir


class KallistoArgs:
    def __init__(self, config):
        self.fastq_patterns = ','.join(config['fastq_patterns'])

def run_config(config):
    args = KallistoArgs(config)
    includes = rnaseq_data_folder_list(config)
    for folder in includes:
        data_folder = os.path.join(config["input_dir"], folder)
        results_folder = os.path.join(config["output_dir"], folder)
        genome_dir = config["genome_dir"]
        transcriptome_file = config['genome_fasta']
        run_pipeline(data_folder, results_folder, genome_dir, transcriptome_file,
                     args)

DESCRIPTION = """run_kallisto.py - run Kallisto pipeline"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('genomedir', help='genome directory')
    parser.add_argument('dataroot', help="parent of input directory")
    parser.add_argument('indir', help="input directory without the input dir part")
    parser.add_argument('transcriptome_file', help="path to transcriptome_file")
    parser.add_argument('outdir', help='output directory')
    parser.add_argument('--fastq_patterns', help="FASTQ file patterns", default="*_{{pairnum}}.fq.*")
    # if config is set, this will replace everything
    parser.add_argument('--config', help="config file")
    args = parser.parse_args()

    if not args.config is None:
        with open(args.config) as infile:
            config = json.load(infile)
            run_config(config)
    else:
        print("Processing directory %s" % args.indir)
        data_folder = os.path.join(args.dataroot, args.indir)
        data_trimmed_dir, fastqc_dir, results_dir = run_pipeline(data_folder, args.outdir, args.genomedir,
                                                                 args.transcriptome_file, args)
