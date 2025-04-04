#!/usr/bin/env python3

#############################################################
##### RNASeq Analysis Pipeline with STAR                #####
##### Last update: 11/15/2022 Serdar Turkarslan         #####
##### Institute for Systems Biology                     #####
############################################################
import glob, sys, os, string, datetime, re, shutil
import argparse
import subprocess
import json

from .find_files import find_fastq_files, rnaseq_data_folder_list
from .trim_galore import trim_galore, collect_trimmed_data, create_result_dirs

DESCRIPTION = """run_STAR_SALMON.py - run STAR and Salmon"""


class HtseqArgs:
    """Class that fakes a command line argument object"""
    def __init__(self, config):
        htseq_options = config['htseq_options']
        self.htseqStranded = htseq_options['stranded']
        self.htseqFeatureType = htseq_options['feature_type']
        self.htseqID = htseq_options['id_attribute']
        self.htseqOrder = htseq_options['order']


class STARSalmonArgs:
    """Class that fakes a command line argument object"""
    def __init__(self, config):
        star_options = config['star_options']
        self.dedup = config['deduplicate_bam_files']
        self.use_htseq = config["rnaseq_algorithm"] == "star_htseq"
        self.fastq_patterns = ','.join(config['fastq_patterns'])
        self.salmon_genome_fasta = config['genome_fasta']
        self.runThreadN = config['star_index_options']['runThreadN']
        self.limitBAMsortRAM = 5784458574
        self.outFilterMatchNmin = star_options['outFilterMatchNmin']
        self.outFilterMismatchNmax = star_options['outFilterMismatchNmax']
        self.outFilterMismatchNoverLmax = star_options['outFilterMismatchNoverLmax']
        self.outFilterScoreMinOverLread = star_options['outFilterScoreMinOverLread']
        self.twopassMode = star_options['twopassMode']
        self.genome_gff = config['genome_gff']
        self.outSAMattributes = ','.join(config['star_options']['outSAMattributes'])
        try:
            self.sjdbGTFfeatureExon = star_options['sjdbGTFfeatureExon']
        except:
            self.sjdbGTFfeatureExon = None
        try:
            self.sjdbGTFtagExonParentGene = star_options['sjdbGTFtagExonParentGene']
        except:
            self.sjdbGTFtagExonParentGene = None
        try:
            self.sjdbGTFtagExonParentTranscript = star_options['sjdbGTFtagExonParentTranscript']
        except:
            self.sjdbGTFtagExonParentTranscript = None
        try:
            self.quantMode = star_options['quantMode']
        except:
            self.quantMode = None
        try:
            self.limitSjdbInsertNsj = star_options['limitSjdbInsertNsj']
        except:
            self.limitSjdbInsertNsj = None
        try:
            self.sjdbOverhang = star_options['sjdbOverhang']
        except:
            self.sjdbOverhang = None

        try:
            self.tmp = star_options['tmp']
        except:
            self.tmp = "/tmp"


        dedup_prefix = '_dedup' if config['deduplicate_bam_files'] else ''
        self.starPrefix = 'star_%s_%s_%s_%s%s' % (star_options['outFilterMismatchNmax'],
                                                  star_options['outFilterMismatchNoverLmax'],
                                                  star_options['outFilterScoreMinOverLread'],
                                                  star_options['outFilterMatchNmin'],
                                                  dedup_prefix)
        self.salmonPrefix = 'salmon_%s_%s_%s_%s%s' % (star_options['outFilterMismatchNmax'],
                                                      star_options['outFilterMismatchNoverLmax'],
                                                      star_options['outFilterScoreMinOverLread'],
                                                      star_options['outFilterMatchNmin'],
                                                      dedup_prefix)

####################### Run STAR #####################################
### We need to add Read GRoup info
### --outSAMattrRGline ID:${i%_TF_R1_val_1.fq.gz}
### https://github.com/BarshisLab/danslabnotebook/blob/main/CBASSAS_GenotypeScreening.md

def run_star(first_pair_group, second_pair_group, results_dir, folder_name,
             genome_dir, is_gzip, args):
    print('\033[33mRunning STAR! \033[0m', flush=True)
    outfile_prefix = '%s/%s_%s_' % (results_dir, folder_name, args.starPrefix)
    star_options = ["--runThreadN", str(args.runThreadN),
                    "--outFilterType", "Normal",
                    "--outSAMstrandField", "intronMotif",
                    "--outFilterIntronMotifs", "RemoveNoncanonical",
                    "--outSAMtype", "BAM", "Unsorted",
                    "--limitBAMsortRAM", str(args.limitBAMsortRAM)]
    if is_gzip:
        star_options.extend(["--readFilesCommand", "zcat"])
    star_options.extend(["--outReadsUnmapped", "Fastx",
                         "--outFilterMismatchNmax", str(args.outFilterMismatchNmax),
                         "--outFilterMismatchNoverLmax", str(args.outFilterMismatchNoverLmax),
                         "--outFilterScoreMinOverLread", str(args.outFilterScoreMinOverLread),
                         "--outFilterMatchNmin", str(args.outFilterMatchNmin)])

    genome_load = "LoadAndKeep"  # This is the default, for efficiency
    if args.twopassMode:
        star_options.extend(["--twopassMode", "Basic"])
        genome_load = "NoSharedMemory"  # two-pass has to run without shared memory

    command = ["STAR", "--genomeDir", genome_dir]
    command += star_options
    if args.outSAMattributes != "Standard" and len(args.outSAMattributes) > 0:
        out_sam_attrs = args.outSAMattributes.split()
        command.append('--outSAMattributes')
        command += out_sam_attrs

    # Handling for GFF files
    if not args.genome_gff is None and os.path.exists(args.genome_gff):
        genome_load = "NoSharedMemory"  # can't use GFF with a shared genome memory
        gff_args = [
            '--sjdbGTFfile', args.genome_gff,
            '--sjdbGTFtagExonParentTranscript', args.sjdbGTFtagExonParentTranscript,
            '--limitSjdbInsertNsj', str(args.limitSjdbInsertNsj)
        ]
        if args.sjdbOverhang is not None:
            gff_args += ['--sjdbOverhang', str(args.sjdbOverhang)]
        command += gff_args

    command += [ "--readFilesIn", first_pair_group,
                 second_pair_group,
                 "--outFileNamePrefix", outfile_prefix]
    command += ["--genomeLoad", genome_load]

    # add more optional arguments
    if args.sjdbGTFfeatureExon is not None:
        command += ["--sjdbGTFfeatureExon", args.sjdbGTFfeatureExon]
    if args.sjdbGTFtagExonParentGene is not None:
        command += ["--sjdbGTFtagExonParentGene", args.sjdbGTFtagExonParentGene]
    if args.quantMode is not None:
        command += ["--quantMode"] + args.quantMode
    #if args.tmp is not None:
    #    command += ["--outTmpDir", args.tmp]

    cmd = ' '.join(command)
    compl_proc = subprocess.run(command, check=True, capture_output=False, cwd=results_dir)

####################### Samtools sorting and indexing ##########
#
def run_samtools_sort_and_index(results_dir):
    bam_files = glob.glob(os.path.join(results_dir, "*.bam"))
    if len(bam_files) == 0:
        print("ERROR: could not sort and index - bam file not found")
        return
    sorted_bam_path = None
    for f in bam_files:
        if f.endswith("Sorted.out.bam"):
            sorted_bam_path = f

    filename = os.path.basename(bam_files[0])
    if sorted_bam_path is None:
        print("Using samtools to sort STAR BAM result")
        stem = bam_files[0].replace(".bam", "").replace(".out", "")
        sorted_bam_path = os.path.join(results_dir, "%s_Sorted.out.bam" % stem)
        command = ['samtools', 'sort', bam_files[0], '-o',
                   sorted_bam_path]
        compl_proc = subprocess.run(command, check=True,
                                    capture_output=False, cwd=results_dir)
    if not os.path.exists(sorted_bam_path + ".bai"):
        print("Using samtools to index sorted STAR BAM result")
        command = ["samtools", "index", sorted_bam_path]
        compl_proc = subprocess.run(command, check=True,
                                    capture_output=False, cwd=results_dir)


####################### Deduplication (not in _old) ###############################
def dedup(results_dir, folder_name, args):
    print('\033[33mRunning Deduplication! \033[0m', flush=True)
    outfile_prefix = '%s/%s_%s_' %(results_dir, folder_name, args.starPrefix)

    aligned_bam = '%sAligned.out.bam' % (outfile_prefix)
    fixmate_bam = '%sFixmate.out.bam' % (outfile_prefix)
    ordered_bam = '%sOrdered.out.bam' % (outfile_prefix)
    markdup_bam = '%sMarkedDup.out.bam' % (outfile_prefix)
    markdupSTAR_bam = '%sProcessed.out.bam' % (outfile_prefix)
    nosingleton_bam = '%sNoSingleton.out.bam' % (outfile_prefix)
    nosingletonCollated_bam = '%sNoSingletonCollated.out.bam' % (outfile_prefix)

    # STAR mark duplicates
    star_markdup_command = ['STAR', '--runThreadN', str(args.runThreadN),
                            '--runMode',
                            'inputAlignmentsFromBAM',
                            '--bamRemoveDuplicatesType', 'UniqueIdenticalNotMulti',
                            '--inputBAMfile', aligned_bam,
                            '--outFileNamePrefix', outfile_prefix]
    star_markdup_cmd = ' '.join(star_markdup_command)

    # remove duplicates marked by STAR with 0x400
    rmsingletonsSTAR_command = ['samtools', 'view', '-@', '8',
                                '-b', '-F', '0x400', markdupSTAR_bam,
                                '>', nosingleton_bam]
    rmsingletonsSTAR_cmd = ' '.join(rmsingletonsSTAR_command)

    # Collate reads by name
    collatereadsSTAR_command = ['samtools', 'sort', '-o',
                                nosingletonCollated_bam,
                                '-n', '-@', '8', nosingleton_bam]
    collatereadsSTAR_cmd = ' '.join(collatereadsSTAR_command)

    ## STAR based BAM duplicate removal
    # Mark duplicates with STAR
    print('STAR mark duplicates run command:%s' % star_markdup_cmd, flush=True)
    compl_proc = subprocess.run(star_markdup_command, check=True, capture_output=False, cwd=results_dir)

    # Remove marked duplicates withh samtools
    print('Samtools  STAR Dedup Remove run command:%s' % rmsingletonsSTAR_cmd, flush=True)
    compl_proc = subprocess.run(rmsingletonsSTAR_cmd, shell=True, check=True, capture_output=False, cwd=results_dir)

    # Remove marked duplicates withh samtools
    print('Samtools  Collate reads by read name run command:%s' % collatereadsSTAR_cmd, flush=True)
    compl_proc = subprocess.run(collatereadsSTAR_cmd, shell=True, check=True, capture_output=False, cwd=results_dir)


####################### Run Salmon Count ###############################

def get_final_bam_name(results_dir, folder_name, args):
    """determine the name of the input BAM file that comes out of STAR
    This is to preserve it upon cleanup
    """
    outfile_prefix = '%s/%s_%s_' %(results_dir, folder_name, args.starPrefix)
    print("OUTFILE PREFIX: ", outfile_prefix, flush=True)
    # check if we are performing deduplication
    if args.dedup:
        salmon_input = '%sNoSingletonCollated.out.bam' % (outfile_prefix)
    else:
        salmon_input = '%sAligned.out.bam' % (outfile_prefix)

        # Use BAM file aligned to transcriptome for salmon input if it exists
        salmon_transcriptome_input = "%sAligned.toTranscriptome.out.bam" % outfile_prefix
        if os.path.exists(salmon_transcriptome_input):
            salmon_input = salmon_transcriptome_input
    return salmon_input

def get_salmon_result_dir(results_dir, args):
    """create name of salmon result directory so we know where to cleanup"""
    return '%s/%s_salmon_quant' % (results_dir, args.salmonPrefix)


# WW: Check the names of the input files they will be different from _out
def run_salmon_quant(final_bam, results_dir, folder_name, genome_fasta, args):
    print('\033[33mRunning salmon-quant! \033[0m', flush=True)
    command = ['salmon', 'quant', '-t', genome_fasta,
               '-l', 'A',  '-a', final_bam,
               '-o', get_salmon_result_dir(results_dir, args)]
    cmd = ' '.join(command)
    print("Salmon quant command: '%s'" % cmd, flush=True)
    # run as a joined string
    compl_proc = subprocess.run(cmd, check=True, capture_output=False, cwd=results_dir, shell=True)


####################### Run HTSEq Count ###############################

def run_htseq_count(final_bam, htseq_resultdir,
                    folder_name, genome_gff, args):
    if not os.path.exists(htseq_resultdir):
        os.makedirs(htseq_resultdir)
    # sort and index the final_bam file for htseq, it is required
    # 1. check if exists
    htseq_inputfile = os.path.basename(final_bam).replace(".out.bam",
                                                          ".sortedByCoord.out.bam")
    htseq_inputpath = os.path.join(os.path.dirname(final_bam), htseq_inputfile)
    if not os.path.exists(htseq_inputpath):
        # SORT
        cmd = ["samtools", "sort", final_bam,
               "-o", htseq_inputpath]
        command = " ".join(cmd)
        print("Running '%s'" % command)
        os.system(command)
        index_file = htseq_inputfile.replace("bam", "bai")
        index_path = os.path.join(os.path.dirname(final_bam), index_file)
        if not os.path.exists(index_path):
            cmd = ["samtools", "index", htseq_inputpath]
            command = " ".join(cmd)
            print("Running '%s'" % command)
            os.system(command)

    htseq_resultfile = os.path.join(htseq_resultdir,
                                    "%s_htseqcounts.txt" % folder_name)
    if not os.path.exists(htseq_resultfile):
        cmd = ["htseq-count",
               "-s", args.htseqStranded,
               "-t", args.htseqFeatureType,
               "-i", args.htseqID,
               "-r", args.htseqOrder,
               "--max-reads-in-buffer", "60000000",
               "-f", "bam",
               htseq_inputpath,
               genome_gff,
               ">",
               htseq_resultfile]
        command = " ".join(cmd)
        print("Running '%s'" % command)
        os.system(command)
    else:
        print("'%s' exists -> skipping" % htseq_resultfile)

##### Cleanup step

def cleanup_after_run(cleanup_config, data_trimmed_dir, results_dir,
                      folder_name, args):
    if cleanup_config["intermediate_bam_files"]:
        final_bam = get_final_bam_name(results_dir, folder_name,
                                       args)
        print("final_bam: ", final_bam)
        if os.path.exists(final_bam):
            print("FINAL_BAM EXISTS")
            final_bam_name = os.path.basename(final_bam)
            #for f in os.listdir(results_dir):
            #    if f.endswith(".bam") and f != final_bam_name:
            #        print("delete BAM file: '%s'" % f)
            #        os.remove(os.path.join(results_dir, folder_name, f))
        else:
            print("FINAL_BAM DOES NOT EXIST")

    if cleanup_config["trim_files"]:
        print("trimmed_dir: ", data_trimmed_dir)
        if os.path.exists(data_trimmed_dir):
            print("TRIMMED_DIR EXISTS")
            #shutil.rmtree(data_trimmed_dir)
        else:
            print("TRIMMED_DIR DOES NOT EXIST")

    if cleanup_config["salmon_log"]:
        salmon_dir = get_salmon_result_dir(results_dir, args)
        if os.path.exists(salmon_dir):
            print("SALMON_DIR EXISTS")
            salmon_log_dir = os.path.join(salmon_dir, "logs")
            if os.path.exists(salmon_log_dir):
                #for f in glob(os.path.join(salmon_log_dir, "*.log")):
                #    os.remove(f)
                print("SALMON LOG DIR EXISTS !!! -> REMOVE")
                #shutil.rmtree(salmon_log_dir)
            else:
                print("SALMON LOG DIR DOES NOT EXIST")
        else:
            print("SALMON_DIR DOES NOT EXIST")


####################### Running the Pipeline ###############################

def run_pipeline(data_folder, results_folder, genome_dir, genome_fasta, args,
                 config):

    # Loop through each data folder
    folder_name = data_folder.split('/')[-1]
    print('\033[33mProcessing Folder: %s\033[0m' % folder_name, flush=True)

    # Get the list of first file names in paired end sequences
    ## We need to make sure we capture fastq data files
    pair_files = find_fastq_files(data_folder, args.fastq_patterns.split(','))

    # Program specific results directories
    data_trimmed_dir = "%s/%s/trimmed" % (results_folder,folder_name)
    fastqc_dir = "%s/%s/fastqc_results" % (results_folder,folder_name)

    results_dir = "%s/%s/results_STAR_Salmon" %(results_folder, folder_name)
    htseq_dir = "%s/htseq_counts" % (results_dir)

    # Run create directories function to create directory structure
    create_result_dirs(data_trimmed_dir, fastqc_dir, results_dir, htseq_dir)

    print("PAIR_FILES: ", pair_files, flush=True)

    # Loop through each file and create filenames
    file_count = 1

    is_gzip = True
    is_paired_end = True

    for pair_file in pair_files:
        first_pair_file, second_pair_file = pair_file
        if second_pair_file is None:
            is_paired_end = False

        fastq_fname = os.path.basename(first_pair_file)
        is_gzip = fastq_fname.endswith("gz")
        if is_gzip:
            file_ext = '.'.join(fastq_fname.split('.')[-2:])
        else:
            file_ext = fastq_fname.split('.')[-1]

        print('\033[32m Processing file set: %s of %s (first is "%s")\033[0m' % (file_count, len(pair_files),
                                                                                 first_pair_file),
              flush=True)

        # Collect Sample attributes
        sample_id = fastq_fname.replace(file_ext, "")
        print("sample_id: %s" % sample_id, flush=True)

        # Run TrimGalore
        trim_galore(first_pair_file, second_pair_file, folder_name,sample_id, data_trimmed_dir, fastqc_dir)
        file_count += 1

    # Collect Trimmed data for input into STAR
    first_pair_group, second_pair_group = collect_trimmed_data(data_trimmed_dir, is_gzip, is_paired_end)

    # Run STAR
    run_star(first_pair_group, second_pair_group, results_dir, folder_name, genome_dir, is_gzip, args)

    # Run samtools, sorting and indexing
    run_samtools_sort_and_index(results_dir)

    # Run Deduplication
    if args.dedup:
        print('\033[33mRunning Deduplication: \033[0m', flush=True)
        dedup(results_dir,folder_name, args)

    final_bam = get_final_bam_name(results_dir, folder_name, args)

    if args.use_htseq:
        htseq_args = HtseqArgs(config)
        run_htseq_count(final_bam, htseq_dir, folder_name, args.genome_gff,
                        htseq_args)
    else:
        # Run Salmon Quant
        if args.salmon_genome_fasta is not None:
            genome_fasta = args.salmon_genome_fasta

        run_salmon_quant(final_bam, results_dir, folder_name, genome_fasta, args)

    try:
        cleanup_config = config["cleanup_after_run"]
        print("CLEANUP AFTER RUN...")
        cleanup_after_run(cleanup_config, data_trimmed_dir, results_dir,
                          folder_name, args)
    except:
        print("CLEANUP AFTER RUN - SKIPPED")
        pass

    return data_trimmed_dir, fastqc_dir, results_dir


def aws_s3_sync(result_dir, bucket_url):
    parent_dir = os.path.dirname(result_dir)
    compl_proc = subprocess.run("aws s3 sync %s %s" % (parent_dir, bucket_url),
                                shell=True, check=True, capture_output=False)
    compl_proc = subprocess.run("rm -rf %s" % result_dir,
                                shell=True, check=True, capture_output=False)

def run_config(configfile):
    """Run from config file"""
    with open(configfile) as infile:
        config = json.load(infile)

    starsalmon_args = STARSalmonArgs(config)
    data_folders = sorted(rnaseq_data_folder_list(config))
    for data_folder in data_folders:
        run_pipeline(os.path.join(config['input_dir'], data_folder),
                     config['output_dir'],
                     config['genome_dir'],
                     config['genome_fasta'],
                     starsalmon_args,
                     config)
        # post run action, currently always an S3 sync
        try:
            actions = config["nocluster_post_run"]
            result_dir = os.path.join(config['output_dir'], data_folder)
            for action in actions:
                action_type = action["action_type"]
                url = action["url"]
                aws_s3_sync(result_dir, url)
        except:
            # ignore
            pass


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('genomedir', help='genome directory')
    parser.add_argument('dataroot', help="parent of input directory")
    parser.add_argument('indir', help="input directory (R<somenumber>)")
    parser.add_argument('outdir', help='output directory')
    parser.add_argument("--use_htseq", help="htseq instead of salmon", action="store_true")
    parser.add_argument('--fastq_patterns', help="FASTQ file patterns", default="*_{{pairnum}}.fq.*")
    parser.add_argument('--genome_gff', help='genome GFF file')
    parser.add_argument('--genome_fasta', help='genome FASTA file')
    parser.add_argument('--dedup', action='store_true', help='should we deduplicate bam files (True or False)')
    parser.add_argument('--twopassMode', action='store_true', help='run STAR in two-pass mode')
    parser.add_argument('--starPrefix', help="STAR output file name prefix")
    parser.add_argument('--salmonPrefix', help="Salmon output folder name prefix")
    parser.add_argument('--outFilterMismatchNmax', nargs='?', const=10, type=int)
    parser.add_argument('--outFilterMismatchNoverLmax', nargs='?', const=0.3, type=float)
    parser.add_argument('--outFilterScoreMinOverLread', nargs='?', const=0.66, type=float)
    parser.add_argument('--outFilterMatchNmin', nargs='?', const=0, type=int)
    parser.add_argument('--outSAMattributes', nargs='?', type=str, default="Standard")
    parser.add_argument('--runThreadN', type=int, default=32)
    parser.add_argument('--limitBAMsortRAM', type=int, default=5784458574)
    parser.add_argument('--sjdbGTFtagExonParentTranscript', default="Parent")
    parser.add_argument('--sjdbOverhang', type=int, default=None)
    parser.add_argument('--limitSjdbInsertNsj', type=int, default=1602710)

    parser.add_argument('--sjdbGTFfeatureExon')
    parser.add_argument('--sjdbGTFtagExonParentGene')
    parser.add_argument('--quantMode', nargs="+")
    parser.add_argument('--salmon_genome_fasta')
    parser.add_argument('--config', help="config file, override everything")
    parser.add_argument("--tmp", default="/tmp")

    # htseq options
    parser.add_argument("--htseqStranded", default='no')
    parser.add_argument("--htseqFeatureType", default="exon")
    parser.add_argument("--htseqID", default="Parent")
    parser.add_argument("--htseqOrder", default="pos")

    args = parser.parse_args()

    now = datetime.datetime.now()
    timeprint = now.strftime("%Y-%m-%d %H:%M")

    if args.config is not None:
        run_config(args.config)
    else:
        data_folder = "%s/%s" % (args.dataroot, args.indir)
        if args.genome_fasta is not None and os.path.exists(args.genome_fasta):
            genome_fasta = args.genome_fasta
        else:
            genome_fasta = glob.glob('%s/*.fasta' % (args.genomedir))[0]

        # make htseq config from comman line args
        config = {
            "htseq_options": {
                "stranded": args.htseqStranded,
                "feature_type": args.htseqFeatureType,
                "id_attribute": args.htseqID,
                "order": args.htseqOrder
            },
            # added cleanup config
            "cleanup_after_run": {
                "intermediate_bam_files": True,
                "trim_files": True,
                "salmon_log": True
            }
        }
        data_trimmed_dir,fastqc_dir,results_dir = run_pipeline(data_folder, args.outdir, args.genomedir, genome_fasta, args, config)
