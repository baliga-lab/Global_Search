import os
from fs.osfs import OSFS
import fs

TRIMGALORE_SUFFIX_SINGLE = "_trimmed.fq"
TRIMGALORE_SUFFIX_PAIRED = "_val_1.fq"
PAIRED_END_FQ_PATTERN = '/*_val_%d.fq'
SINGLE_END_FQ_PATTERN = '/*_trimmed.fq'

"""
trim_galore --fastqc_args "--outdir /proj/omics4tb2/wwu/Global_Search/redsea-output/R1/fastqc_results/" --paired --output_dir /proj/omics4tb2/wwu/Global_Search/redsea-output/R1/trimmed/ /proj/omics4tb2/wwu/GlobalSearch.old/Pilot_Fail_Concat/rawdata/R1/R1_concat_1.fq.gz /proj/omics4tb2/wwu/GlobalSearch.old/Pilot_Fail_Concat/rawdata/R1/R1_concat_2.fq.gz"""

def trim_galore(first_pair_file, second_pair_file, folder_name, sample_id, data_trimmed_dir,
                fastqc_dir, num_cores=1):
    paired = second_pair_file is not None  # make sure we are not single end

    # check whether the result already exists and skip if it does
    file_base = os.path.basename(first_pair_file).replace(".gz", "").replace(".fastq", "").replace(".fq", "")
    if paired:
        trimmed_file_name = file_base + TRIMGALORE_SUFFIX_PAIRED
    else:
        trimmed_file_name = file_base + TRIMGALORE_SUFFIX_SINGLE

    trimmed_path = os.path.join(data_trimmed_dir, trimmed_file_name)
    trimmed_path_gz = trimmed_path + ".gz"
    if os.path.exists(trimmed_path) or os.path.exists(trimmed_path_gz):
        print("Trimmed file '%s' found, skipping trim_galore" % trimmed_file_name, flush=True)
        return

    print ("\033[34m Running TrimGalore \033[0m")
    # create sample specific trimmed and fastqc directories
    if not os.path.exists(data_trimmed_dir):
        os.makedirs(data_trimmed_dir)
    if not os.path.exists(fastqc_dir):
        os.makedirs(fastqc_dir)

    # run Command
    command = ['trim_galore', '--fastqc_args', '"--outdir %s/"' % fastqc_dir]
    if num_cores > 1:
        command += ['--j', '%d' % num_cores]

    if paired:
        command.append('--paired')
    command += ['--output_dir', data_trimmed_dir, first_pair_file]
    if second_pair_file is not None:
        command.append(second_pair_file)
    cmd = ' '.join(command)
    print( '++++++ Trimgalore Command:', cmd)
    os.system(cmd)  # run with os.system(), since you have that funny outdir parameter


####################### Collect trimmed data files ###############################


def collect_trimmed_data(data_trimmed_dir, is_gzip, is_paired_end, rootfs=OSFS("/")):
    """
    Collect trimmed data files from the specified trimmed directory
    """
    filesys = rootfs.opendir(data_trimmed_dir)
    # define result files
    pattern = PAIRED_END_FQ_PATTERN if is_paired_end else SINGLE_END_FQ_PATTERN
    if is_gzip:
        pattern += ".gz"

    if is_paired_end:
        first_pair_trimmed, second_pair_trimmed = [filesys.glob(pattern % i) for i in [1, 2]]
    else:
        first_pair_trimmed = filesys.glob(pattern)
        second_pair_trimmed = None

    first_pair_trimmed = [fs.path.combine(data_trimmed_dir, match.path)
                          for match in first_pair_trimmed]
    if second_pair_trimmed is not None:
        second_pair_trimmed = [fs.path.combine(data_trimmed_dir, match.path)
                               for match in second_pair_trimmed]
    else:
        second_pair_trimmed = []

    first_pair_group = ' '.join(first_pair_trimmed)
    second_pair_group = ' '.join(second_pair_trimmed)

    return first_pair_group, second_pair_group


def create_result_dirs(data_trimmed_dir, fastqc_dir, results_dir, htseq_dir):
    dirs = [data_trimmed_dir, fastqc_dir, results_dir, htseq_dir]
    for dir in dirs:
        # create results folder
        if not os.path.exists('%s' %(dir)):
            os.makedirs('%s' %(dir))
        else:
            print('\033[31m %s directory exists. Not creating. \033[0m' %(dir))
