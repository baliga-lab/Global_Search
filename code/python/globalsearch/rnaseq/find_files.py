#!/usr/bin/env python3
"""
find_files.py - module for flexible finding of FASTQ files
"""
from collections import defaultdict
import glob
import logging
import os, re, fs
import jinja2
from fs.osfs import OSFS  # make sure we install the fs package !!!


def _get_inner_dirs(filesys, curpath):
    subdirs = [d for d in filesys.listdir(curpath)
               if filesys.isdir(fs.path.combine(curpath, d))]
    if len(subdirs) == 0:
        return [curpath]

    sdds = [_get_inner_dirs(filesys, fs.path.combine(curpath, sd))
            for sd in subdirs]
    flattened = []
    for elem in sdds:
        if isinstance(elem, list):
            flattened.extend(elem)
        else:
            flattened.append(elem)

    return flattened


def rnaseq_data_folder_list(config, filesys=OSFS('/')):
    """Function to determine the list of directories that are to be submitted to
    the cluster for RNA sequencing analysis
    """
    result = []
    #pattern = re.compile('R[E]?\d+.*')
    try:
        includes_list = config['includes']
    except KeyError:
        includes_list = []
    try:
        include_path = config['include_file']
        if len(include_path) > 0 and filesys.exists(include_path):
            with filesys.open(include_path) as infile:
                for line in infile:
                    includes_list.append(line.strip())
        elif len(include_path) > 0 and not filesys.exists(include_path):
            logging.warning("include file '%s' does not exist", include_path)
    except:
        # ignore any errors in the include_file
        pass

    input_dir = config['input_dir']
    inner_dirs = _get_inner_dirs(filesys, input_dir)
    inner_dirs = [d.replace(input_dir, "").lstrip('/') for d in inner_dirs]

    if len(includes_list) > 0:
        includes_set = set(includes_list)
        result = [d for d in inner_dirs
                  if os.path.basename(d) in includes_set]
    else:
        # we simply return inner dirs, but we might think about
        # checking the last part against the pattern
        # take the top level directories in the input directory
        # that match the pattern
        #print("ADD TOPLEVEL FILES: %s" % config['input_dir'])
        #result = [d for d in filesys.listdir(config['input_dir']) if re.match(pattern, d)]
        #result = [d for d in filesys.listdir(config['input_dir'])
        #        if filesys.isdir(fs.path.combine(config['input_dir'], d))]
        result = inner_dirs
    return result


def _find_fastq_files(data_folder, patterns, rootfs):
    """Stable version of file finder, search multiple patterns"""
    logger = logging.getLogger("rnaseq")
    result = []
    filesys = rootfs.opendir(data_folder)
    for pat in patterns:
        templ = jinja2.Template(pat)
        patmatch = defaultdict(list)
        for read_num in [1, 2]:
            # keep "pairnum" for legacy reasons
            pattern_replace_readnum = templ.render({'pairnum': read_num, 'readnum': read_num})
            pattern = fs.path.combine('/**', pattern_replace_readnum)
            #print("SEARCHING PATTERN: '%s'" % pattern, flush=True)
            if pattern_replace_readnum == pat and read_num == 2:  # no readnum provided in pattern
                patmatch[read_num].append(None)
            else:
                for match in filesys.glob(pattern):
                    patmatch[read_num].append(fs.path.combine(data_folder, match.path))
        #print("PATMATCH: ", patmatch)
        try:
            first = patmatch[1][0]
        except:
            # no match -> skip this pattern
            continue
        try:
            second = patmatch[2][0]
        except:
            second = None
        result.append((first, second))
    return result


def find_fastq_files(data_folder, patterns, filesys=OSFS('/')):
    """
    This function finds the FASTQ files according to the specified patterns. It will
    start at data_folder and try to find all FASTQ files, if possible as pairs.
    If there is no second file, the second component of the result will be None.
    For the most part, patterns will follow the glob format, you can specify the position
    of the pair number in Jinja2 format
    Example:

    A pattern of

    *_{{readnum}}.fastq.gz

    with a readnum of 1 would "match myfile_1.fastq.gz"


    :param data_folder: the top level folder to start searching from
    :param patterns: list of patterns to use in glob searching
    :param filesys: the PyFS file system to use in glob searching
    :return the list of matching paths
    """
    return _find_fastq_files(data_folder, patterns, filesys)


# Base pattern for searching
#DATA_SEARCH1 = '%s/RNA/*R1*.fastq*' % data_folder # test
#DATA_SEARCH1 = '%s/*_1.fq*' % data_folder
DATA_SEARCH1 = "/RNA/*[_.R]1*.fastq.gz"
DATA_SEARCH2 = "/RNA/*/*[_.R]1*.fastq.gz"

def find_fastq_files_experimental(data_folder, pattern1=DATA_SEARCH1, pattern2=DATA_SEARCH2):
    # get logger object
    logger = logging.getLogger("rnaseq")

    # Get the list of first file names in paired end sequences
    # Search first set of folders
    first_pair_files1 = glob.glob(data_folder + pattern1, recursive=True)
    # search additional folders if there are different levels
    first_pair_files2 = glob.glob(data_folder + pattern2, recursive=True)

    # combine search results
    first_pair_files = first_pair_files1 + first_pair_files2

    if len(first_pair_files) == 0:
        logger.error("Error: I didnt find any file with '*fastq*' extension in %s\n" % (data_folder))

    """
    # Program specific results directories
    data_trimmed_dir = "%s/trimmed" % (results_folder)
    fastqc_dir = "%s/fastqc_results" % (results_folder)

    results_dir = "%s/results_STAR" %(results_folder)
    rsem_results_dir = "%s/results_RSEM" %(results_folder) 

    # Run create directories function to create directory structure
    create_dirs(data_trimmed_dir, fastqc_dir, results_dir, rsem_results_dir)

    print("FIRST_PAIR_FILES: ", first_pair_files)

    # Loop through each file and create filenames
    file_count = 1
    for first_pair_file in first_pair_files:
        first_file_name_full = first_pair_file.split('/')[-1]

        if re.search("_R1_001.fastq.gz", first_pair_file):
            print("Found a match for: '_R1_001.fastq.gz'\n")
            second_pair_file = re.sub("_R1_001.fastq.gz", "_R2_001.fastq.gz", first_pair_file)
        elif re.search(".R1.fastq.gz", first_pair_file):
            print("Found a match for: '.R1.fastq.gz'\n")
            second_pair_file = re.sub(".R1.fastq.gz", ".R2.fastq.gz", first_pair_file)
        elif re.search("_1.fastq.gz", first_pair_file):
            print("Found a match for: '_1.fastq.gz'\n")
            second_pair_file = re.sub("_1.fastq.gz", "_2.fastq.gz", first_pair_file)
        elif re.search(".1.fastq.gz", first_pair_file):
            print("Found a match for: '.1.fastq.gz'\n")
            second_pair_file = re.sub(".1.fastq.gz", ".2.fastq.gz", first_pair_file)
        else:
            print("No match found for any defined types\n")
            with open(error_file,'w') as f:
                f.write("No match found for any defined types\n")
                f.write(str(first_pair_file))

        print(first_pair_file)
        print(second_pair_file)

        second_file_name_full = second_pair_file.split('/')[-1]
        file_ext = first_pair_file.split('.')[-1]
        print()
        print ('\033[32m Processing File: %s of %s (%s)\033[0m' %(file_count, len(first_pair_files), first_file_name_full ))

        first_file_name = re.split('.fq|.fq.gz',first_file_name_full)[0]
        second_file_name = re.split('.fq|.fq.gz',second_file_name_full)[0]
        print('first_file_name:%s, second_file_name:%s' %(first_file_name,second_file_name))
"""

DATA_FOLDER = "/proj/omics4tb2/wwu/GlobalSearch.old/Pilot_Fail_Concat/rawdata"
SEARCH_PATTERN = '/*_1.fq*'
if __name__ == '__main__':
    files = find_fastq_files(os.path.join(DATA_FOLDER, 'R1'), SEARCH_PATTERN)
    print(files)
 
