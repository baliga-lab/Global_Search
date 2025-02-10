#!/usr/bin/env python3 

"""
Script to check parameters and prepare for submission
"""
import argparse
import json
import os, sys, glob, subprocess
from rpy2.robjects.packages import importr, PackageNotInstalledError
import rpy2


DESCRIPTION = """gs_prepare.py - prepare data for workflow submission"""

OUTSAM_ATTRS_STD = {
    "NH", "HI", "AS", "nM", "NM", "MD", "jM", "jI", "XS", "MC", "ch"
}
OUTSAM_ATTRS_EXT = { "vA", "vG", "vW", "CR", "CY", "UR", "UY" }
OUTSAM_ATTRS_EXT2 = { "rB", "vR" }
OUTSAM_ATTRS_SPECIAL = { "None", "Standard", "All" }
OUTSAM_ATTRS_MULTI = OUTSAM_ATTRS_STD | OUTSAM_ATTRS_EXT | OUTSAM_ATTRS_EXT2
OUTSAM_ATTRS_SINGLE = OUTSAM_ATTRS_STD | OUTSAM_ATTRS_EXT | OUTSAM_ATTRS_EXT2 | OUTSAM_ATTRS_SPECIAL

def silent_rpy2_print(o):
    pass

rpy2.rinterface_lib.callbacks.consolewrite_print = silent_rpy2_print
rpy2.rinterface_lib.callbacks.consolewrite_warnerror = silent_rpy2_print


def check_star_options(star_options):
    try:
        outsam_attrs = star_options['outSAMattributes']
        if not isinstance(outsam_attrs, list):
            raise TypeError("STAR options: outSAMattributes expected type: 'list' but was '%s'" % type(outsam_attrs))
        if len(outsam_attrs) == 1:
            if not outsam_attrs[0] in OUTSAM_ATTRS_SINGLE:
                raise ValueError("STAR options (single): outSAMattributes: '%s'" % [outsam_attrs[0]])
        else:
            for attr in outsam_attrs:
                if not attr in OUTSAM_ATTRS_MULTI:
                    raise ValueError("STAR options (multiple): outSAMattributes: '%s'" % attr)
    except KeyError:
        # not specified -> doesn't matter
        pass


def check_params(config, rna_algo):
    """ensures integrity of the parameters"""
    print("checking integrity of parameters...", end="", flush=True)
    if not os.path.exists(config["input_dir"]):
        sys.exit("Input directory '%s' does not exist" % config["input_dir"])

    try:
        if len(config['genome_fasta'].strip()) > 0 and not os.path.exists(config["genome_fasta"]):
            sys.exit("ERROR: Specified Genome FASTA '%s' does not exist" % config["genome_fasta"])
    except:
        print("WARNING: genome_gff not specified (htseq will be skipped)", flush=True)

    if not os.path.exists(config["genome_dir"]):
        sys.exit("Genome directory '%s' does not exist" % config["genome_dir"])
        genome_fasta = glob.glob('%s/*.fasta' % (args.genomedir))
        if len(genome_fasta) == 0:
            sys.exit("Genome directory '%s' does not contain any FASTA files" % config["genome_dir"])

    if len(config["includes"]) > 0:
        # check the existence of the included sample directories.
        # those directories are potentially several levels deep
        for incl in config["includes"]:
            matches = glob.glob(os.path.join(config['input_dir'], "**", incl),
                                recursive=True)
            if len(matches) == 0:
                sys.exit("Input directory '%s' does not exist" % incl)
            else:
                print(matches)
    try:
        check_star_options(config['star_options'])
    except ValueError as e:
        sys.exit(str(e))
    except KeyError:
        pass

    print("done.", flush=True)


def create_dirs(config):
    if not os.path.exists(config["output_dir"]):
        print("creating output directory '%s'" % config["output_dir"], flush=True)
        os.makedirs(config["output_dir"])
    if not os.path.exists(config["log_dir"]):
        print("creating log directory '%s'" % config["log_dir"], flush=True)
        os.makedirs(config["log_dir"])


def __check_command(command, num_info_components=1, version_index=-1,
                    check_version=None, version_switch='--version',
                    multiline=False, info_line=0, fail_if_not_exists=True):
    """Generic command checker, can check for different version info formats and
    restrict version numbers
    """
    try:
        print("checking for %s... " % command, end="", flush=True)
        compl_proc = subprocess.run([command, version_switch], check=True, capture_output=True)
        if multiline:
            info_string = compl_proc.stdout.decode('utf-8').split('\n')[info_line].strip()
        else:
            info_string = compl_proc.stdout.decode('utf-8').strip()
        comps = info_string.split()
        version = comps[version_index]

        print("(found version '%s') ..." % version, end="", flush=True)
        if check_version is not None and check_version != version:
            sys.exit("Unsupported version %s. Currently, only %s %s is supported" % (command, check_version))
        print("done.", flush=True)
    except FileNotFoundError:
        if fail_if_not_exists:
            sys.exit("Can not find %s (not installed or not in PATH)" % command)
        else:
            print("WARN: %s does not exist, but is optional" % command, flush=True)


def check_salmon():
    #__check_command("salmon"), num_info_components=2, check_version="0.13.1")
    __check_command("salmon")


def check_star():
    __check_command("STAR")


def check_htseq():
    __check_command("htseq-count")


def check_samtools():
    __check_command("samtools", multiline=True, num_info_components=2)


def check_kallisto():
    __check_command("kallisto", version_switch="version", num_info_components=3,
                    fail_if_not_exists=False)

def check_trim_galore():
    __check_command("trim_galore", num_info_components=2, multiline=True, info_line=3)


def check_rlibrary_installed(libname):
    try:
        lib = importr(libname)
        return True
    except PackageNotInstalledError:
        return False


def check_rlibraries_installed():
    print("checking for installed R libraries...", end="", flush=True)
    not_installed = []
    for rlib in ['GlobalSearch']:
        if not check_rlibrary_installed(rlib):
            not_installed.append(rlib)
    print("done", flush=True)

    return not_installed

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('configfile', help="configuration file")
    args = parser.parse_args()

    with open(args.configfile) as infile:
        config = json.load(infile)
    rna_algo = config['rnaseq_algorithm']

    if rna_algo == 'star_salmon':
        check_star()
        check_salmon()
    if rna_algo == 'kallisto':
        check_kallisto()

    check_htseq()
    check_samtools()
    check_trim_galore()
    not_installed = check_rlibraries_installed()
    if len(not_installed) > 0:
        for libname in not_installed:
            print("ERROR: Please install the R library '%s'" % libname)
        sys.exit()

    check_params(config, rna_algo)
    create_dirs(config)
    print(rna_algo)
