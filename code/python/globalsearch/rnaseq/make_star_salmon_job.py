#!/usr/bin/env python3

import jinja2
import os
import argparse
import json

from globalsearch.rnaseq.find_files import rnaseq_data_folder_list

TEMPLATE = """#!/bin/bash

#SBATCH -J star_salmon_{{genome}}
#SBATCH -o {{log_dir}}/"%A"."%a".out
#SBATCH -e {{log_dir}}/"%A"."%a".out
#SBATCH --array={{array_range}}

{{sbatch_options_comments}}

echo "ARRAY TASK ID: $SLURM_ARRAY_TASK_ID"
data_folders=({{data_folders}})
data_folder=${data_folders[$SLURM_ARRAY_TASK_ID]}

{{sbatch_extras}}

time python3 -m globalsearch.rnaseq.run_star_salmon {{config_file}} $data_folder
"""

DESCRIPTION = """make_star_salmon_job.py - Create STAR Salmon job file for Slurm"""

def make_sbatch_options(config):
    result = ""
    try:
        for option in config['sbatch_options']['options']:
            result += "#SBATCH %s\n" % option
    except KeyError:
        pass
    return result

def make_sbatch_extras(config):
    result = ""
    try:
        for extra in config['sbatch_options']['extras']:
            result += "%s\n" % extra
    except KeyError:
        pass
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('configfile', help="configuration file")
    args = parser.parse_args()
    with open(args.configfile) as infile:
        config = json.load(infile)

    templ = jinja2.Template(TEMPLATE)
    genome = os.path.basename(os.path.normpath(config['genome_dir']))
    config['config_file'] = args.configfile
    config['genome'] = genome
    config['sbatch_extras'] = make_sbatch_extras(config)
    config['sbatch_options_comments'] = make_sbatch_options(config)

    data_folders = ['"%s"' % f for f in rnaseq_data_folder_list(config)]
    config["data_folders"] = ' '.join(data_folders)

    # Array specification
    try:
        array_max_tasks = config['sbatch_options']['array_max_tasks']
    except:
        array_max_tasks = 0
    if array_max_tasks > 0:
        array_max_task_spec = "%%%d" % array_max_tasks
    else:
        array_max_task_spec = ""

    config["array_range"] = "0-%d%s" % (len(data_folders) - 1, array_max_task_spec)
    print(templ.render(config))
