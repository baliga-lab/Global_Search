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
star_prefix="star_{{star_options.outFilterMismatchNmax}}_{{star_options.outFilterMismatchNoverLmax}}_{{star_options.outFilterScoreMinOverLread}}_{{star_options.outFilterMatchNmin}}{{dedup_prefix}}"
salmon_prefix="salmon_{{star_options.outFilterMismatchNmax}}_{{star_options.outFilterMismatchNoverLmax}}_{{star_options.outFilterScoreMinOverLread}}_{{star_options.outFilterMatchNmin}}{{dedup_prefix}}"

{{sbatch_extras}}

python3 -m globalsearch.rnaseq.run_star_salmon {{star_extra_options}} {{salmon_extra_options}} {{twopass_mode}} {{fastq_patterns}} {{runThreadN}} {{out_sam_attributes}} --outFilterMismatchNmax {{star_options.outFilterMismatchNmax}} --outFilterMismatchNoverLmax {{star_options.outFilterMismatchNoverLmax}} --outFilterScoreMinOverLread {{star_options.outFilterScoreMinOverLread}} --outFilterMatchNmin {{star_options.outFilterMatchNmin}} {{dedup_option}} --starPrefix $star_prefix --salmonPrefix $salmon_prefix {{genome_gff_option}} {{genome_fasta_option}} {{genome_dir}} {{input_dir}} $data_folder {{output_dir}}
"""

DESCRIPTION = """make_star_salmon_job.py - Create STAR Salmon job file for Slurm"""

def make_sbatch_options(config):
    result = ""
    for option in config['sbatch_options']['star_salmon']['options']:
        result += "#SBATCH %s\n" % option
    return result

def make_sbatch_extras(config):
    result = ""
    for extra in config['sbatch_options']['star_salmon']['extras']:
        result += "%s\n" % extra
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
    config['genome'] = genome
    config['sbatch_extras'] = make_sbatch_extras(config)
    config['sbatch_options_comments'] = make_sbatch_options(config)

    config['dedup_prefix'] = '_dedup' if config['deduplicate_bam_files'] else ''
    config['dedup_option'] = '--dedup' if config['deduplicate_bam_files'] else ''
    config['dedup_option'] = '--dedup' if config['deduplicate_bam_files'] else ''
    config['twopass_mode'] = '--twopassMode' if config['star_options']['twopassMode'] else ''

    # override runThreadN
    try:
        config["runThreadN"] = "--runThreadN %d" % config['star_options']['runThreadN']
    except:
        config["runThreadN"] = ''

    # override outSAMattributes
    try:
        out_sam_attributes = config['star_options']['outSAMattributes']
        config['out_sam_attributes'] = '--outSAMattributes %s' % ' '.join(out_sam_attributes) if len(out_sam_attributes) > 0 else ''
    except:
        config['out_sam_attributes'] = ''

    config['fastq_patterns'] = '--fastq_patterns "%s"' % ','.join(config['fastq_patterns']) if len(config['fastq_patterns']) > 0 else ''

    # see if optional genome_gff exists
    try:
        config['genome_gff_option'] = ''
        genome_gff = config['genome_gff']
        if os.path.exists(genome_gff):
            config['genome_gff_option'] = '--genome_gff %s' % genome_gff
        try:
            sjdb_gtf_tag_exon_parent_transcript = config['star_options']['sjdbGTFtagExonParentTranscript']
            config['genome_gff_option'] += (' --sjdbGTFtagExonParentTranscript %s' % sjdb_gtf_tag_exon_parent_transcript)
        except KeyError:
            pass  # ignore if doesn't exist
        try:
            sjdb_overhang = config['star_options']['sjdbOverhang']
            config['genome_gff_option'] += (' --sjdbOverhang %s' % str(sjdb_overhang))
        except KeyError:
            pass  # ignore if doesn't exist
        try:
            limit_sjdb_insert_nsj = config['star_options']['limitSjdbInsertNsj']
            config['genome_gff_option'] += (' --limitSjdbInsertNsj %s' % str(limit_sjdb_insert_nsj))
        except KeyError:
            pass  # ignore if doesn't exist
    except:
        pass

    # More extra options
    star_extra_options = []
    try:
        star_extra_options += ["--sjdbGTFfeatureExon", config['star_options']['sjdbGTFfeatureExon']]
    except:
        pass
    try:
        star_extra_options += ["--sjdbGTFtagExonParentTranscript", config['star_options']['sjdbGTFtagExonParentTranscript']]
    except:
        pass
    try:
        star_extra_options += ["--sjdbGTFtagExonParentGene", config['star_options']['sjdbGTFtagExonParentGene']]
    except:
        pass
    try:
        quantmode = config['star_options']['quantMode']
        star_extra_options += ["--quantMode"]
        star_extra_options += quantmode
    except:
        pass
    config['star_extra_options'] = ' '.join(star_extra_options)

    salmon_extra_options = []
    try:
        salmon_extra_options += ["--salmon_genome_fasta", config['salmon_options']['genome_fasta']]
    except:
        pass
    config['salmon_extra_options'] = ' '.join(salmon_extra_options)

    # see if optional genome_fasta exists
    try:
        config['genome_fasta_option'] = ''
        genome_fasta = config['genome_fasta']
        if os.path.exists(genome_fasta):
            config['genome_fasta_option'] = '--genome_fasta %s' % genome_fasta
    except:
        pass

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
