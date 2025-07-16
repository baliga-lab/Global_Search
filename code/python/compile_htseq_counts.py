#!/usr/bin/env python3

import pandas
import sys, os
import fileinput
import json
from collections import defaultdict

"""
Example:
find . -name "*_htseqcounts.txt" | xargs realpath | /proj/omics4tb2/wwu/Global_Search/code/python/compile_htseq_counts.py ../../../gbr-starhtseq-Ahya_Cgor.json

for Plob_Cgor
find . -name "*_count.txt" | xargs realpath | /proj/omics4tb2/wwu/Global_Search/code/python/compile_htseq_counts.py ../../../pacific-starhtseq-Plob_Cgor.json
"""
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: {0} <config-file>".format(sys.argv[0]))
        exit(1)

    with open(sys.argv[1]) as infile:
        config = json.load(infile)

    organisms = config['organisms']
    samples = {}
    gene2sample = {}
    #organism_samples = defaultdict(list)
    org_groups = set()
    organism_genes = defaultdict(list)

    for line in sys.stdin:
        path = line.strip()
        fname = os.path.basename(path)
        if fname == 'host_count.txt' or fname == 'sym_count.txt':
            sample = path.split('/')[-4]
        else:
            sample = fname.replace("_htseqcounts.txt", "")
        samples[sample] = {}
        print("Process sample ", sample)
        #print(path)
        with open(path) as infile:
            for data_row in infile:
                gene, htseq_count = data_row.strip().split()
                samples[sample][gene] = int(htseq_count)
                gene2sample[gene] = sample
                if gene.startswith("Host_") or gene.startswith("Sym_"):
                    organism = '_'.join(gene.split("_")[:2])
                    #organism_samples[organism].append(sample)
                    org_groups.add(organism)
                    organism_genes[organism].append(gene)

    samples_sorted = sorted(samples.keys())
    genes_sorted = sorted(gene2sample.keys())

    # we probably don't need the organism list
    # first organism is the host
    # the symbionts are all the following
    organism_string = '_'.join(organisms)
    header = ["gene_id"]
    header.extend(list(samples_sorted))
    num_key_errors = 0
    missing_sample_genes = defaultdict(list)

    with open("STAR_htseq_%s_NumReads_Matrix_Merged.csv" % organism_string, "w") as outfile0:
        outfile0.write(",".join(header) + "\n")
        for org_group in org_groups:
            with open("STAR_htseq_%s_NumReads_Matrix.csv" % org_group, "w") as outfile:
                outfile.write(",".join(header) + "\n")
                #genes = sorted(organism_genes[org_group])
                genes = sorted(set(organism_genes[org_group]))
                for gene in genes:
                    out_row = [gene]
                    for sample in samples_sorted:
                        # append the count
                        try:
                            out_row.append(str(samples[sample][gene]))
                        except KeyError:
                            # gene not found, so we append 0
                            num_key_errors += 1
                            missing_sample_genes[sample].append(gene)
                            out_row.append("0")
                    outfile.write(",".join(out_row) + "\n")
                    outfile0.write(",".join(out_row) + "\n")  # write to merged

    print("# key errors: ", num_key_errors)
    with open("missing_sample_genes.json", "w") as outfile:
        json.dump(missing_sample_genes, outfile)
