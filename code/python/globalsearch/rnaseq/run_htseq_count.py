"""run_htseq - module for running htseq-count"""

import os

class HtseqArgs:
    """Class that fakes a command line argument object"""
    def __init__(self, config):
        """create args from configuration file"""
        htseq_options = config['htseq_options']
        self.htseqStranded = htseq_options['stranded']
        self.htseqFeatureType = htseq_options['feature_type']
        self.htseqID = htseq_options['id_attribute']
        self.htseqOrder = htseq_options['order']



def run_htseq_count(final_bam, htseq_resultdir,
                    folder_name, genome_gff, args):
    """ Run htseq-count """
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
