#!/usr/bin/env python3
"""
gff3_to_gtf.py

This script parses a GFF3 file and converts supported feature types
(gene, mRNA/transcript, exon, CDS) into GTF format while maintaining
the correct gene_id and transcript_id relationships.

Functions:
    - parse_attributes(attr_string)
    - format_gtf_attributes(attributes, feature_type, gene_id=None, transcript_id=None)
    - convert_gff3_to_gtf(input_file, output_file)
    - main()

Usage:
    python3 gff3_to_gtf.py <input_gff3> <output_gtf>
"""

import sys
import re
from pathlib import Path

def parse_attributes(attr_string):
    """
    GFF3 attributes are formatted as key-value pairs separated by
    semicolons. This function converts the attribute string into a
    Python dictionary for easier access.

    Args:
        - attr_string (str): Attribute column from a GFF3 file.

    Returns:
        - attributes (dict): Dictionary mapping attribute keys to values.
    """
    attributes = {}
    if not attr_string or attr_string == '.':
        return attributes
    
    for item in attr_string.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            attributes[key] = value
    return attributes

def format_gtf_attributes(attributes, feature_type, gene_id=None, transcript_id=None):
    """
    Constructs the GTF attribute string using gene_id and transcript_id
    along with optional gene_name or transcript_name fields if present.

    Args:
        - attributes (dict): Parsed GFF3 attributes.
        - feature_type (str): Feature type (gene, transcript, exon, CDS).
        - gene_id (str, optional): Gene identifier.
        - transcript_id (str, optional): Transcript identifier.

    Returns:
        - str: Properly formatted GTF attribute string.
    """
    gtf_attrs = []
    
    if gene_id:
        gtf_attrs.append(f'gene_id "{gene_id}"')
    
    if transcript_id:
        gtf_attrs.append(f'transcript_id "{transcript_id}"')
    
    if 'Name' in attributes:
        if feature_type == 'gene':
            gtf_attrs.append(f'gene_name "{attributes["Name"]}"')
        elif feature_type in ['transcript', 'mRNA']:
            gtf_attrs.append(f'transcript_name "{attributes["Name"]}"')
    
    return '; '.join(gtf_attrs) + ';'

def convert_gff3_to_gtf(input_file, output_file):
    """
    Reads a GFF3 file line by line, parses feature attributes, and writes
    supported features (gene, transcript/mRNA, exon, CDS) to a GTF file.
    Transcript-to-gene relationships are tracked to correctly assign
    gene_id and transcript_id values.

    Args:
        - input_file (str): Path to the input GFF3 file.
        - output_file (str): Path to the output GTF file.

    Returns:
        - N/A
    """
    transcript_to_gene = {}
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line_num, line in enumerate(infile, 1):
            if line_num % 100000 == 0:
                print(f"Processed {line_num:,} lines...")
            
            line = line.strip()
            
            if line.startswith('#') or not line:
                continue
            
            fields = line.split('\t')
            if len(fields) != 9:
                continue
            
            seqname, source, feature, start, end, score, strand, frame, attributes = fields
            attr_dict = parse_attributes(attributes)
            
            if feature == 'gene':
                gene_id = attr_dict.get('ID', '')
                
                gtf_attributes = format_gtf_attributes(attr_dict, feature, gene_id=gene_id)
                gtf_line = f"{seqname}\t{source}\tgene\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{gtf_attributes}"
                outfile.write(gtf_line + '\n')
            
            elif feature in ['mRNA', 'transcript']:
                transcript_id = attr_dict.get('ID', '')
                parent_gene = attr_dict.get('Parent', '')
                
                transcript_to_gene[transcript_id] = parent_gene
                
                gtf_attributes = format_gtf_attributes(attr_dict, feature, gene_id=parent_gene, transcript_id=transcript_id)
                gtf_line = f"{seqname}\t{source}\ttranscript\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{gtf_attributes}"
                outfile.write(gtf_line + '\n')
            
            elif feature in ['exon', 'CDS']:
                parent_id = attr_dict.get('Parent', '')
                
                # Check if parent is a transcript/mRNA or directly a gene
                if parent_id in transcript_to_gene:
                    # Parent is a transcript
                    parent_gene = transcript_to_gene[parent_id]
                    transcript_id = parent_id
                else:
                    # Parent might be a gene (some exons directly reference genes)
                    parent_gene = parent_id
                    transcript_id = parent_id  # Use gene ID as transcript ID for these cases
                
                gtf_attributes = format_gtf_attributes(attr_dict, feature, gene_id=parent_gene, transcript_id=transcript_id)
                gtf_line = f"{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{gtf_attributes}"
                outfile.write(gtf_line + '\n')

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 gff3_to_gtf.py <input_gff3> <output_gtf>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not Path(input_file).exists():
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    try:
        convert_gff3_to_gtf(input_file, output_file)
        print(f"Successfully converted {input_file} to {output_file}")
    except Exception as e:
        print(f"Error during conversion: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()