# GFF3 to GTF Fixes

This directory contains small helper scripts used to make manual fixes to GTF files after converting annotations from **GFF3 → GTF**.

While the main conversion script performs the bulk of the format conversion, some reference annotations required small adjustments to ensure compatibility with downstream tools such as **TPMCalculator**. These fixes were typically related to feature naming or attribute formatting.

Examples of fixes performed:
- Converting CDS to exon label
- Fixing GTF chromosome names to match BAM headers

Each script documents the specific transformation it performs. These scripts were run **after the initial GFF3 → GTF conversion step** and before running downstream expression quantification.