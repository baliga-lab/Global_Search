#!/bin/bash
# Fix GTF chromosome names to match BAM headers for TPMCalculator
# Original issue: BAM and GTF scaffold names did not match, causing TPMCalculator to run without producing gene output.

INPUT_GTF="Valid_Phar_Cgor_Dtre_Smic_concat.gtf"
OUTPUT_GTF="Valid_Phar_Cgor_Dtre_Smic_concat_fixed.gtf"

awk 'BEGIN{FS=OFS="\t"} 
{
  match($1, /JBDLLT010000[0-9]+/)
  if (RSTART) {
    id=substr($1, RSTART, RLENGTH)
    sub(/JBDLLT010000/, "", id)
    $1="Host_Phar_gnl|WGS:JBDLLT|Phar_scaff_" id "|gb|JBDLLT010000" id
  }
  print
}' "$INPUT_GTF" > "$OUTPUT_GTF"

echo "Fixed GTF written to $OUTPUT_GTF"