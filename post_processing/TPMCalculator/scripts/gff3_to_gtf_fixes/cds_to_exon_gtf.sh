#!/bin/bash
# Convert CDS to exon for host so TPMCalculator can count them

GTF="Valid_Plob_Cgor_concat.gtf"
FIXED_GTF="Valid_Plob_Cgor_concat_fixed.gtf"

awk 'BEGIN{OFS="\t"}
{
    # Split line into first 8 fields + rest as attributes
    attrs = ""
    if(NF>8) {
        attrs = $9
        for(i=10;i<=NF;i++) attrs = attrs " " $i
    }

    if($3=="CDS") {
        $3="exon"
        # Extract transcript_id
        match(attrs,/transcript_id "[^"]+"/)
        tid=substr(attrs,RSTART+15,RLENGTH-16)  # get only the value inside quotes
        # Prepend gene_id
        attrs="gene_id \"" tid "\"; " attrs
    }

    # Rebuild the line: columns 1-8 + fixed attrs
    print $1,$2,$3,$4,$5,$6,$7,$8,attrs
}' "$GTF" > "$FIXED_GTF"

echo "Fixed GTF created: $FIXED_GTF"