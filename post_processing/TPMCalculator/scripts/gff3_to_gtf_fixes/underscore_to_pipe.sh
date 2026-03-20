#!/bin/bash
# Convert underscore naming convention to piped for Host

GTF="Valid_Pdae_Cgor_Dtre_Smic_concat.gtf"
FIXED_GTF="Valid_Pdae_Cgor_Dtre_Smic_concat_fixed.gtf"

# Replace '_size' to '|size'
awk -F'\t' -v OFS='\t' '
    /^#/ {print; next} 
    {
        # fix first column
        sub(/_size/, "|size", $1)

        # combine everything after column 8 into column 9
        if (NF > 9) {
            $9 = "";
            for(i=9; i<=NF; i++) {
                $9 = $9 $i
                if(i<NF) $9 = $9 ";"
            }
        }

        print
    }' "$GTF" > "$FIXED_GTF"

echo "Fixed GTF written to $FIXED_GTF"