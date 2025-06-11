#!/bin/bash\
#
# pfam hmmr hits
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file: {run}/pfam/All_bins_pfam_sql.out
echo $1
# database file: config["pfam_hmmr"]["database_names"]
echo $2
# output file: {run}/sql/pfam.tsv
echo $3

cat $2 | awk -F"\t" 'BEGIN{OFS="\t"}
 FNR==NR{gsub("\r","",$0);h[$4]=$1;next}
 {
 print $1, h[$2], $3, $4, $5, $6;
 }' - $1 > $3
