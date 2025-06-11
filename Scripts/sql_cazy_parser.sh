#!/bin/bash\
#
# cazy hmmr hits
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file: {run}/cazy/All_bins_cazy_sql.out
echo $1
# database file: config["cazy_hmmr"]["database_names"]
echo $2
# database id: 03
echo $3
# output file: {run}/sql/cazy.tsv
echo $4

cat $2 | awk -v database_id=$3 -F"\t" 'BEGIN{OFS="\t"}
 FNR==NR{gsub("\r","",$0);h[$1]=$2;next}
 {
 gsub("'\''","\\'\''",h[$2]);
 print database_id, $1, $3, $4, $2, h[$2], $5, $6, "\\N";
 }' - $1 > $4
