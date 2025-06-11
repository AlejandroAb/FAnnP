#!/bin/bash\
#
# arCOG hmmr hits
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file: {run}/arCOG/All_bins_arCOG_sql.out
echo $1
# database file: config["arCOG_hmmr"]["database_names"]
echo $2
# database id: 01
echo $3
# output file: {run}/sql/arCOG.tsv
echo $4

cat $2 | awk -v database_id=$3 -F"\t" 'BEGIN{OFS="\t"}
 FNR==NR{gsub("\r","",$0);h[$1]=$2" - "$3" - "$4;next}
 {
 split($2,a,".")
 gsub("'\''","\\'\''",h[a[1]])
 print database_id, $1, $3, $4, $2, h[a[1]], $5, $6, "\\N";
 }' - $1 > $4
