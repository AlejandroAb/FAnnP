#!/bin/bash\
#
# transporter blastp hits
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file: {run}/transporter/All_bins_transporter.out
echo $1
# database file: config["transporterDB_blast"]["database_names"]
echo $2
# database id: 06
echo $3
# output file: {run}/sql/transporter.tsv
echo $4

cat $2 | iconv -c -f utf-8 -t ascii | awk -v database_id=$3 -F"\t" 'BEGIN{OFS="\t"}
 FNR==NR{gsub("\r","",$0);a[$1]=$2" - "$3;next}
 {
 gsub("'\''","\\'\''",a[$2]);
 print database_id, $1, $7, $8, $2, a[$2], $11, $12, $3;
 }' - $1 > $4


