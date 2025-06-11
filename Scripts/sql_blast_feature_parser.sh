#!/bin/bash\
#
# custom database blast hits
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file: {run}/db{number}_blast/All_bins_db{number}_blast.out
echo $1
# database id: 10 or 11
echo $2
# output file: {run}/sql/db{number}_blast.tsv
echo $3

cat $1 | awk -v database_id=$2 -F"\t" 'BEGIN{OFS="\t"}
 {
 print database_id, $1, $7, $8, $2, "\\N", $11, $12, $3;
 }' > $3
