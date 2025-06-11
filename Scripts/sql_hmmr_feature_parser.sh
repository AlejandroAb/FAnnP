#!/bin/bash\
#
# custom database hmmr hits
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file: {run}/db{number}_hmmr/All_bins_db{number}_hmmr_sql.out
echo $1
# database id: 08 or 09
echo $2
# output file: {run}/sql/db{number}_hmmr.tsv
echo $3

cat $1 | awk -v database_id=$2 -F"\t" 'BEGIN{OFS="\t"}
 {
 print database_id, $1, $3, $4, $2, "\\N", $5, $6, "\\N";
 }' > $3
