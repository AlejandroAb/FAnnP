#!/bin/bash\
#
# ko hmmr hits
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file: {run}/kfam/All_bins_ko_sql.out
echo $1
# output file: {run}/sql/ko.tsv
echo $2

cat $1 | awk -F"\t" 'BEGIN{OFS="\t"}
 {
 print $1, $2, $3, $4, $5, $6;
 }' > $2
