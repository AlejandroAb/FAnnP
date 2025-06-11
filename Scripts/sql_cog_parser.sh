#!/bin/bash\
#
# cog hmmr hits
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file: {run}/cog/All_bins_cog_sql.out
echo $1
# output file: {run}/sql/cog.tsv
echo $2

cat $1 | awk -F"\t" 'BEGIN{OFS="\t"}
 {
 print $1, $2, $3, $4, $5, $6;
 }' > $2
