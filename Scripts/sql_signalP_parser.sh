#!/bin/bash\
#
# signalP hits
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file: {run}/signalP/prediction_summary.txt
echo $1
# database id: 04
echo $2
# output file: {run}/sql/signalP.tsv
echo $3

cat $1 | awk -v database_id=$2 -F"\t" 'BEGIN{OFS="\t"}
 {
 split($4, a, ".");
 split(a[1], b, ":");
 gsub(" ", "", b[2]);
 split(b[2], c, "-");
 split($4, d, "Pr:");
 gsub(" ", "", d[2]);
 print database_id, $1, c[1], c[2], $2, "\\N", $3, "\\N", d[2];
 }' > $3

