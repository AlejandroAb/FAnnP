#!/bin/bash\
#
# This script uses a metadata file.
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file (metadata custom information): config["meta_data"]["bin_meta_data"]
echo $1
# param: Run ID
echo $2
# output file: {run}/sql/metadata.tsv
echo $3

cat $1 | awk -v run_name=$2 -F"\t" 'BEGIN{OFS="\t"}
 {
 gsub("\r","",$0);
 if(NR==1){
   gsub("#","",$0);
   split($0,column_names,"\t");
   }
 else{
   split($0,row,"\t");
   for(i=2;i<length(column_names)+1;i++){
     if(length(row[i])==0){
       row[i]="\\N";
       }
	 print "0", row[1], run_name, column_names[i], row[i];
     }
   }
 }' > $3
