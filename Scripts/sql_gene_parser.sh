#!/bin/bash\
#
# This script takes a list of bins and grabs all gene information from the .gff files used for the MySQL database. 
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file (gff files): expand("{run}/prokka/{bin}/{bin}.gff", bin=config["bin_list"], run=run)
echo $1
# param (T or F): config["rename_bins"]
echo $2
# output file: {run}/sql/feature.tsv
echo $3

for file in $1
do {
	echo $file
    cat $file | egrep -v "^#" | awk -v file_loc=$file -v renamed=$2 -F"\t" 'BEGIN{OFS="\t"}
	 {
	 split($9,a,";");
	 split(a[1],b,"_");
	 split(file_loc,c,"/");
	 sub("\\.gff","",c[length(c)]);
	 if(renamed=="T"){
	   print c[length(c)]"|"$1"_"b[2], c[length(c)]"|"$1, c[length(c)], "\\N", $4, $5, $7, $3, $9;
	   }else{
	   print $1"_"b[2], $1, c[length(c)], "\\N", $4, $5, $7, $3, $9;
	   }
	 }'
}
done | egrep -v "^/" > $3




