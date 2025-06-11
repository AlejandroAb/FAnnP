#!/bin/bash\
#
# This script takes a list of bins and grabs all bin information from the .faa proteins used for the MySQL database. 
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file (pipeline input bins): expand(config["bin_dir"]+"{bin}."+config["bin_ext"], bin=config["bin_list"])
echo $1
# param (file extension type of bins): config["bin_ext"]
echo $2
# param: Run ID
echo $3
# output file: {run}/sql/bin.tsv
echo $4

for file in $1
do {
	echo $file
	cat $file | egrep -v "^#" | awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' | awk -v file_loc=$file -v bin_ext=$2 -v run_name=$3 'BEGIN{OFS="\t"}
	 NR%2{gsub("\r","",$0);a=$0;next}
	 {
	 tot+=length($1);
	 sub(">","",a);
	 split(file_loc,c,"/");
	 split(a,q,"( |\\t)");
	 sub("\\."bin_ext,"",c[length(c)]);
	 } END {
	 print c[length(c)], run_name, file_loc, bin_ext, "\\N", tot*3, "\\N", "\\N";
	 }'
}
done | egrep -v "^/" > $4
