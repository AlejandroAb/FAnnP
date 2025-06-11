#!/bin/bash\
#
# This script takes a list of bins and grabs all gene information from the .faa proteins used for the MySQL database. 
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file (protein files): expand("{run}/prokka/{bin}/{bin}.faa", bin=config["bin_list"], run=run)
echo $1
# param (file extension type of bins): config["bin_ext"]
echo $2
# param (T or F): config["rename_bins"]
echo $3
# output file: {run}/sql/feature.tsv
echo $4

for file in $1
do {
	echo $file
	cat $file | egrep -v "^#" | awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' | awk -v file_loc=$file -v bin_ext=$2 -v renamed=$3 'BEGIN{OFS="\t"}
	 NR%2{gsub("\r","",$0);a=$0;next}
	 {
	 tot=length($1)*3;
	 all_tot+=tot;
	 gsub("(A|T)","",$1);
	 cg=length($1);
	 sub(">","",a);
	 split(file_loc,c,"/");
	 split(a,q,"( |\\t)");
	 sub("\\."bin_ext,"",c[length(c)]);
	 if(renamed=="T"){
       print c[length(c)]"|"q[1], c[length(c)], c[length(c)], "\\N", all_tot-tot+1, all_tot, "\\N", "\\N", "\\N";
	   }else{
	   print q[1], c[length(c)], c[length(c)], "\\N", all_tot-tot+1, all_tot, "\\N", "\\N", "\\N";
	   }
	 }'
}
done | egrep -v "^/" > $4




