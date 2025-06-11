#!/bin/bash\
#
# This script takes a list of bins and grabs all contig information from the .fasta files used for the MySQL database. 
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file (pipeline input bins): expand(config["bin_dir"]+"{bin}."+config["bin_ext"], bin=config["bin_list"])
echo $1
# param (file extension type of bins): config["bin_ext"]
echo $2
# param (T or F): config["rename_bins"]
echo $3
# output file: {run}/sql/contig.tsv
echo $4
# Meta-Cascabel output file (contig coverage): config["meta_cascabel"]["contig_coverage"]
echo $5

for file in $1
do {
	echo $file
	if [ -z $5 ]; then {
		cat $file | egrep -v "^#" | awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' | awk -v file_loc=$file -v bin_ext=$2 -v renamed=$3 'BEGIN{OFS="\t"}
		 NR%2{gsub("\r","",$0);a=$0;next}
		 {
		 tot=length($1);
		 gsub("(A|T)","",$1);
		 cg=length($1);
		 sub(">","",a);
		 split(file_loc,c,"/");
		 split(a,q,"( |\\t)");
		 sub("\\."bin_ext,"",c[length(c)]);
		 if(renamed=="T"){
	       print c[length(c)]"|"q[1], c[length(c)], "\\N", tot, cg/tot*100, "\\N";
	       }else{
	       print q[1], c[length(c)], "\\N", tot, cg/tot*100, "\\N";
	       }
		 }'
	} else {
		cat $5 | awk -F"\t" 'BEGIN{OFS="\t"}
		 FNR==NR{gsub("\r","",$0);h[$1]=$2;next}
		 {
		 print $1, $2, "\\N", $3, $4, h[$5];
		 }' - <(cat $file | egrep -v "^#" | awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' | awk -v file_loc=$file -v bin_ext=$2 -v renamed=$3 'BEGIN{OFS="\t"}
			 NR%2{gsub("\r","",$0);a=$0;next}
			 {
			 tot=length($1);
			 gsub("(A|T)","",$1);
			 cg=length($1);
			 sub(">","",a);
			 split(file_loc,c,"/");
			 split(a,q,"( |\\t)");
			 sub("\\."bin_ext,"",c[length(c)]);
			 if(renamed=="T"){
	           print c[length(c)]"|"q[1], c[length(c)], tot, cg/tot*100, q[1];
	           }else{
	           print q[1], c[length(c)], tot, cg/tot*100, q[1];
	           }
			 }')
	}
	fi
}
done | egrep -v "^/" > $4




