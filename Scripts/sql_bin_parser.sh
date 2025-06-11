#!/bin/bash\
#
# This script takes a list of bins and grabs all bin information from the .fasta files used for the MySQL database. 
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
# Meta-Cascabel output file (final bins): config["meta_cascabel"]["location_final_bins"]
echo $5
# Meta-Cascabel output file (contig coverage): config["meta_cascabel"]["contig_coverage"]
echo $6

for file in $1
do {
	echo $file
	if [ -z $5 ] || [ -z $6 ]; then {
		cat $file | egrep -v "^#" | awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' | awk -v file_loc=$file -v bin_ext=$2 -v run_name=$3 'BEGIN{OFS="\t"}
		 NR%2{gsub("\r","",$0);a=$0;next}
		 {
		 tot+=length($1);
		 gsub("(A|T)","",$1);
		 cg+=length($1);
		 sub(">","",a);
		 split(file_loc,c,"/");
		 split(a,q,"( |\\t)");
		 sub("\\."bin_ext,"",c[length(c)]);
		 } END {
		 print c[length(c)], run_name, file_loc, bin_ext, "\\N", tot, cg/tot*100, "\\N";
		 }'
	} else {
		cat $5 | awk -F"\t" 'BEGIN{OFS="\t"}
		 ARGIND == 1 {gsub("\r","",$0);a[$1]=$17;next}
		 ARGIND == 2 {gsub("\r","",$0);
		 split($1,c,"\\.");
		 gsub("\\."c[length(c)], "", $1);
		 b[$1] += $2;
		 d[$1]++;next}
		 END {
		 print $1, $2, $3, $4, a[$1], $5, $6, b[$1]/d[$1];
		 }' - $6 - <(cat $file | egrep -v "^#" | awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' | awk -v file_loc=$file -v bin_ext=$2 -v run_name=$3 'BEGIN{OFS="\t"}
			 NR%2{gsub("\r","",$0);a=$0;next}
			 {
			 tot+=length($1);
			 gsub("(A|T)","",$1);
			 cg+=length($1);
			 sub(">","",a);
			 split(file_loc,c,"/");
			 split(a,q,"( |\\t)");
			 sub("\\."bin_ext,"",c[length(c)]);
			 } END {
			 print c[length(c)], run_name, file_loc, bin_ext, tot, cg/tot*100;
			 }')
	}
	fi
}
done | egrep -v "^/" > $4
