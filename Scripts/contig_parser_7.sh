#!/bin/bash\
#
# Combines 7 rules from the old FAnnP version.
# Steps clean_concatenated_bins to add_header_to_file are all combined in this script.
#

# input file (proteins): {run}/prokka/All_bins.faa
echo $1
# input file (bin/contig/protein list): {run}/contig_mapping/Bin_Contig_Protein_list_merged.txt
echo $2
# output file: {run}/contig_mapping/B_GenomeInfo.txt
echo $3
# Meta-Cascabel param (T or F): config["meta_cascabel"]["run"]
echo $4
# Meta-Cascabel output file (final bins): config["meta_cascabel"]["location_final_bins"]
echo $5
# Meta-Cascabel output file (contig coverage): config["meta_cascabel"]["contig_coverage"]
echo $6

# Add headers
printf "accession\tBinID\tTaxString_Bin\tBin_avg_cov\tBin_avg_gc\tNewContigID\tOldContigId\tContigIdMerge\tContigNewLength\tContigOldLength\tContig_GC\tContig_cov\tProteinID\tProteinLength\tProkka\n" > $3

# Look through all files with ARGIND for mapping
if [ $4 = "T" -a -n $5 -a -n $6 ]
then {
	cat $1 | awk 'NR==0 {print ; next} {printf /^>/ ? "\n"$1"\t" : $1} END {printf "\n"}' | tail -n +2 | sed 's/>//g' | awk -F"\t" 'BEGIN{OFS="\t"}
		 ARGIND == 1 {gsub("\r","",$0);a[$1]=$2;next}
		 ARGIND == 2 {gsub("\r","",$0);b[$1]=$17;c[$1]=$23;d[$1]=$24;if(length($17)==0){b[$1]="NA"};if(length($23)==0){c[$1]="NA"};if(length($24)==0){d[$1]="NA"};next}
		 ARGIND == 3 {gsub("\r","",$0);e[$1]=$2;if(length($2)==0){e[$1]="NA"};next}
		 ARGIND == 4 {gsub("\r","",$0);
		 print $1,$2,b[$2],c[$2],d[$2],$3,$4,$5,$6,$7,$8*100,e[$4],$9,length(a[$1]),$10;
		 }' - $5 $6 $2 >> $3
} else {
	cat $1 | awk 'NR==0 {print ; next} {printf /^>/ ? "\n"$1"\t" : $1} END {printf "\n"}' | tail -n +2 | sed 's/>//g' | awk -F"\t" 'BEGIN{OFS="\t"}
		 ARGIND == 1 {gsub("\r","",$0);a[$1]=$2;next}
		 ARGIND == 2 {gsub("\r","",$0);
		 print $1,$2,"NA","NA","NA",$3,$4,$5,$6,$7,$8*100,"NA",$9,length(a[$1]),$10;
		 }' - $2 >> $3
}
fi


