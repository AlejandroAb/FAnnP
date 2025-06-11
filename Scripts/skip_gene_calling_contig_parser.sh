#!/bin/bash\
#
# Uses the faa. format to format it into the Contig_Protein_list.txt format Prokka and Prodigal creates. 
#

# input file (protein files): expand("{run}/prokka/renamed/{bin}.faa", bin=config["bin_list"], run=run)
echo $1
# output file (contig/protein list): {run}/contig_mapping/Contig_Protein_list.txt
echo $2

for target_faa_file in ${1} 
do {
	target_faa_name=$(printf "${target_faa_file}" | awk -F"/" '{print $NF}' | awk -F"." '{NF--; print}')
	contig_length=$(($(cat $target_faa_file | egrep -v "^>" | tr -d "\n" | wc -c)*3))
	cat $target_faa_file | awk 'NR==0 {print ; next} {printf /^>/ ? "\n"$0"\t" : $1} END {printf "\n"}' | tail -n+2 | awk -v target_faa_str=$target_faa_name -v contig_int=$contig_length -F"\t" 'BEGIN{OFS="\t"};{sub(">", "", $1);print "1", target_faa_str, target_faa_str, contig_int, $1, "NA"}'
}
done > $2
