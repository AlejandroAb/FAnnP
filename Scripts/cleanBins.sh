#!/bin/bash
#needs 5 arguments
#$1 Bin directory;
#$2 Bin extension: {run}/contig_mapping/B_GenomeInfo.txt
#$3 Coverage/GC rejected
#$4 Taxonomy rejected
#$5 Output dir

if  [ $# -lt 5 ]; then
    echo -e "This script needs 5 arguments\n"
    exit 1
fi
#init log file with header
echo -e "BinID\tNumContigs\tSurvivingContigs\n" > $5/clean.log
#concat uniq contigs 
cat $3 $4 |  cut -f1,2 | sort | uniq > $5/ref.ids
for file in $1/*.$2; do
 #retain only the bin name
 bin=$(echo $file | awk -F"/" '{print $NF}' | sed 's/\.$2$//')
 #retain only contigs belonging to the target bin in case more than one bin have the same contig name 
 #this is not possible if the bin comes from the metagenomic pipeline but possible if it comes from some other source
 cat $5/ref.ids | awk -F"\t" -v bid=${bin} '$1 ~ bid{print $2}' | grep -F -w -v -f - $file | grep "^>" | sed 's/>//' > $5/${bin}.ok
 #subseq only passing contigs
 seqtk subseq $file $5/${bin}.ok > $5/${bin}.$2
 oric=$(grep -c  "^>" $file)
 newc=$(grep -c  "^>" $5/${bin}.$2)
 echo -e $bin"\t"$oric\t$newc >> $5/clean.log 
 rm $5/${bin}.ok
done

rm $5/ref.ids

exit 0 

