#!/bin/bash
#needs 4 arguments
#$1 Full path to bin
#$2 Bin ext
#$3 Bin Coverage/GC rejected
#$4 Taxonomy rejected
#$5 Output dir

if  [ $# -lt 4 ]; then
    echo -e "This script needs 4 arguments\n"
    exit 1
fi
#init log file with header
#echo -e "BinID\tNumContigs\tSurvivingContigs" > $5/clean.log
#concat uniq contigs 
 #retain only the bin name
 bin=$(echo $1 | awk -F"/" '{print $NF}' | sed "s/\.$2$//")
 #retain only contigs belonging to the target bin in case more than one bin have the same contig name 
 #this is not possible if the bin comes from the metagenomic pipeline but possible if it comes from some other source
 cat  $3 $4 |  cut -f1,2 | sort | uniq | awk -F"\t" -v bid=${bin} '$1 ~ bid{print $2}' | grep -F -w -v -f - $1 | grep "^>" | sed 's/>//' > $5/${bin}.ok
 #subseq only passing contigs
 seqtk subseq $1 $5/${bin}.ok > $5/${bin}.$2
 oric=$(grep -c  "^>" $1)
 newc=$(grep -c  "^>" $5/${bin}.$2)
 echo -e $bin"\t"$oric"\t"$newc > $5/clean.${bin}.log 
 rm $5/${bin}.ok

exit 0 

