#!/bin/bash
#needs 3 arguments
#$1 Run directory;
#$2 Base file: {run}/contig_mapping/B_GenomeInfo.txt
#$3 Outputile

if  [ $# -lt 3 ]; then
    echo -e "This script needs 3 arguments\n"
    echo -e "Run directory\n"
    echo -e "Base file. File with Genome information. {run}/contig_mapping/B_GenomeInfo.txt"
    echo -e "Outputfile"
   exit 1
fi

cp $2 $3
for file in $1/*/*_map.tsv; do
  cat $file | awk -F "\t" 'NR==FNR{if(NR==1){header=$2; for(i=3;i<=NF;i++){header=header"\t"$i} ; br="-"; \
              for(i=2;i<=NF;i++){br=br"\t-"}} else{row=$2; for(i=3;i<=NF;i++){row=row"\t"$i}; r[$1]=row};next} \
              {if(FNR==1){print $0"\t"header}else{if(r[$1]){print $0"\t"r[$1]}else{print $0"\t"br}}}' \
   - $3 > $3.tmp
 mv $3.tmp  $3
done

exit 0

