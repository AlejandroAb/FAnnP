#!/bin/bash\
#
# Part of the FAnnP pipeline to map '_sql.out' files (processed hmmr output file).
# Takes all entries from the second, third, and fourth column (domains, evals, and scores).
# Cycles through these domain entries to map them against database files.
# Once the entries are mapped, the program puts the domains together again for each gene.
#
# tigr
#

cat $2 | awk -F"\t" 'BEGIN{OFS="\t"}
 FNR==NR{h[$1]=$2;a[$1]=$3;b[$1]=$4;c[$1]=$5;next}
 found_id"";
 found_name"";
 found_desc="";
 found_ec="";
 found_gene="";
 found_evals="";
 found_scores="";
 {n=split($2,domains,",");
 n=split($3,evals,",");
 n=split($4,scores,",");
  for(i=1;i<=n;i++){
      found_id=found_id domains[i]",";
	  found_name=found_name h[domains[i]]",";
      found_desc=found_desc a[domains[i]]",";
      found_ec=found_ec b[domains[i]]",";
	  found_gene=found_gene c[domains[i]]",";
      found_evals=found_evals evals[i]",";
      found_scores=found_scores scores[i]",";
      }
    if(h[domains[1]]){
      print $1,
	    substr( found_id, 1, length(found_id)-1 ),
	    substr( found_name, 1, length(found_name)-1 ),
        substr( found_desc, 1, length(found_desc)-1 ),
        substr( found_ec, 1, length(found_ec)-1 ),
		substr( found_gene, 1, length(found_gene)-1 ),
        substr( found_evals, 1, length(found_evals)-1 ),
        substr( found_scores, 1, length(found_scores)-1 );
      }
    else {
      print $1,$2,"-","-","-","-","-","-";
    }
	found_id="";
	found_name="";
    found_desc="";
    found_ec="";
	found_gene="";
    found_evals="";
    found_scores="";
    }' - $1 > $3

