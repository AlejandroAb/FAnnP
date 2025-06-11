#!/bin/bash\
#
# Part of the FAnnP pipeline to map '_sql.out' files (processed hmmr output file).
# Takes all entries from the second, third, and fourth column (domains, evals, and scores).
# Cycles through these domain entries to map them against database files.
# Once the entries are mapped, the program puts the domains together again for each gene.
#
# pfam
#

cat $2 | awk -F"\t" 'BEGIN{OFS="\t"}
 FNR==NR{h[$4]=$1;a[$4]=$5;next}
 found_pfam_ids"";
 found_domains="";
 found_association="";
 found_evals="";
 found_scores="";
 {n=split($2,domains,",");
 n=split($3,evals,",");
 n=split($4,scores,",");
  for(i=1;i<=n;i++){
      found_pfam_ids=found_pfam_ids h[domains[i]]",";
      found_domains=found_domains domains[i]",";
      found_association=found_association a[domains[i]]",";
      found_evals=found_evals evals[i]",";
      found_scores=found_scores scores[i]",";
      }
    if(h[domains[1]]){
      print $1,
	    substr( found_pfam_ids, 1, length(found_pfam_ids)-1 ),
        substr( found_domains, 1, length(found_domains)-1 ),
        substr( found_association, 1, length(found_association)-1 ),
        substr( found_evals, 1, length(found_evals)-1 ),
        substr( found_scores, 1, length(found_scores)-1 );
      }
    else {
      print $1,"-","-","-","-","-";
    }
	found_pfam_ids="";
    found_domains="";
    found_association="";
    found_evals="";
    found_scores="";
    }' - $1 > $3

