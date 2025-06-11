#!/bin/bash\
#
# Part of the FAnnP pipeline to map '_sql.out' files (processed hmmr output file).
# Takes all entries from the second, third, and fourth column (domains, evals, and scores).
# Cycles through these domain entries to map them against database files.
# Once the entries are mapped, the program puts the domains together again for each gene.
#
# ko
# Also checks if the domain score is higher than the mapped bitscore,
# in which case the domain entry will get marked by the word high_score in the last column. 
#

cat $2 | awk -F"\t" 'BEGIN{OFS="\t"}
 FNR==NR{h[$1]=$2;a[$1]=$12;next}
 found_id"";
 found_evals="";
 found_scores="";
 found_bitscore"";
 found_desc="";
 found_confidence="";
 {n=split($2,domains,",");
 n=split($3,evals,",");
 n=split($4,scores,",");
  for(i=1;i<=n;i++){
      found_id=found_id domains[i]",";
      found_evals=found_evals evals[i]",";
      found_scores=found_scores scores[i]",";
	  found_bitscore=found_bitscore h[domains[i]]",";
      found_desc=found_desc a[domains[i]]",";
	  if (scores[i] > h[domains[i]]) {
	    found_confidence=found_confidence "high_score,";
		}
	  else {
	    found_confidence=found_confidence "-,";
	    }
      }
    if(h[domains[1]]){
      print $1,
	    substr( found_id, 1, length(found_id)-1 ),
        substr( found_evals, 1, length(found_evals)-1 ),
        substr( found_scores, 1, length(found_scores)-1 ),
		substr( found_bitscore, 1, length(found_bitscore)-1 ),
        substr( found_desc, 1, length(found_desc)-1 ),
		substr( found_confidence, 1, length(found_confidence)-1 );
      }
    else {
      print $1,$2,"-","-","-","-","-";
    }
	found_id="";
    found_evals="";
    found_scores="";
	found_bitscore="";
    found_desc="";
    }' - $1 > $3

