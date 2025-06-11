#!/bin/bash\
#
# Part of the FAnnP pipeline to map '_sql.out' files (processed hmmr output file).
# Takes all entries from the second, third, and fourth column (domains, evals, and scores).
# Cycles through these domain entries to map them against database files.
# Once the entries are mapped, the program puts the domains together again for each gene.
#
# cog
# A pathway might not always be found for cog, so if checks if the column is empty,
# in which case a "-" will be printed in its place.
#

cat $2 | awk -F"\t" 'BEGIN{OFS="\t"}
 FNR==NR{h[$1]=$3;a[$1]=$2;b[$1]=$5;next}
 found_id"";
 found_desc"";
 found_pathway_id="";
 found_pathway="";
 found_evals="";
 {n=split($2,domains,",");
 n=split($3,evals,",");
  for(i=1;i<=n;i++){
      found_id=found_id domains[i]",";
	  found_desc=found_desc h[domains[i]]",";
      found_pathway_id=found_pathway_id a[domains[i]]",";
	  if (b[domains[i]]){
        found_pathway=found_pathway b[domains[i]]",";
	    }
	  else {
	    found_pathway=found_pathway "-,";
	    }
      found_evals=found_evals evals[i]",";
      }
    if(h[domains[1]]){
      print $1,
	    substr( found_id, 1, length(found_id)-1 ),
	    substr( found_desc, 1, length(found_desc)-1 ),
        substr( found_pathway_id, 1, length(found_pathway_id)-1 ),
        substr( found_pathway, 1, length(found_pathway)-1 ),
        substr( found_evals, 1, length(found_evals)-1 );
      }
    else {
      print $1,$2,"-","-","-","-";
      }
	found_id="";
	found_desc="";
    found_pathway_id="";
    found_pathway="";
    found_evals="";
    }' - $1 > $3

