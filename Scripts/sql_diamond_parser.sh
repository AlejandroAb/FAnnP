#!/bin/bash\
#
# NCBI diamond hits
# The script constructs a tsv table for inserting it into MySQL tables.
#

# input file: {run}/diamond/All_bins.tsv
echo $1
# database file: config["diamond"]["taxid_to_taxonomy"]
echo $2
# output file: {run}/sql/diamond.tsv
echo $3

cat $2 | awk -v hit="False" -v counter=0 -F"\t" 'BEGIN{OFS="\t"}
 FNR==NR{gsub("\r","",$0);h[$1]=$2;next}
 {
 gsub("\r","",$0)
 id=$1;
 if(id==old_id){
   for(locs in store_locs){
      split(store_locs[locs],loc_list,",");
	  min=int((loc_list[3]<$13?loc_list[3]:$13)/4)
      if(loc_list[1]+min<$8 && loc_list[2]-min>$7){
	    hit="False";
		}
      }
   }
 else{
   old_id=id;
   hit="True";
   counter=1;
   delete store_locs;
   }
 if(hit=="True"){
   found_tax_ids="";
   found_tax_descs="";
   n=split($15,tax_ids,";");
   for(i=1;i<=n;i++){
      found_tax_ids=found_tax_ids tax_ids[i]",";
      found_tax_descs=found_tax_descs h[tax_ids[i]]",";
      }
   print "0", $1, $4, $5, substr( found_tax_ids, 1, length(found_tax_ids)-1 ), substr( found_tax_descs, 1, length(found_tax_descs)-1 ), $7, $8, $9, $10, $11, $12, $14;
   store_locs[counter]=$7","$8","$13;
   counter++;
   hit="True";
   }
 }' - $1 > $3
 
