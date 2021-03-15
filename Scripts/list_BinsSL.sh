#!/bin/bash
#needs 3 arguments
#$1 Bin directory;
#$2 Bin extension
#$3 Configfile
#$4 Optional new configfile

if  [ $# -lt 3 ]; then
    echo -e "This script needs 3 arguments\n"
    echo -e "Bin directory\n"
    echo -e "Extension of the bins"
    echo -e "Config file"
   exit 1
fi

#bins=$(ls -l $1/*.$2 | awk 'BEGIN{ORS=","} NR>1{print "\""$NF"\""}' | sed "s/,$/]/"
bins=$(ls -l $1/*.$2 | sed "  s/ ->.*$// ;  s/\.$2$//" | awk -F'/' '{print $NF}' | awk 'BEGIN{ORS=","} {print "\""$NF"\""}' | sed "s/,$/]/" | sed 's/\"/\\\"/g ; s/\//\\\//g')

if [ -z "$4" ]; then
   sed  -i "s/bin_list:.*/bin_list: [${bins}/" $3 
else
   sed  "s/bin_list:.*/bin_list: [${bins}/" $3 > $4
fi
exit 0

