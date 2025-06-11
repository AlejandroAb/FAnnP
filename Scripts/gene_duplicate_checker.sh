#!/bin/bash\
#
# gene file
# The script checks if the names for every gene are not present twice.
# If there is a gene name that appears more than one time the fannp will crash (intentionally)!!!
#

# input file (all proteins): {run}/prokka/All_bins.faa
echo $1
# param (T or F): config["rename_bins"]
echo $2
# output file: {run}/prokka/duplications_checked_log.txt
echo $3

if [ -z $(cat $1 | egrep '^>' | awk '{print $1}' | sort | uniq -c | egrep -v '     1 ') ]
then {
	echo "Success!" > $3
} else {
	if [ $2 = "T" ]
	then {
		echo "Error, there are duplicate gene names!!! Renaming is already on so there must be a duplicate bin and gene name combination!"
		echo "Error, there are duplicate gene names!!! Renaming is already on so there must be a duplicate bin and gene name combination!" > $3
	} else {
		echo "Error, there are duplicate gene names!!! Turning rename_bins: 'T' in config might resolve this since there could be a duplicate bin or gene name!"
		echo "Error, there are duplicate gene names!!! Turning rename_bins: 'T' in config might resolve this since there could be a duplicate bin or gene name!" > $3
	}
	fi
	exit 1
}
fi

