#!/bin/bash\
#
# gene file
# The script checks if the names for every gene is 100 characters or lower for the MySQL database.
# If there is a gene name that has more than 100 character the fannp will crash (intentionally)!!!
#

# input file (protein files): expand("{run}/prokka/renamed/{bin}.faa", bin=config["bin_list"], run=run)
echo $1
# output file: {run}/prokka/renamed_checked_log.txt
echo $2

for file in $1
do {
	if [ -z $(cat $1 | awk -v a="False" 'BEGIN{OFS="\t"}{gsub("\r","",$0);if(substr($1,1,1) == ">"){if(length($1) > 100 && a == "False"){a="True";print $1;}}}') ]
	then {
		echo "Success!" > $2
	} else {
		echo "Error, there are gene names with a name length of > 100!!! This cannot be possible for the MySQL database!"
		echo "Error, there are gene names with a name length of > 100!!! This cannot be possible for the MySQL database!" > $2
		exit 1
	}
	fi
}
done
