#!/bin/bash\
#
# Takes important data specific to the FAnnP run.
# The script constructs a tsv table for inserting it into MySQL tables.
#

# param: Run ID
echo $1
# param: Run Name
echo $2
# param: Current Date Time
echo $3
# output file: {run}/sql/run.tsv
echo $4

run_name=$(echo $2 | awk -F"/" '{print $NF}')
work_dir=$(pwd)
version=$(cat Snakefile | egrep "^Version: " | awk '{print $2}')
user=$(whoami)
if [[ $2 == /* ]]; then output_dir=$2; else output_dir=$(printf "${work_dir}/${2}"); fi

echo "" | awk -v run_id=$1 -v run_name=$run_name -v start_time="$3" -v work_dir=$work_dir -v version=$version -v user=$user -v output_dir=$output_dir 'BEGIN{OFS="\t"}
 {
 print run_id, run_name, start_time, work_dir, version, user, output_dir;
 }' > $4
