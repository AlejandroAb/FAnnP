#!/bin/bash\
#
# Creates a Marimo Notebook for the specific run being executed.
#

# input file: {run}/sql/run_to_sql.txt
echo $1
# param: Run ID
echo $2
# script file: Scripts/MySQL_FAnnP_Analyzer.py
echo $3
# param: MySQL username
echo $4
# param: MySQL password
echo $5
# param: MySQL database
echo $6
# output file: {run}/sql/run_to_notebook.txt
echo $7

user_name=$(whoami)
cat $3 | sed "s/User_Name_Placeholder/${user_name}/g" | sed "s/Run_ID_Placeholder/${2}/g" | sed "s/MySQL_Name_Placeholder/${4}/g" | sed "s/MySQL_Password_Placeholder/${5}/g" | sed "s/MySQL_Database_Placeholder/${6}/g" > "MySQL_FAnnP_Analyzer_${2}.py"
echo "Creating Marimo Notebook: MySQL_FAnnP_Analyzer_${2}.py" > $7
