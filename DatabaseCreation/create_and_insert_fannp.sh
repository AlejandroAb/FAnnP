#!/bin/bash\
#
# Populate MySQL tables with database information using create_fannp_database.py
# For COG, Pfam, and ...
#

#
#
# Settings

# MySQL account information
# MySQL Username
user=$1
# MySQL Password
password=$2

# Database name, if we want to create a new one we can change this from fannp to something else
database=$3
# Create database or not (turn to False if already exists and otherwise True)
create_database=$4
# Load MySQL map (from which map you CAN load MySQL files)
mysql_loadmap=$5
# Extra params, fill it in the extra params:
# --reload_databases	= Reload databases.
# --skip_reload_kegg	= Skip reloading kegg since reloading KEGG will around 30 mins more time than Pfam and COG.
extra_params="${6} ${7}"


# DB files
declare cog_info="DatabaseFiles/cog-24.def.tab"
declare cog_family_info="DatabaseFiles/CATEGORIES.txt"

declare pfam_info="DatabaseFiles/Pfam-A.clans.tsv"
declare pfam_type_info="DatabaseFiles/Pfam-A.hmm.dat"

# TSV files
declare tsv_files_location=$(printf $mysql_loadmap"static_"$database"/")
mkdir $tsv_files_location

#
#
# Static named files
declare cog_family_sql=$(printf $tsv_files_location"cog_family.tsv")
declare cog_sql=$(printf $tsv_files_location"cog.tsv")
declare cog_has_cog_family_sql=$(printf $tsv_files_location"cog_has_cog_family.tsv")

declare pfam_clan_sql=$(printf $tsv_files_location"pfam_clan.tsv")
declare pfam_sql=$(printf $tsv_files_location"pfam.tsv")

declare ko_sql=$(printf $tsv_files_location"ko.tsv")
declare module_sql=$(printf $tsv_files_location"module.tsv")
declare pathway_sql=$(printf $tsv_files_location"pathway.tsv")
declare ec_sql=$(printf $tsv_files_location"ec.tsv")
declare module_has_ko_sql=$(printf $tsv_files_location"module_has_ko.tsv")
declare pathway_has_ko_sql=$(printf $tsv_files_location"pathway_has_ko.tsv")
declare pathway_has_module_sql=$(printf $tsv_files_location"pathway_has_module.tsv")
declare ko_has_ec_sql=$(printf $tsv_files_location"ko_has_ec.tsv")

#
#
# Create Database (We will need to create this before using the python connector!)
set +H
if [ ${create_database} = "True" ]
then {
  declare create_database_file=$(printf $tsv_files_location"create_database.sql")
	echo "CREATE DATABASE ${database};" > $create_database_file
  echo "CREATE USER IF NOT EXISTS 'faanp_user'@'%' IDENTIFIED BY 'FaanP2025!';" >> $create_database_file
	echo "GRANT SELECT ON ${database}.* TO   'faanp_user'@'%';" >> $create_database_file
	echo "GRANT UPDATE ON ${database}.* TO   'faanp_user'@'%';" >> $create_database_file
	echo "GRANT INSERT ON ${database}.* TO   'faanp_user'@'%';" >> $create_database_file
	echo "GRANT ALL PRIVILEGES ON ${database}.* TO   'faanp_user'@'%';" >> $create_database_file
	echo "UPDATE mysql.user  SET file_priv='Y'   WHERE user='faanp_user';" >> $create_database_file
	mysql -u $user -p$password < $create_database_file  # Might have to use 'set +H' if '!' symbol is an issue.
}
fi

#
#
# Generate TSV tables

#
# COG
# Parse information for cog_family
cat $cog_family_info | awk -F"\t" 'BEGIN{OFS="\t"}
 {
 gsub("\r","",$0);
 print $1, $2;
 }' > $cog_family_sql

# Parse information for cog
cat $cog_info | iconv -c -f utf-8 -t ascii | awk -F"\t" 'BEGIN{OFS="\t"}
 {
 gsub("\r","",$0);
 gsub("'\''","\\'\''",$3);
 print $1, $3;
 }' > $cog_sql

# Parse information for cog_has_cog_family
cat $cog_info | awk -F"\t" 'BEGIN{OFS="\t"}
 {
 gsub("\r","",$0);
 for(i=1; i<=length($2);i++){
  print substr($2, i, 1), $1;
  }
 }' > $cog_has_cog_family_sql

#
# Pfam
# Parse information for pfam_clan
cat $pfam_info | awk -F"\t" 'BEGIN{OFS="\t"}{print $2,$3}' | sort | uniq | egrep "^C" | awk -F"\t" 'BEGIN{OFS="\t"}
 {
 gsub("\r","",$0);
 print $1, $2;
 }' > $pfam_clan_sql
 
# Parse information for pfam
cat $pfam_info | iconv -c -f utf-8 -t ascii | awk -F"\t" 'BEGIN{OFS="\t"}
 FNR==NR{gsub("\r","",$0);h[$1]=$4;i[$1]=$5;j[$1]=$2;next}
 {
 gsub("\r","",$0);
 gsub("'\''","\\'\''",i[$1]);
 print $1, h[$1], i[$1], $2, j[$1];
 }' - <(cat $pfam_type_info | awk 'BEGIN{OFS="\t"}
 {
 if($2=="AC")
   {
   split($3,a,".");
   }
 if($2=="TP")
   {
   a[2]=$3;
   }
 if($1=="//")
   {
   print a[1],a[2];
   }
 }' | sort) > $pfam_sql

#
# KO
# Parse information for ko
reload=${7:2}
if [ "${reload}" = "skip_reload_kegg" ]
then {
    echo "Skipping KEGG loading."
} else {
	echo "Loading KEGG tables."
	python /export/lv10/user/ttensen/work_dir/Database_scripts/create_kegg_tables.py -l $tsv_files_location
}
fi

# Insert data into database
echo "Creating Database."
python /export/lv10/user/ttensen/work_dir/Database_scripts/create_fannp_database.py -u $user -p $password -d $database -l $tsv_files_location $extra_params


# Remove temp tsv files
echo "Removing files."
rm $cog_family_sql $cog_sql $cog_has_cog_family_sql
rm $pfam_clan_sql $pfam_sql
if [ "${reload}" = "skip_reload_kegg" ]
then {
    echo "No KEGG files to remove."
} else {
	rm $ko_sql $module_sql $pathway_sql $ec_sql $module_has_ko_sql $pathway_has_ko_sql $pathway_has_module_sql $ko_has_ec_sql
}
fi
# Remove SQL file
if [ ${create_database} = "True" ]
then {
	rm $create_database_file
}
fi
