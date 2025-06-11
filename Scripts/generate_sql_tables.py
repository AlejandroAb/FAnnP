#!/usr/bin/python
# Tijn 2024-2025 :)
# Last recorded update: 03-2025
#
# Inserts all the Run specific MySQL tables and inserts all the information tables from the pipeline (the created sql TSV files).
# Script runs by using the Functional Annotation Pipeline.
#

# Libraries
import sys
import os
import argparse
import shutil
import yaml
import mysql.connector
from mysql.connector import errorcode


# Arguments you can use when running this script.
parser = argparse.ArgumentParser(description='Inserts tsv files into MySQL database')
parser.add_argument('-r','--run_name', metavar='str', required = True, type=str, help='the id of the run')
parser.add_argument('-c','--config_name', metavar='file', required = True, type=str, help='config file used by the snakemake pipeline')
parser.add_argument('-o','--output', metavar='file', required = True, type=str, help='txt file')
parser.add_argument('--store_metadata', action="store_true", help='to declare metadata from the bin will be stored')
parser.add_argument('--store_ko', action="store_true", help='to declare ko information will be stored')
parser.add_argument('--store_cog', action="store_true", help='to declare cog information will be stored')
parser.add_argument('--store_arCOG', action="store_true", help='to declare atCOG information will be stored')
parser.add_argument('--store_pfam', action="store_true", help='to declare pfam information will be stored')
parser.add_argument('--store_tigr', action="store_true", help='to declare tigr information will be stored')
parser.add_argument('--store_cazy', action="store_true", help='to declare cazy information will be stored')
parser.add_argument('--store_signalP', action="store_true", help='to declare signalP information will be stored')
parser.add_argument('--store_merops', action="store_true", help='to declare merops information will be stored')
parser.add_argument('--store_transporter', action="store_true", help='to declare transporterDB information will be stored')
parser.add_argument('--store_hydDB', action="store_true", help='to declare hydDB information will be stored')
parser.add_argument('--store_hmmr1', action="store_true", help='to declare custom hmmr database 1 information will be stored')
parser.add_argument('--store_hmmr2', action="store_true", help='to declare custom hmmr database 2 information will be stored')
parser.add_argument('--store_blast1', action="store_true", help='to declare custom blast database 1 information will be stored')
parser.add_argument('--store_blast2', action="store_true", help='to declare custom blast database 2 information will be stored')
args = parser.parse_args()


# Check if the arguments required are supplied by the user (will do it automatically in the pipeline).
if args.run_name is None:
    print("Run name is missing!")
    sys.exit(-1)
if args.config_name is None and not os.path.exists(args.config_name):
    print("Config name is missing!")
    sys.exit(-1)
if args.output is None and not os.path.exists(args.output):
    print("Output file is missing!")
    sys.exit(-1)

# Open target config file
with open(args.config_name, "r") as config_file:
    config = yaml.safe_load(config_file)

# Check if maps exist
if not os.path.isdir(config["mysql_loadmap"] + config["mysql_database"] +"/"):
    os.mkdir(config["mysql_loadmap"] + config["mysql_database"] +"/")

# Create current and /tmp map locations
sql_map = config["RUN"] +"/sql/"
sql_tmp_map = config["mysql_loadmap"] + config["mysql_database"] +"/"+ args.run_name +"/"
if not os.path.isdir(sql_tmp_map):
    os.mkdir(sql_tmp_map)

# Checks which params are enabled through the contig settings to store different annotations and data from the run into the database.
def run_sql_queries(cnx_handle):
    with cnx_handle.cursor() as cursor:
        # Always creates required tables and insert information into the required/static tables.
        static_insert(cnx_handle)
        # Insert metadata information about the bins.
        if args.store_metadata:
            metadata_insert(cnx_handle)
        # Only insert annotation table if at least one of the tools stored in the table has run.
        if args.store_arCOG or args.store_tigr or args.store_cazy or args.store_signalP or args.store_merops or args.store_transporter or \
            args.store_hydDB or args.store_hmmr1 or args.store_hmmr2 or args.store_blast1 or args.store_blast2:
            annotation_table(cnx_handle)
        # Create and insert information connecting genes to kegg database.
        if args.store_ko:
            ko_table(cnx_handle)
        # Create and insert information connecting genes to cog database.
        if args.store_cog:
            cog_table(cnx_handle)
        # Create and insert information connecting genes to pfam database.
        if args.store_pfam:
            pfam_table(cnx_handle)


# Stores run, bin, contig, feature and diamond data into existing (run/bin) and new tables.
def static_insert(cnx_handle):
    with cnx_handle.cursor() as cursor:
        # Create required tables, contig, feature, and diamond
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ config["mysql_database"] +"`.`Contig_run_"+ args.run_name +"` ("
            "    `contig_id` VARCHAR(100) NOT NULL,"
            "    `bin_id` VARCHAR(45) NOT NULL,"
            "    `oldname` VARCHAR(255) NULL,"
            "    `length` INT NULL,"
            "    `gc` FLOAT NULL,"
            "    `coverage` FLOAT NULL,"
            "    PRIMARY KEY (`contig_id`),"
            "    INDEX `fk_Contig_Bin1_idx` (`bin_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Contig_Bin1`"
            "        FOREIGN KEY (`bin_id`)"
            "        REFERENCES `"+ config["mysql_database"] +"`.`Bin` (`bin_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")
        
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ config["mysql_database"] +"`.`Feature_run_"+ args.run_name +"` ("
            "    `feature_id` VARCHAR(100) NOT NULL,"
            "    `contig_id` VARCHAR(100) NOT NULL,"
            "    `bin_id` VARCHAR(45) NOT NULL,"
            "    `name` VARCHAR(45) NULL,"
            "    `c_start` INT NULL,"
            "    `c_end` INT NULL,"
            "    `direction` VARCHAR(1) NULL,"
            "    `f_type` VARCHAR(45) NULL,"
            "    `definition` VARCHAR(255) NULL,"
            "    PRIMARY KEY (`feature_id`),"
            "    INDEX `fk_Feature_Contig1_idx` (`contig_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Feature_Contig1`"
            "        FOREIGN KEY (`contig_id`)"
            "        REFERENCES `"+ config["mysql_database"] +"`.`Contig_run_"+ args.run_name +"` (`contig_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")

        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ config["mysql_database"] +"`.`Diamond_Annotation_run_"+ args.run_name +"` ("
            "    `feature_annotation_id` INT NOT NULL AUTO_INCREMENT,"
            "    `feature_id` VARCHAR(100) NOT NULL,"
            "    `acc` VARCHAR(20) NOT NULL,"
            "    `description` MEDIUMTEXT BINARY NULL,"
            "    `tax_id` MEDIUMTEXT NULL,"
            "    `taxonomy` MEDIUMTEXT NULL,"
            "    `query_start` INT NULL,"
            "    `query_end` INT NULL,"
            "    `subject_start` INT NULL,"
            "    `subject_end` INT NULL,"
            "    `evalue` FLOAT NULL,"
            "    `bitscore` FLOAT NULL,"
            "    `identity` FLOAT NULL,"
            "    INDEX `fk_Diamond_Annotation_Feature1_idx` (`feature_id` ASC) VISIBLE,"
            "    PRIMARY KEY (`feature_annotation_id`),"
            "    CONSTRAINT `fk_Diamond_Annotation_Feature1`"
            "        FOREIGN KEY (`feature_id`)"
            "        REFERENCES `"+ config["mysql_database"] +"`.`Feature_run_"+ args.run_name +"` (`feature_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")

        # Copy files to tmp map
        shutil.copy(sql_map +"bin.tsv", sql_tmp_map)
        shutil.copy(sql_map +"contig.tsv", sql_tmp_map)
        shutil.copy(sql_map +"feature.tsv", sql_tmp_map)
        shutil.copy(sql_map +"diamond.tsv", sql_tmp_map)
        # Extract information for run only to log the time with NOW().
        with open(sql_map +"run.tsv") as run_file:
            run_file_info = run_file.read().strip("\n").split("\t")
        # Static tables inserts
        insert_table_list = []
        insert_table_list.append("INSERT INTO Run VALUES('"+ run_file_info[0] +"', '"+ run_file_info[1] +"', '"+ run_file_info[2] +"', NOW(), '"+ run_file_info[3] +"', '"+ run_file_info[4] +"', '"+ run_file_info[5] +"', '"+ run_file_info[6] +"');")
        insert_table_list.append("LOAD DATA INFILE '"+ sql_tmp_map +"bin.tsv' into table Bin FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        # Required tables inserts
        insert_table_list.append("LOAD DATA INFILE '"+ sql_tmp_map +"contig.tsv' into table Contig_run_"+ args.run_name +" "
                                    "FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ sql_tmp_map +"feature.tsv' into table Feature_run_"+ args.run_name +" "
                                    "FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ sql_tmp_map +"diamond.tsv' into table Diamond_Annotation_run_"+ args.run_name +" "
                                    "FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        for insert_table in insert_table_list:
            print(insert_table)
            cursor.execute(insert_table)

# Stores metadata into the already existing metadata table.
def metadata_insert(cnx_handle):
    with cnx_handle.cursor() as cursor:
        # Copy file to tmp map
        shutil.copy(sql_map +"metadata.tsv", sql_tmp_map)
        # Insert metadata table
        insert_table_list = ["LOAD DATA INFILE '"+ sql_tmp_map +"metadata.tsv' into table Metadata FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';"]
        for insert_table in insert_table_list:
            print(insert_table)
            cursor.execute(insert_table)

# Stores annotations from several different databases into the created general annotation feature table.
def annotation_table(cnx_handle):
    with cnx_handle.cursor() as cursor:
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ config["mysql_database"] +"`.`Annotation_has_Feature_run_"+ args.run_name +"` ("
            "    `annotation_id` INT NOT NULL,"
            "    `feature_id` VARCHAR(100) NOT NULL,"
            "    `start` INT NOT NULL,"
            "    `end` INT NOT NULL,"
            "    `acc` VARCHAR(45) NOT NULL,"
            "    `description` MEDIUMTEXT NULL,"
            "    `evalue` FLOAT NULL,"
            "    `score` FLOAT NULL,"
            "    `identity` FLOAT NULL,"
            "    PRIMARY KEY (`annotation_id`, `feature_id`, `start`),"
            "    INDEX `fk_Annotation_has_Feature_Feature1_idx` (`feature_id` ASC) VISIBLE,"
            "    INDEX `fk_Annotation_has_Feature_Annotation1_idx` (`annotation_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Annotation_has_Feature_Feature1`"
            "        FOREIGN KEY (`feature_id`)"
            "        REFERENCES `"+ config["mysql_database"] +"`.`Feature_run_"+ args.run_name +"` (`feature_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION,"
            "    CONSTRAINT `fk_Annotation_has_Feature_Annotation1`"
            "        FOREIGN KEY (`annotation_id`)"
            "        REFERENCES `"+ config["mysql_database"] +"`.`Annotation` (`annotation_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")
        
        # Copy file to tmp map
        shutil.copy(sql_map +"annotation.tsv", sql_tmp_map)
        # Annotation tables inserts
        insert_table_list = ["LOAD DATA INFILE '"+ sql_tmp_map +"annotation.tsv' into table Annotation FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';"]
        insert_entry_list = []
        # Looks through all the arguments to see which databases are enabled, to store the TSV files into the database.
        if args.store_arCOG:
            insert_entry_list.append("arCOG.tsv")
        if args.store_tigr:
            insert_entry_list.append("tigr.tsv")
        if args.store_cazy:
            insert_entry_list.append("cazy.tsv")
        if args.store_signalP:
            insert_entry_list.append("signalP.tsv")
        if args.store_merops:
            insert_entry_list.append("merops.tsv")
        if args.store_transporter:
            insert_entry_list.append("transporter.tsv")
        if args.store_hydDB:
            insert_entry_list.append("hydDB.tsv")
        if args.store_hmmr1:
            insert_entry_list.append("db1_hmmr.tsv")
        if args.store_hmmr2:
            insert_entry_list.append("db2_hmmr.tsv")
        if args.store_blast1:
            insert_entry_list.append("db1_blast.tsv")
        if args.store_blast2:
            insert_entry_list.append("db2_blast.tsv")
        for insert_enty in insert_entry_list:
            # Copy files to tmp map
            shutil.copy(sql_map + insert_enty, sql_tmp_map)
            insert_table_list.append("LOAD DATA INFILE '"+ sql_tmp_map + insert_enty +"' into table Annotation_has_Feature_run_"+ args.run_name +" " \
                                        "FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        for insert_table in insert_table_list:
            print(insert_table)
            cursor.execute(insert_table)

# Stores KOs into the created feature KO table.
def ko_table(cnx_handle):
    with cnx_handle.cursor() as cursor:
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ config["mysql_database"] +"`.`Feature_has_KO_run_"+ args.run_name +"` ("
            "    `feature_id` VARCHAR(100) NOT NULL,"
            "    `ko_id` VARCHAR(6) NOT NULL,"
            "    `start` INT NOT NULL,"
            "    `end` INT NOT NULL,"
            "    `evalue` FLOAT NULL,"
            "    `score` FLOAT NULL,"
            "    PRIMARY KEY (`feature_id`, `ko_id`, `start`),"
            "    INDEX `fk_Feature_has_KO_KO1_idx` (`ko_id` ASC) VISIBLE,"
            "    INDEX `fk_Feature_has_KO_Feature1_idx` (`feature_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Feature_has_KO_Feature1`"
            "        FOREIGN KEY (`feature_id`)"
            "        REFERENCES `"+ config["mysql_database"] +"`.`Feature_run_"+ args.run_name +"` (`feature_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION,"
            "    CONSTRAINT `fk_Feature_has_KO_KO1`"
            "        FOREIGN KEY (`ko_id`)"
            "        REFERENCES `"+ config["mysql_database"] +"`.`KO` (`ko_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")
        # Copy file to tmp map
        shutil.copy(sql_map +"ko.tsv", sql_tmp_map)
        # KO tables inserts
        insert_table_list = ["LOAD DATA INFILE '"+ sql_tmp_map +"ko.tsv' into table Feature_has_KO_run_"+ args.run_name +" FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';"]
        for insert_table in insert_table_list:
            print(insert_table)
            cursor.execute(insert_table)

# Stores COGs into the created feature COG table.
def cog_table(cnx_handle):
    with cnx_handle.cursor() as cursor:
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ config["mysql_database"] +"`.`Feature_has_Cog_run_"+ args.run_name +"` ("
            "    `feature_id` VARCHAR(100) NOT NULL,"
            "    `cog_id` VARCHAR(8) NOT NULL,"
            "    `start` INT NOT NULL,"
            "    `end` INT NOT NULL,"
            "    `evalue` FLOAT NULL,"
            "    `score` FLOAT NULL,"
            "    PRIMARY KEY (`feature_id`, `cog_id`, `start`),"
            "    INDEX `fk_Feature_has_Cog_Cog1_idx` (`cog_id` ASC) VISIBLE,"
            "    INDEX `fk_Feature_has_Cog_Feature1_idx` (`feature_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Feature_has_Cog_Feature1`"
            "        FOREIGN KEY (`feature_id`)"
            "        REFERENCES `"+ config["mysql_database"] +"`.`Feature_run_"+ args.run_name +"` (`feature_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION,"
            "    CONSTRAINT `fk_Feature_has_Cog_Cog1`"
            "        FOREIGN KEY (`cog_id`)"
            "        REFERENCES `"+ config["mysql_database"] +"`.`Cog` (`cog_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")
        # Copy file to tmp map
        shutil.copy(sql_map +"cog.tsv", sql_tmp_map)
        # COG tables inserts
        insert_table_list = ["LOAD DATA INFILE '"+ sql_tmp_map +"cog.tsv' into table Feature_has_Cog_run_"+ args.run_name +" FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';"]
        for insert_table in insert_table_list:
            print(insert_table)
            cursor.execute(insert_table)

# Stores Pfams into the created feature Pfam table.
def pfam_table(cnx_handle):
    with cnx_handle.cursor() as cursor:
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ config["mysql_database"] +"`.`Feature_has_Pfam_run_"+ args.run_name +"` ("
            "    `feature_id` VARCHAR(100) NOT NULL,"
            "    `pfam_id` VARCHAR(7) NOT NULL,"
            "    `start` INT NOT NULL,"
            "    `end` INT NOT NULL,"
            "    `evalue` FLOAT NULL,"
            "    `score` FLOAT NULL,"
            "    PRIMARY KEY (`feature_id`, `pfam_id`, `start`),"
            "    INDEX `fk_Feature_has_Pfam_Pfam1_idx` (`pfam_id` ASC) VISIBLE,"
            "    INDEX `fk_Feature_has_Pfam_Feature1_idx` (`feature_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Feature_has_Pfam_Pfam1`"
            "        FOREIGN KEY (`pfam_id`)"
            "        REFERENCES `"+ config["mysql_database"] +"`.`Pfam` (`pfam_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION,"
            "    CONSTRAINT `fk_Feature_has_Pfam_Feature1`"
            "        FOREIGN KEY (`feature_id`)"
            "        REFERENCES `"+ config["mysql_database"] +"`.`Feature_run_"+ args.run_name +"` (`feature_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")
        # Copy file to tmp map
        shutil.copy(sql_map +"pfam.tsv", sql_tmp_map)
        # Pfam tables inserts
        insert_table_list = ["LOAD DATA INFILE '"+ sql_tmp_map +"pfam.tsv' into table Feature_has_Pfam_run_"+ args.run_name +" FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';"]
        for insert_table in insert_table_list:
            print(insert_table)
            cursor.execute(insert_table)

# This function will try connecting to the FAnnP MySQL database.
def connect_to_mysql():
    try:
        return mysql.connector.connect(user=config["mysql_username"], password=config["mysql_password"], host=config["mysql_host"], database=config["mysql_database"])
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password")
            return None
        elif err.errno == errorcode.ER_BAD_DB_ERROR:
            print("Database does not exist")
            return None
        else:
            print(err)
            return None

# Take the MySQL database object and use it to run SQL queries to store all the annotation from the used databases.
cnx = connect_to_mysql()
if cnx and cnx.is_connected():
    print("Connection with the database established.")
    run_sql_queries(cnx)
    cnx.close()
    print("Connection with the database closed.")
elif cnx:
    print("Error, no connection with the server established.")
    cnx.close()

# Makes sure we have an output since it is required for the Snakemake workflow.
with open(args.output, "w") as output_handle:
    output_handle.write("Script execution completed.")

# Remove /tmp/sql/fannp/run_id folder and tsv files
shutil.rmtree(sql_tmp_map)
