#!/usr/bin/python
# Tijn 2024-2025 :)
# Last recorded update: 29-10-2024
#
# !!! -> Make sure to run this using the create_and_insert_fannp.sh script. <- !!!
# This script is necessary for the creation of the tsv files and stores them in the /tmp/sql/static folder.
# Currently only in the /tmp folder data can be loaded into tables with the command 'LOAD DATA'.
#

# Libraries
import sys
import os
import argparse
import mysql.connector
from mysql.connector import errorcode

# Arguments you can use when running this script.
parser = argparse.ArgumentParser(description='Inserts tsv files into MySQL database')
parser.add_argument('-u','--user_name', metavar='str', required=True, type=str, help='username for MySQL account')
parser.add_argument('-p','--password', metavar='str', required=True, type=str, help='password for MySQL account')
parser.add_argument('-d','--database_name', metavar='str', required=True, type=str, help='database name for MySQL database')
parser.add_argument('-l','--location_name', metavar='folder', required=True, type=str, help='from where the files will be loaded into the database')
parser.add_argument('--reload_databases', action="store_true", help='to reload the database tables KEGG, COG, and Pfam')
parser.add_argument('--skip_kegg', action="store_false", help='to skip kegg when reloading the databases')
args = parser.parse_args()

print("")
# Check if the arguments required are supplied by the user.
if args.user_name is None:
    print("Username is missing!")
    sys.exit(-1)
if args.password is None:
    print("Password is missing!")
    sys.exit(-1)
if args.database_name is None:
    print("Database name is missing!")
    sys.exit(-1)
if args.location_name is None and not os.path.exists(args.location_name):
    print("Location name is missing!")
    sys.exit(-1)

# Only empties all static tables and then repopulates them
def reload_static_tables(cnx_handle):
    with cnx_handle.cursor() as cursor:
        # Delete all entries
        # Used to be "TRUNCATE TABLE `"+ args.database_name +"`.`Table`;";
        cursor.execute("DELETE FROM `"+ args.database_name +"`.`Cog_has_Cog_family`;")
        cursor.execute("DELETE FROM `"+ args.database_name +"`.`Cog`;")
        cursor.execute("DELETE FROM `"+ args.database_name +"`.`Cog_family`;")

        # Delete all entries
        cursor.execute("DELETE FROM `"+ args.database_name +"`.`Pfam_clan`;")
        cursor.execute("DELETE FROM `"+ args.database_name +"`.`Pfam`;")

        # Delete all entries
        if args.skip_kegg:
            cursor.execute("DELETE FROM `"+ args.database_name +"`.`KO`;")
            cursor.execute("DELETE FROM `"+ args.database_name +"`.`Module`;")
            cursor.execute("DELETE FROM `"+ args.database_name +"`.`Pathway`;")
            cursor.execute("DELETE FROM `"+ args.database_name +"`.`EC`;")
            cursor.execute("DELETE FROM `"+ args.database_name +"`.`Module_has_KO`;")
            cursor.execute("DELETE FROM `"+ args.database_name +"`.`Pathway_has_KO`;")
            cursor.execute("DELETE FROM `"+ args.database_name +"`.`Pathway_has_Module`;")
            cursor.execute("DELETE FROM `"+ args.database_name +"`.`KO_has_EC`;")

        # Add all entries back
        insert_cog_information(cnx_handle)
        insert_pfam_information(cnx_handle)
        if args.skip_kegg:
            insert_ko_information(cnx_handle)


# Creates or deletes the structure of all static tables and then populates them
def create_static_tables(cnx_handle):
    with cnx_handle.cursor() as cursor:
        # Drop run, bin, metadata and annotation tables if exist
        """cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Run`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Bin`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Metadata`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Annotation`;")"""

        # Drop COG tables if exist
        """cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Cog_has_Cog_family`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Cog`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Cog_family`;")"""

        # Drop Pfam tables if exist
        """cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Pfam_clan`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Pfam`;")"""

        # Drop KEGG tables if exist
        """cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`KO`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Module`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Pathway`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`EC`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Module_has_KO`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Pathway_has_KO`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Pathway_has_Module`;")
        cursor.execute("DROP TABLE IF EXISTS `"+ args.database_name +"`.`KO_has_EC`;")"""

        # Create run, bin, metadata and annotation tables
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Run` ("
            "    `run_id` VARCHAR(8) NOT NULL,"
            "    `name` VARCHAR(45) NOT NULL,"
            "    `start` DATETIME NULL,"
            "    `end` DATETIME NULL DEFAULT CURRENT_TIMESTAMP,"
            "    `wd` VARCHAR(255) NULL,"
            "    `version` VARCHAR(6) NULL,"
            "    `user` VARCHAR(45) NULL,"
            "    `od` VARCHAR(255) NULL,"
            "    PRIMARY KEY (`run_id`))"
            "ENGINE = MyISAM;")

        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Bin` ("
            "    `bin_id` VARCHAR(45) NOT NULL,"
            "    `run_id` VARCHAR(8) NOT NULL,"
            "    `bin_path` VARCHAR(255) NULL,"
            "    `bin_ext` VARCHAR(45) NULL,"
            "    `taxonomy` MEDIUMTEXT NULL,"
            "    `bin_length` INT NULL,"
            "    `bin_avg_gc` FLOAT NULL,"
            "    `bin_avg_cov` FLOAT NULL,"
            "    PRIMARY KEY (`bin_id`, `run_id`),"
            "    INDEX `fk_Bin_Run1_idx` (`run_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Bin_Run1`"
            "        FOREIGN KEY (`run_id`)"
            "        REFERENCES `"+ args.database_name +"`.`Run` (`run_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")

        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Metadata` ("
            "    `metadata_id` INT NOT NULL AUTO_INCREMENT,"
            "    `bin_id` VARCHAR(45) NOT NULL,"
            "    `run_id` VARCHAR(8) NOT NULL,"
            "    `name` VARCHAR(45) NULL,"
            "    `value` MEDIUMTEXT NULL,"
            "    PRIMARY KEY (`metadata_id`),"
            "    INDEX `fk_Metadata_Bin1_idx` (`bin_id` ASC, `run_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Metadata_Bin1`"
            "        FOREIGN KEY (`bin_id` , `run_id`)"
            "        REFERENCES `"+ args.database_name +"`.`Bin` (`bin_id` , `run_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")

        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Annotation` ("
            "    `annotation_id` INT NOT NULL,"
            "    `run_id` VARCHAR(8) NOT NULL,"
            "    `database` VARCHAR(45) NOT NULL,"
            "    `db_version` VARCHAR(45) NOT NULL,"
            "    `def_file` MEDIUMTEXT NULL,"
            "    `db_file` MEDIUMTEXT NULL,"
            "    INDEX `fk_Annotation_Run1_idx` (`run_id` ASC) VISIBLE,"
            "    PRIMARY KEY (`annotation_id`, `run_id`),"
            "    CONSTRAINT `fk_Annotation_Run1`"
            "        FOREIGN KEY (`run_id`)"
            "        REFERENCES `"+ args.database_name +"`.`Run` (`run_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")

        # Create COG tables
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Cog_family` ("
            "    `cog_family_id` VARCHAR(1) NOT NULL,"
            "    `cog_family_name` MEDIUMTEXT NULL,"
            "    PRIMARY KEY (`cog_family_id`))"
            "ENGINE = MyISAM;")
        
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Cog` ("
            "    `cog_id` VARCHAR(7) NOT NULL,"
            "    `cog_desc` MEDIUMTEXT NULL,"
            "    PRIMARY KEY (`cog_id`))"
            "ENGINE = MyISAM;")
        
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Cog_has_Cog_family` ("
            "    `cog_family_id` VARCHAR(1) NOT NULL,"
            "    `cog_id` VARCHAR(7) NOT NULL,"
            "    PRIMARY KEY (`cog_family_id`, `cog_id`),"
            "    INDEX `fk_Cog_has_Cog_family_Cog1_idx` (`cog_id` ASC) VISIBLE,"
            "    INDEX `fk_Cog_has_Cog_family_Cog_family1_idx` (`cog_family_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Cog_has_Cog_family_Cog_family1`"
            "        FOREIGN KEY (`cog_family_id`)"
            "        REFERENCES `"+ args.database_name +"`.`Cog_family` (`cog_family_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION,"
            "    CONSTRAINT `fk_Cog_has_Cog_family_Cog1`"
            "        FOREIGN KEY (`cog_id`)"
            "        REFERENCES `"+ args.database_name +"`.`Cog` (`cog_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")

        # Create Pfam tables
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Pfam_clan` ("
            "    `pfam_clan_id` VARCHAR(7) NOT NULL,"
            "    `pfam_clan_name` VARCHAR(45) NULL,"
            "    PRIMARY KEY (`pfam_clan_id`))"
            "ENGINE = MyISAM;")
        
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Pfam` ("
            "    `pfam_id` VARCHAR(7) NOT NULL,"
            "    `pfam_acc` VARCHAR(45) NULL,"
            "    `pfam_desc` VARCHAR(255) NULL,"
            "    `pfam_type` VARCHAR(45) NULL,"
            "    `pfam_clan_id` VARCHAR(7) NULL,"
            "    PRIMARY KEY (`pfam_id`),"
            "    INDEX `fk_Pfam_Pfam_clan1_idx` (`pfam_clan_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Pfam_Pfam_clan1`"
            "        FOREIGN KEY (`pfam_clan_id`)"
            "        REFERENCES `"+ args.database_name +"`.`Pfam_clan` (`pfam_clan_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")
        
        # Create KEGG tables
        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`KO` ("
            "    `ko_id` VARCHAR(6) NOT NULL,"
            "    `ko_symbol` VARCHAR(255) NULL,"
            "    `ko_desc` VARCHAR(255) NULL,"
            "    PRIMARY KEY (`ko_id`))"
            "ENGINE = MyISAM;")

        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Module` ("
            "    `module_id` VARCHAR(6) NOT NULL,"
            "    `module_desc` VARCHAR(255) NULL,"
            "    `module_order` MEDIUMTEXT NULL,"
            "    PRIMARY KEY (`module_id`))"
            "ENGINE = MyISAM;")

        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Pathway` ("
            "    `pathway_id` VARCHAR(8) NOT NULL,"
            "    `pathway_desc` VARCHAR(255) NULL,"
            "    PRIMARY KEY (`pathway_id`))"
            "ENGINE = MyISAM;")

        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`EC` ("
            "    `ec_id` VARCHAR(16) NOT NULL,"
            "    `ec_name` MEDIUMTEXT NULL,"
            "    PRIMARY KEY (`ec_id`))"
            "ENGINE = MyISAM;")

        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Module_has_KO` ("
            "    `module_id` VARCHAR(6) NOT NULL,"
            "    `ko_id` VARCHAR(6) NOT NULL,"
            "    PRIMARY KEY (`ko_id`, `module_id`),"
            "    INDEX `fk_KO_has_Module_Module1_idx` (`module_id` ASC) VISIBLE,"
            "    INDEX `fk_KO_has_Module_KO1_idx` (`ko_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_KO_has_Module_KO1`"
            "        FOREIGN KEY (`ko_id`)"
            "        REFERENCES `"+ args.database_name +"`.`KO` (`ko_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION,"
            "    CONSTRAINT `fk_KO_has_Module_Module1`"
            "        FOREIGN KEY (`module_id`)"
            "        REFERENCES `"+ args.database_name +"`.`Module` (`module_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")

        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Pathway_has_KO` ("
            "    `pathway_id` VARCHAR(8) NOT NULL,"
            "    `ko_id` VARCHAR(6) NOT NULL,"
            "    PRIMARY KEY (`pathway_id`, `ko_id`),"
            "    INDEX `fk_Pathway_has_KO_KO1_idx` (`ko_id` ASC) VISIBLE,"
            "    INDEX `fk_Pathway_has_KO_Pathway1_idx` (`pathway_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Pathway_has_KO_Pathway1`"
            "        FOREIGN KEY (`pathway_id`)"
            "        REFERENCES `"+ args.database_name +"`.`Pathway` (`pathway_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION,"
            "    CONSTRAINT `fk_Pathway_has_KO_KO1`"
            "        FOREIGN KEY (`ko_id`)"
            "        REFERENCES `"+ args.database_name +"`.`KO` (`ko_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")

        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`Pathway_has_Module` ("
            "    `pathway_id` VARCHAR(8) NOT NULL,"
            "    `module_id` VARCHAR(6) NOT NULL,"
            "    PRIMARY KEY (`pathway_id`, `module_id`),"
            "    INDEX `fk_Pathway_has_Module_Module1_idx` (`module_id` ASC) VISIBLE,"
            "    INDEX `fk_Pathway_has_Module_Pathway1_idx` (`pathway_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_Pathway_has_Module_Pathway1`"
            "        FOREIGN KEY (`pathway_id`)"
            "        REFERENCES `"+ args.database_name +"`.`Pathway` (`pathway_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION,"
            "    CONSTRAINT `fk_Pathway_has_Module_Module1`"
            "        FOREIGN KEY (`module_id`)"
            "        REFERENCES `"+ args.database_name +"`.`Module` (`module_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")

        cursor.execute(
            "CREATE TABLE IF NOT EXISTS `"+ args.database_name +"`.`KO_has_EC` ("
            "    `ko_id` VARCHAR(6) NOT NULL,"
            "    `ec_id` VARCHAR(16) NOT NULL,"
            "    PRIMARY KEY (`ko_id`, `ec_id`),"
            "    INDEX `fk_KO_has_EC_EC1_idx` (`ec_id` ASC) VISIBLE,"
            "    CONSTRAINT `fk_KO_has_EC_KO1`"
            "        FOREIGN KEY (`ko_id`)"
            "        REFERENCES `"+ args.database_name +"`.`KO` (`ko_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION,"
            "    CONSTRAINT `fk_KO_has_EC_EC1`"
            "        FOREIGN KEY (`ec_id`)"
            "        REFERENCES `"+ args.database_name +"`.`EC` (`ec_id`)"
            "        ON DELETE NO ACTION"
            "        ON UPDATE NO ACTION)"
            "ENGINE = MyISAM;")

        # Add all entries back
        insert_cog_information(cnx_handle)
        insert_pfam_information(cnx_handle)
        insert_ko_information(cnx_handle)


# Load in COG data to their respective tables
def insert_cog_information(cnx_handle):
    with cnx_handle.cursor() as cursor:
        insert_table_list = []
        # COG table inserts
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"cog_family.tsv' into table Cog_family FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"cog.tsv' into table Cog FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"cog_has_cog_family.tsv' into table Cog_has_Cog_family FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        for insert_table in insert_table_list:
            print(insert_table)
            cursor.execute(insert_table)


# Load in Pfam data to their respective tables
def insert_pfam_information(cnx_handle):
    with cnx_handle.cursor() as cursor:
        insert_table_list = []
        # Pfam table inserts
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"pfam_clan.tsv' into table Pfam_clan FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"pfam.tsv' into table Pfam FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        for insert_table in insert_table_list:
            print(insert_table)
            cursor.execute(insert_table)

# Load in KEGG data to their respective tables
def insert_ko_information(cnx_handle):
    with cnx_handle.cursor() as cursor:
        insert_table_list = []
        # KEGG table inserts
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"ko.tsv' into table KO FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"module.tsv' into table Module FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"pathway.tsv' into table Pathway FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"ec.tsv' into table EC FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"module_has_ko.tsv' into table Module_has_KO FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"pathway_has_ko.tsv' into table Pathway_has_KO FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"pathway_has_module.tsv' into table Pathway_has_Module FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        insert_table_list.append("LOAD DATA INFILE '"+ args.location_name +"ko_has_ec.tsv' into table KO_has_EC FIELDS TERMINATED BY '\\t' LINES TERMINATED BY '\\n';")
        for insert_table in insert_table_list:
            print(insert_table)
            cursor.execute(insert_table)

# This function will try connecting to the FAnnP (designated) MySQL database.
def connect_to_mysql():
    try:
        return mysql.connector.connect(user=args.user_name, password=args.password, host='localhost', database=args.database_name)
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

# Take the MySQL database object and use it to run SQL queries to store all target information tables into the database.
cnx = connect_to_mysql()
if cnx and cnx.is_connected():
    print("Connection with the database established.")
    if args.reload_databases:
        print("Reloading all static databases...")
        reload_static_tables(cnx)
    else:
        print("Constructing database structure...")
        create_static_tables(cnx)
    cnx.close()
    print("Connection with the database closed.")
elif cnx:
    print("Error, no connection with the server established.")
    cnx.close()
print("")
