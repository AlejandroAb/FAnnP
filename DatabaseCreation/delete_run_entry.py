#!/usr/bin/python
#
# Deletes all entries and tables connected to a specific run ID.
#

# Libraries
import sys
import argparse
import mysql.connector
from mysql.connector import errorcode

# Arguments you can use when running this script.
parser = argparse.ArgumentParser(description='Delete a specific run')
parser.add_argument('-r','--run_name', metavar='str', required = True, type=str, help='the id of the run')
parser.add_argument('-d','--database_name', metavar='str', required=True, type=str, help='database name for MySQL database')
args = parser.parse_args()

# Check if the arguments required are supplied by the user.
if args.run_name is None:
    print("Run name is missing!")
    sys.exit(-1)

# Delete all tables that is linked to our target run ID.
def delete_all_run_entries(cnx_handle):
    with cnx_handle.cursor() as cursor:
        delete_table_list = []
        # Delete dynamic tables
        delete_table_list.append("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Feature_has_Pfam_run_"+ args.run_name +"`;")
        delete_table_list.append("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Feature_has_Cog_run_"+ args.run_name +"`;")
        delete_table_list.append("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Feature_has_KO_run_"+ args.run_name +"`;")
        delete_table_list.append("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Annotation_has_Feature_run_"+ args.run_name +"`;")
        # Delete required tables
        delete_table_list.append("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Diamond_Annotation_run_"+ args.run_name +"`;")
        delete_table_list.append("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Feature_run_"+ args.run_name +"`;")
        delete_table_list.append("DROP TABLE IF EXISTS `"+ args.database_name +"`.`Contig_run_"+ args.run_name +"`;")
        for delete_table in delete_table_list:
            print(delete_table)
            cursor.execute(delete_table)

        delete_entry_list = []
        # Delete all run name entries in static tables
        delete_entry_list.append("DELETE FROM Annotation WHERE run_id='"+ args.run_name +"';")
        delete_entry_list.append("DELETE FROM Metadata WHERE run_id='"+ args.run_name +"';")
        delete_entry_list.append("DELETE FROM Bin WHERE run_id='"+ args.run_name +"';")
        delete_entry_list.append("DELETE FROM Run WHERE run_id='"+ args.run_name +"';")
        for delete_entry in delete_entry_list:
            print(delete_entry)
            cursor.execute(delete_entry)

# This function will try connecting to the FAnnP MySQL database.
def connect_to_mysql():
    try:
        return mysql.connector.connect(user="faanp_user", password="FaanP2025!", host='localhost', database=args.database_name)
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

# Take the MySQL database object and use it to run SQL queries to remove all target run ID tables from the database.
cnx = connect_to_mysql()
if cnx and cnx.is_connected():
    print("Connection with the database established.")
    delete_all_run_entries(cnx)
    cnx.close()
    print("Connection with the database closed.")
elif cnx:
    print("Error, no connection with the server established.")
    cnx.close()
