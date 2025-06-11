#!/usr/bin/python
# Tijn 2024-2025 :)
# Last recorded update: 02-2025
#
# python for annotation to read config file
# !!! -> New databases added to the config file might need to be added to the potential tables list <- !!!
#

# Libraries
import sys
import os
import argparse
import yaml


# Arguments you can use when running this script.
parser = argparse.ArgumentParser(description='removes bad domains, copied domains, and overlapped domains')
parser.add_argument('-r','--run_name', metavar='str', required = True, type=str, help='the id of the run')
parser.add_argument('-c','--config_name', metavar='file', required = True, type=str, help='config file used by the snakemake pipeline')
parser.add_argument('-o','--output', metavar='file', required = True, type=str, help='tsv file containing database information from the config file')
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

# !!! -> All relevant database entries, if a new one is added it should also be added to this list <- !!!
potential_tables = {"01": "arCOG_hmmr", "02": "tigr_hmmr", "03": "cazy_hmmr", "04": "signalP",
                    "05": "merops_blast", "06": "transporterDB_blast", "07": "hydDB_blast",
                    "08": "db1_hmmr", "09": "db2_hmmr", "10": "db1_blast", "11": "db2_blast"}

# Add missing config entries in the config dictionary
required_entry_vars = ["database_name", "database_version", "database_names", "database"]

# Store all the different database entries with relevant information into the output file
with open(args.output, "w") as output_handle1:
    for num in (potential_tables.keys()):
        if config[potential_tables[num]]["annotate"] == "T":
            for required_entry in required_entry_vars:
                if required_entry not in config[potential_tables[num]].keys():
                    config[potential_tables[num]][required_entry] = ""
            output_handle1.write(num + "\t" +
                        args.run_name + "\t" +
                        config[potential_tables[num]]["database_name"] + "\t" +
                        config[potential_tables[num]]["database_version"] + "\t" +
                        config[potential_tables[num]]["database_names"] + "\t" +
                        config[potential_tables[num]]["database"] + "\n")
