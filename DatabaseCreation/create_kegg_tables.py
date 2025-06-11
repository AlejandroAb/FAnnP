#!/usr/bin/python
# Tijn 2024-2025 :) (Reference from Nina: https://github.com/ndombrowski/Custom_Code/blob/main/KEGG_to_Pathway_mapping/workflow.qmd)
# Last recorded update: 06-01-2025
#
# !!! -> Make sure to run this using the create_and_insert_fannp.sh script. <- !!!
# This script is necessary for the creation of the tsv files and stores them in the /tmp/sql/static folder.
# Currently only in the /tmp folder data can be loaded into tables with the command 'LOAD DATA'.
#

# Libraries
import pandas as pd
import requests as req
import argparse

# Arguments you can use when running this script.
parser = argparse.ArgumentParser(description='Create KEGG tables for MySQL database storage')
parser.add_argument('-l','--location_name', metavar='folder', required=True, type=str, help='from where the files will be loaded into the database')
args = parser.parse_args()

# Check if the arguments required are supplied by the user.
if args.location_name is None and not os.path.exists(args.location_name):
    print("Location name is missing!")
    sys.exit(-1)

# Extract all target database KEGG tables from the site the information is stored on.
print("Extracting tables from site, estimated time 20 minutes.")
ko_df = pd.read_csv("https://rest.kegg.jp/list/ko", sep='\t', names = ['KO_id', 'KO_desc'])
ko_df[['Symbol', 'KO_desc']] = ko_df['KO_desc'].str.split(';', n=1, expand=True)
ko_df = ko_df[['KO_id', 'Symbol', 'KO_desc']]
module_df = pd.read_csv("https://rest.kegg.jp/list/module", sep='\t', names = ['module_id', 'module_desc'])
module_definitions = []
for i, module_number in enumerate(module_df["module_id"]):
    module_definitions.append(req.get("https://rest.kegg.jp/get/"+ module_number).text.split("DEFINITION")[1].split("\n")[0].lstrip(" ").strip("\r"))
    if i % 100 == 0:
        print(f"Module {i}")
module_df["module_order"] = module_definitions
pathway_df = pd.read_csv("https://rest.kegg.jp/list/pathway", sep='\t', names = ['pathway_map', 'pathway_desc'])
ec_df = pd.read_csv("https://rest.kegg.jp/list/ec", sep='\t', names = ['EC_id', 'EC_desc'])

# Extract linked tables.
print("Extracting linked tables from site.")
module_has_ko_df = pd.read_csv("https://rest.kegg.jp/link/ko/module", sep='\t', names = ['module_id', 'KO_id'])
pathway_has_ko_df = pd.read_csv("https://rest.kegg.jp/link/ko/pathway", sep='\t', names = ['pathway_map', 'KO_id'])
pathway_has_module_df = pd.read_csv("https://rest.kegg.jp/link/module/pathway", sep='\t', names = ['pathway_map', 'module_id'])
ko_has_ec_df = pd.read_csv("https://rest.kegg.jp/link/ec/ko", sep='\t', names = ['KO_id', 'EC_id'])

# Modify the tables to fit it with our ERD.
print("Modifying tables for database.")
module_has_ko_df["module_id"] = module_has_ko_df["module_id"].str.replace("md:", "")
module_has_ko_df["KO_id"] = module_has_ko_df["KO_id"].str.replace("ko:", "")
module_has_ko_df = module_has_ko_df.drop_duplicates()
pathway_has_ko_df = pathway_has_ko_df.query("not pathway_map.str.contains('path:ko.+')")
pathway_has_ko_df["pathway_map"] = pathway_has_ko_df["pathway_map"].str.replace("path:", "")
pathway_has_ko_df["KO_id"] = pathway_has_ko_df["KO_id"].str.replace("ko:", "")
pathway_has_ko_df = pathway_has_ko_df.drop_duplicates()
pathway_has_module_df["pathway_map"] = pathway_has_module_df["pathway_map"].str.replace("path:", "")
pathway_has_module_df["module_id"] = pathway_has_module_df["module_id"].str.replace("md:", "")
ko_has_ec_df["KO_id"] = ko_has_ec_df["KO_id"].str.replace("ko:", "")
ko_has_ec_df["EC_id"] = ko_has_ec_df["EC_id"].str.replace("ec:", "")

# Store the tables inside the target map location.
print("Exporting tables.")
ko_df.to_csv(args.location_name +"ko.tsv", header=None, index=False, sep='\t')
module_df.to_csv(args.location_name +"module.tsv", header=None, index=False, sep='\t')
pathway_df.to_csv(args.location_name +"pathway.tsv", header=None, index=False, sep='\t')
ec_df.to_csv(args.location_name +"ec.tsv", header=None, index=False, sep='\t')
module_has_ko_df.to_csv(args.location_name +"module_has_ko.tsv", header=None, index=False, sep='\t')
pathway_has_ko_df.to_csv(args.location_name +"pathway_has_ko.tsv", header=None, index=False, sep='\t')
pathway_has_module_df.to_csv(args.location_name +"pathway_has_module.tsv", header=None, index=False, sep='\t')
ko_has_ec_df.to_csv(args.location_name +"ko_has_ec.tsv", header=None, index=False, sep='\t')

