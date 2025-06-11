#!/usr/bin/python
# Tijn 2024-2025 :)
# Last recorded update: 02-2025
#
# hmmr --domtblout output parser
#

# Libraries
import sys
import os
import argparse
import parse_hmmr_domains_function


# Arguments you can use when running this script.
parser = argparse.ArgumentParser(description='removes bad domains, copied domains, and overlapped domains')
parser.add_argument('-i','--input', metavar='file', required = True, type=str, help='processed hmmscan --domtblout file which has to be sorted on gene and then highest independent domain E-value')
parser.add_argument('-p','--percentage', metavar='int', required = True, type=str, help='var containing the percentage overlap we allow from smallest domain')
parser.add_argument('-c','--cutoff', metavar='float', required = True, type=str, help='var containing the cutoff value of the independent e-value of domains')
parser.add_argument('-d','--domain_index', metavar='int', required = True, type=str, help='var containing the index of the domain name')
parser.add_argument('-o','--output', metavar='file', required = True, type=str, help='file that will contain a list of proteins validated domains')
# These two arguments are both required to be filled in or -m == "F" (generate_sql_database in config.yaml)
parser.add_argument('-m', '--make_sql_database', action="store_true", help='whether a database is going to get created')
# Required if make_sql_database = True
parser.add_argument('--sql', metavar='file', default='', required = False, type=str, help='file that will contain a list of proteins validated domains used as tsv input for database')
args = parser.parse_args()


# Check if the arguments required are supplied by the user (will do it automatically in the pipeline).
if args.input is None and not os.path.exists(args.input):
	print("Input file is missing!")
	sys.exit(-1)
if args.percentage is None or args.cutoff is None or args.domain_index is None:
	print("Variable is missing!")
	sys.exit(-1)
if args.output is None and not os.path.exists(args.output):
	print("Output file is missing!")
	sys.exit(-1)
if args.make_sql_database and args.sql is None and not os.path.exists(args.sql):
	print("SQL output file is missing!")
	sys.exit(-1)


# Execute parse function.
parse_hmmr_domains_function.parse_hmmr_domains(args.input, args.percentage, args.cutoff, args.domain_index, args.output, args.make_sql_database, args.sql)


