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
from io import StringIO
import pandas as pd
from matplotlib import pyplot as plt


# Arguments you can use when running this script.
parser = argparse.ArgumentParser(description='removes bad domains, copied domains, and overlapped domains')
parser.add_argument('-b','--benchmark_map', metavar='str', required = True, type=str, help='the map where the benchmarks of a run are stored')
parser.add_argument('-o','--output', metavar='file', required = True, type=str, help='file which contains the pie chart')
args = parser.parse_args()


# Check if the arguments required are supplied by the user (will do it automatically in the pipeline).
if args.benchmark_map is None:
    print("Run name is missing!")
    sys.exit(-1)
if args.output is None and not os.path.exists(args.output):
    print("Output file is missing!")
    sys.exit(-1)


# List all benchmark files.
dir_list = os.listdir(args.benchmark_map)


# Make a Pandas table with all the benchmarks.
all_benchmarks_string = "rule\ts\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time"
for benchmark_file in dir_list:
    with open(f"{args.benchmark_map}/{benchmark_file}") as input_handle:
        name = benchmark_file.replace('.benchmark', '')
        info = input_handle.readlines()[-1].strip("\n").replace('NA', '0')
        all_benchmarks_string += f"\n{name}\t{info}"
benchmark_df = pd.read_csv(StringIO(all_benchmarks_string), sep ="\t")

# Grab the rules with the highest run time and display the rest in the other rules row.
top_benchmarks_listed = benchmark_df.set_index(["rule"]).sort_values(["s"], ascending=False).index[:7]
benchmark_df.loc[~benchmark_df["rule"].isin(top_benchmarks_listed), "rule"] = "other rules"
benchmark_df = benchmark_df[["rule", "s"]].groupby(["rule"]).aggregate({"rule": "first", "s": "sum"})
print("\nIf you see that the pie chart cannot be displayed here, don't mind it.\nThe pie chart should be stored in the directory you run the script in.\nOutput file: {args.output}\n")

# Create the pie chart.
ax = plt.subplot(111)
patches, texts = ax.pie(benchmark_df["s"], startangle=90, radius=1.2)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height * 0.6])
ax.legend(patches, benchmark_df["rule"], loc='center left', bbox_to_anchor=(-0.3, 1.5), fontsize=8)
total_hours = int(benchmark_df["s"].sum()/3600)
total_minutes = int(benchmark_df["s"].sum()%3600/60)
total_seconds = int(benchmark_df["s"].sum()%60)
ax.set_title(f"Total run time in h:m:s = {total_hours}:{total_minutes}:{total_seconds}")
plt.savefig(args.output)
