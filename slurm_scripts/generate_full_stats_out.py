#!/usr/bin/env python3

import pandas as pd
import argparse as ap
import os
#import sys

parser = ap.ArgumentParser()

parser.add_argument("-d", help="directory of coverage directories")

args = parser.parse_args()

coverage_dir = "coverage"
flagstat_dir = "flagstat"

cov_df_list = []
flag_df_list = []

for dir in os.scandir(args.d):
    cov = os.path.basename(dir)
    print(cov)
    for stat_dir in os.scandir(dir):
        if os.path.basename(stat_dir) == "coverage":
            try:
                for entry in os.scandir(stat_dir):
                    if entry.is_file and entry.name.endswith("abridged_coverage.csv"):
                        df = pd.read_csv(entry.path, delimiter=",", header=0)
                        filename = entry.name.split("_abridged_coverage.csv")[0] 
                        df.insert(0, "coverage_target", cov)
                        df.insert(0, "Filename", filename)
                        cov_df_list.append(df)
            except OSError as e:
                print(f"Issue with coverage data import {e}")
        if os.path.basename(stat_dir) == "flagstat":
            try:
                for entry in os.scandir(stat_dir):
                    if entry.is_file:
                        df = pd.read_csv(entry.path, delimiter="\t")
                        filename = entry.name.split("_aligned_reads_sorted.flagstat")[0]
                        num = df.iloc[6,0]
                        data = {"Filename": [filename], "coverage_target": [cov], "read_Count": [num]}
                        df = pd.DataFrame(data)
                        flag_df_list.append(df)
            except OSError as e:
                print(f"Issue with flagstat data import {e}")

flagstat_df = pd.concat(flag_df_list, ignore_index=True)
coverage_df = pd.concat(cov_df_list, ignore_index=True)
df = coverage_df.merge(flagstat_df, left_on=["Filename", "coverage_target"], right_on=["Filename", "coverage_target"])
df["scaffold_type"] = df["#rname"].str[:2]
print(df)
df.to_csv("out.csv", index=False)
