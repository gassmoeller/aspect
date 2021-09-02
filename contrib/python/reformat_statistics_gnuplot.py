""" A minimal python script to reformat one or more ASPECT statistics files for gnuplot."""

import sys
import pandas
import csv

import aspect_data

assert(len(sys.argv) > 1), "Please provide at least one ASPECT statistics file to process"

for file in sys.argv[1:]:
    print("Processing file: ", file)
    statistics = aspect_data.read_statistics(file)
    statistics.to_csv(file + ".dat", sep=" ", index=False, quoting=csv.QUOTE_NONNUMERIC)

print("Finished processing.")

