#!/bin/python

# This interpolation script takes the HeFESTo table in its original format,
# subsamples it, reorders the columns, and adds a header so that it can be used
# by the AsciiDataLookup class inside ASPECT.

import numpy as np

input_filename = "opxtable_s.aspect.full"
output_filename = "opxtable_s.aspect"

data = np.genfromtxt(input_filename)
table = data.reshape(1401,241,9)

modified = table[::10,:,:]
modified = modified[:,::5,:]
modified = modified.reshape(modified.shape[0]*modified.shape[1],9)

permutation = [2, 1, 3, 4, 5, 6, 7, 8, 0]
idx = np.empty_like(permutation)
idx[permutation] = np.arange(len(permutation))
modified = modified[:, idx]

headerfile = open("header_table.txt","r")
headertext = headerfile.read()
headertext = headertext[:-1]

np.savetxt(output_filename, modified,fmt='%.5g', header=headertext, comments='')
