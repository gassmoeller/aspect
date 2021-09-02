#!/bin/bash

python ../../../contrib/python/reformat_statistics_gnuplot.py output/statistics
gnuplot plot_results.gnuplot