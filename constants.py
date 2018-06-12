#!/usr/bin/env python3

# -------------------------------------------------------------------

# File: constants.py
# Author: Harrison Inocencio
# Date: 06-07-18
# Purpose: Contains all constant definitions needed for the juicer.py
#		   program

# Notes:
# 1.
# 2.
# 3.
# 4.
# 5.

# TODO:
# 1. Figure out how to implement update 
# 2.
# 3.
# 4.
# 5.

# -------------------------------------------------------------------

import argparse

# Default parameters
default_mismatches = 1
default_groupings = 30
default_grep_multi = 3
default_output = "juicer_out"

# Gap configs
one_rep_min = 5
two_rep_min = 10
max_size = 14

# Program configs
disection_toggle = 0 # Toggle disection mode, 1=ON, 0=OFF | Currently NA
start_grp = 1 # Currently Unused
frame_multi = 2 # Don't change this
zfill_amount = 4 # Amount of zero padding in file names
min_motif_len = 5 # Minimum allowed motif length
grep_temp = ".grep_temp_juicer"
out_prefix = "juicer_group"

# Argparse constants
fastq_directory_help = 'Path to directory containing fastq files to be juiced. Searches each file for set motif and groups them based on the number of times this motif repeats from the start of the sequence'
motif_help = 'Nucleic motif to be used as the repeat base'
type_help = "Output in fasta, or fastq format"
group_count_help = 'Specifies the number of groupings to be used (This also codes for the maximum repeat count)(Default = 30)' 
mismatches_help = 'Number of mismatches to be allowed in each instance of the motif (Default = 1)'
output_name_help = 'Name of the output directory for juicer (Default = juicer_out)'
grep_multi_help = 'Juicer uses the grep program to initial extract it\'s motif sequences. This option controls the number of consecutive copies of the motif grep will search for (Default = 3)'
type_choices = ["fasta", "fastq"]

# set_usr_const function
# > uses argparse to set all input parameters
# > returns args datastruct
def set_usr_const():
 # Argument parser
	parser = argparse.ArgumentParser()
	parser.add_argument('fastq_directory', help=fastq_directory_help)
	parser.add_argument('motif', help=motif_help)
	parser.add_argument('type', help=type_help, choices=type_choices)
	parser.add_argument('-g', '--number-of-groupings', default=default_groupings, help=group_count_help, type=int)
	# Update for new method
	parser.add_argument('-m', '--mismatches', default=default_mismatches, help=mismatches_help, type=int)
	parser.add_argument('-o', '--output-name', default=default_output, help=output_name_help)
	parser.add_argument('--grep-multi', default=default_grep_multi, help=grep_multi_help, type=int)
	args = parser.parse_args()

	return args
