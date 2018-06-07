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

# Default parameters
default_mismatches = 1
default_groupings = 30
default_grep_multi = 3
default_output = "juicer_out"

# Program configs
start_grp = 1 # Currently Unused
frame_multi = 2 # Don't change this
zfill_amount = 4 # Amount of zero padding in file names
min_motif_len = 5 # Minimum allowed motif length
grep_temp = ".grep_temp_juicer"
out_prefix = "juicer_group"

def update():
	return 0
