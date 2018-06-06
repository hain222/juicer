#!/usr/bin/env python3

# -------------------------------------------------------------------

# File: juicer.py
# Author: Harrison Inocencio
# Date: 06-01-18
# Purpose: Extract sequences from a read library according to if they contain a specific motif, 
# then sort according to the number of sequence repeats. Inputs a directory containing fastq 
# reads, and outputs n number of files according to how many repeat groupings are desired

# Usage: NA

# Compilation: NA

# Notes:
# 1. ONLY COUNTS BEGINNING REPEATS, i.e with motif AA, the seq AAAAAAGGAA would be sorted into
#	 group 3
# 2. Prefix's that are subsequences of the motif are allowed not counted. i.e with motif AGT
#	 the seq GTAGTAGTAAAA would be sorted into group 2
# 3. Can set the number of allowed mismatches in a repeat, this applies to prefix repeats as well
#	 i.e with motif AGT and mismatch set to 1, the seq GGAGTATTAGTATGGGG would be sorted into 
#	 group 3
# 4. Does not handle insertions/deletions, moves by frame of reference
# 5. prefixes will also be subject to mismatching, complex system for finding starting frame
# 6. sequences less than the 2n base pairs, with n being the length of the motif, will not
#    evaluate correctly since the current function needs at least that much for its frame
#	 evaluation 

# TODO:
# 1. Detail Notes better/Compelte descriptions
# 2. Write parser help descriptions
# 3. Error check all input parameters
# 4. Add option for fasta, or have a specific optional for fasta and fastq
# 5. May need to add clean-up function for errors
# 6. Add option for retaining all generated files not just outputs
# 7. remove hardcoded decoding command
# 8. Find better solution for grep parsing besides writing a file...
# 9. Possibly add start_grp to parameters, (implement a range), so major changes may be
#	 required for this 
# 10. Should handle single base insertions/deletions in the future?
# 11. Complex function to determine starting frame, show to expert
# 12. Make option for a log/dump file that logs the decisions made for each read
# 13. Compartamentalize the mismatch code
# 14. Add constraints to the motif parameter >= 3 etc.
# 15. Frame evaluation reaches a lot, needs more error checking (long reads, etc)
# 16. Possible base the last piece of frame evaluation off of the sum of the target frame and its
#	  sister frame
# 17. Improve frame evaluation error
# 18. Check reads before hand, throw out ones shorter than 2n where n is the length of the motif
# 19. Need a better error handling system

# -------------------------------------------------------------------

#Globals
default_mismatches = 1
default_groupings = 30
start_grp = 1
frame_multi = 2
grep_temp = ".grep_temp_juicer"

import os
import argparse
import subprocess
from Bio import SeqIO

# init function
# > creates initial environment and startup files
def init():

	return 0

# mop function
# > clean up function, deletes temporary files etc..
def mop():

	return 0

# crit_error function
# > handles unrecoverable program errors
def crit_error(error):
	print("CRITICAL ERROR:", error)
	exit(1)

# grep_parse function
# > parses the output from the grep subprocess child
def grep_parse(grep_string, start_dir):
	global grep_temp

	# Split and write to temporary file so SeqIO can handle
	grep_split = grep_string.split("\n")
	with open(start_dir + "/" + grep_temp, "w") as fh:
		for line in grep_split:
			if line != "--":
				fh.write(line+"\n")
		
	all_recs = []
	for seq_rec in SeqIO.parse(start_dir + "/" + grep_temp, "fastq"):
		all_recs.append(seq_rec)

	# Return all sequence records from the original grep_string
	return all_recs

# frame_eval function
# > moves frame by frame to determine the most likely starting frame
# > returns the index of the starting base of the most likely starting frame, or -1 if no
# > sutiable frame could be found.
def frame_eval(multi_seg, full_seq, seq_motif, max_mismatch):
	frame_idx = 0
	frame_shift = []

	# Primary Layer:
	# Generate frame array
	while frame_idx != len(seq_motif):
		frame_seg = multi_seg[frame_idx:frame_idx+len(seq_motif)]
		#print(frame_seg)
		mismatch=0

		for idx, base in enumerate(frame_seg):
			if base != seq_motif[idx]:
				mismatch += 1

		frame_shift.append(mismatch)
		frame_idx += 1
	
	#frame_shift = [5, 3, 3, 6, 5, 4]
	#print(frame_shift)	
	# Determine best matches
	best_mismatch_array = []
	cur_best_mismatch = max_mismatch
	for idx, entry in enumerate(frame_shift):
		if entry == cur_best_mismatch:
			best_mismatch_array.append(idx)
		elif entry < cur_best_mismatch:
			best_mismatch_array.clear()
			best_mismatch_array.append(idx)
			cur_best_mismatch = entry
	#print(best_mismatch_array)	

	# Primary Layer Check:
	# if triggered, none of the best matches in frame_shift were less
	# than the set max mismatch
	if best_mismatch_array == []:
		return(-1)

	# Secondary Layer:
	# Examine best matches sister frames
	final_frames = []
	for frame_idx in best_mismatch_array:
		sister_frame = full_seq[frame_idx+len(seq_motif):frame_idx+(2*len(seq_motif))]
		#print(sister_frame)
		#print(full_seq)
		mismatch=0

		for idx, base in enumerate(sister_frame):
			if base != seq_motif[idx]:
				mismatch+=1

		if mismatch <= max_mismatch:
			final_frames.append(frame_idx)
			
	#print(final_frames)

	# Secondary Layer Check:
	# if triggered, none of sister frames scored less or at the set max mismatch
	if final_frames == []:

		# Temporary bugfix, current algorithm discounts valid 1 rep sequences,
		# this is a workaround but should be replaced with something better as it
		# might cause other issues
		if len(best_mismatch_array) == 1:
			return(best_mismatch_array[0])
		elif len(best_mismatch_array) > 1:
			print("Encountered non-consistent frame array error")
			return(-1)
			
		return(-1)

	# Should be only one final frame, if not, throw an error and send to nogroup
	# This is a case where the frame evaluation failed
	if len(final_frames) > 1:
		print("Multiple matches detected: Frame evaluation failed!")
		return(-1)

	# Else, should only be one viable frame, return this.
	return(final_frames[0])

# prime_seq_eval function
# > 
def prime_seq_eval(seq_rec, seq_motif, max_mismatch):
	global frame_multi
	
	seg_len = len(seq_motif)
	seq_ptr = 0
	grouping = 0

	#tar_seg = seq_rec.seq[seq_ptr:seq_ptr+seg_len]
	#print(tar_seg)
	#temp = seq_rec.seq[seq_ptr+seg_len:seq_ptr+2*seg_len]
	#print(temp)

	# Determine prefix existance
	init_seg = seq_rec.seq[seq_ptr:seq_ptr+(frame_multi*seg_len)]
	mismatch = 0
	
	# Determine best starting frame
	seq_ptr = frame_eval(init_seg, seq_rec.seq, seq_motif, max_mismatch)
	if seq_ptr == -1:
		return grouping
	

	# Plausible starting frame found
	while seq_ptr != len(seq_rec.seq) and len(seq_rec.seq) - seq_ptr >= seg_len:
		tar_seg = seq_rec.seq[seq_ptr:seq_ptr+seg_len]
		mismatch=0
		
		for idx, base in enumerate(tar_seg):
			if base != seq_motif[idx]:
				mismatch += 1
			
		if mismatch > max_mismatch:
			return grouping
		else:
			grouping += 1
			seq_ptr += seg_len

	return grouping

# group_eval function
# > calls sequence evaluation function(s), generates grouping matrix and "nogroup" list
# > element 0 of grouping mtx will always be nogroup, and element 1 is group 1, etc...
def group_eval(seq_recs, seq_motif, grp_count, mismatches):
	global start_grp # Currently unused

	# initialize grouping matrix	
	group_mtx = []
	for cnt in range(grp_count+1):
		group_mtx.append([])

	# Evaluate each sequence, return its grouping
	for seq_rec in seq_recs:
		#print(seq_rec.id)
		grouping = prime_seq_eval(seq_rec, seq_motif, mismatches)
		#print(grouping)
		

# file_eval function
# > file evaluation loop (iterator function)
def file_eval(seq_file, seq_motif, group_count, mismatches, start_dir):
	# Ignores files without fastq suffix
	if seq_file.split(".")[-1] == "fastq":
		# Grep out sequences containing motif
		cmd = "grep -B 1 -A 2 "+seq_motif+" "+seq_file
		ret = subprocess.run(["grep", "-B", "1", "-A", "2", seq_motif, seq_file], stdout=subprocess.PIPE, universal_newlines=True)
		if ret.returncode != 0:
			print("No matches found:", seq_file)
			return

		# Parse grep output
		target_recs = grep_parse(ret.stdout, start_dir)
		group_eval(target_recs, seq_motif, group_count, mismatches)
		#print(len(target_recs))

# main function
# > parses parameters and iterates over files
def main():
	# Set globals
	global default_mistmatches
	global default_groupings

	# Argument parser
	parser = argparse.ArgumentParser()
	parser.add_argument('fastq_directory', help='fastq_directory HELP')
	parser.add_argument('motif', help='motif HELP')
	parser.add_argument('-g', '--number-of-groupings', default=default_groupings, help='number-of-groupings HELP', type=int)
	parser.add_argument('-m', '--mismatches', default=default_mismatches, help='mistmatches HELP', type=int)
	args = parser.parse_args()

	# Set parameters
	fastq_dir = args.fastq_directory
	seq_motif = args.motif
	group_count = args.number_of_groupings
	mismatches = args.mismatches
	start_dir = os.getcwd()

	# Catches problems with directory, not the files it contains
	try:
		os.chdir(fastq_dir)
		for seq_file in os.listdir(os.getcwd()):
			file_eval(seq_file, seq_motif, group_count, mismatches, start_dir)
	except Exception as err:
		crit_error(err)

	return 0

main()
