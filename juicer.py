#!/usr/bin/env python3

# -------------------------------------------------------------------

# File: juicer.py
# Author: Harrison Inocencio
# Date: 06-01-18
# Purpose: Extract sequences from a read library according to if they contain a specific motif, 
# then sort according to the number of sequence repeats. Inputs a directory containing fastq 
# reads, and outputs n number of files according to how many repeat groupings are desired

# Usage: try ./juicer.py -h

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
# 5. 
# 6. sequences less than the 2n base pairs, with n being the length of the motif, will not
#    evaluate correctly since the current function needs at least that much for its frame
#	 evaluation 
# 7. Above examples use a small motif of 3 base pairs, but this program is untested for small
#	 motifs (<6 bp) and probably would not function well
# 8. Currently disregards reads with repeat counts higher than the user's group count variable
# 9. Does not handle single/multi base insertions or deletions, since it uses a frame system
#	 to count the repeats

# TODO:
# 1. Detail Notes better/Complete descriptions
# 2. Write parser help descriptions
# 3. Error check all input parameters
# 4. Add option for fasta, or have a specific optional for fasta and fastq
# 5. +++May need to add clean-up function for errors
# 6. Add option for retaining all generated files not just outputs
# 7. 
# 8. Find better solution for grep parsing besides writing a file...
# 9. Possibly add start_grp to parameters, (implement a range), some major changes may be
#	 required for this 
# 10. 
# 11. Complex function to determine starting frame, show to expert
# 12. Make option for a log/dump file that logs the decisions made for each read
# 13. Compartamentalize the mismatch code
# 14. Add constraints to the motif parameter >= 3 etc.
# 15. Frame evaluation reaches a lot, needs more error checking (long reads, etc)
# 16. 
# 17. Improve frame evaluation error
# 18. Check reads before hand, throw out ones shorter than 2n where n is the length of the motif
# 19. Need a better error handling system, root try should handle more exceptions/ should be more
#	  individual try/except blocks
# 20.
# 21. 
# 22. might want init to create temp files as well, and maybe check variables.
#	  essentially make sure that the program can do everything it needs to before it starts
#	  actually running. Also check for grep install just in case
# 23. +++Consider merging file_eval into main, passing too many variables, bad practice...
# 24. Add compression support
# 25. Currently throws out repeats with more reps than the current grouping, modify so they
#	  write out to their own file
# 26. +++Finish mop, attach to unrecoverable excepts
# 27. Clean up code, make consistent
# 28. Effectively utilize oop
# 29. +++bug with long output paths
# 30. Use repeat multiplier
# 31. Begin implementing a single shifting system, with scoring system and weight based on blast

# NOTES:
# 1. Possibly base the last piece of frame evaluation off of the sum of the target frame and its
#	 sister frame
# 2. Should handle single base insertions/deletions in the future?
# 3. Program currently ignores prefixes, should this change?
# 4. Does grep sequence need a multiplyer?


# -------------------------------------------------------------------

#Globals
default_mismatches = 1
default_groupings = 30
default_output = "juicer_out"
start_grp = 1 # Currently unused
frame_multi = 2	# Probably dont want to change this...
zfill_amount = 4 # Amount of zero padding in file names
min_motif_len = 5 # Minimum allowed motif length
grep_temp = ".grep_temp_juicer"
out_prefix = "juicer_group"

import os
import shutil
import argparse
import subprocess
from Bio import SeqIO

# init function
# > creates initial environment and startup files
# > returns an array of file names, with nogroup as ele 0, 1 rep as ele 1 etc...
def init(group_count, output_dir):
	global out_prefix
	global zfill_amount

	original_dir = os.getcwd()
	# Init output dir, remove if already exists
	try:
		os.mkdir(output_dir)
	except FileExistsError:
		shutil.rmtree(output_dir)
		os.mkdir(output_dir)
	os.chdir(output_dir)	

	# Init output files/array
	name_array = []
	for group in range(group_count+1):
		cur_name = out_prefix+str(group).zfill(zfill_amount)+".fastq" # LINE CHANGE
		open(cur_name, 'w').close()
		name_array.append(output_dir+"/"+cur_name)

	os.chdir(original_dir)
	
	return name_array

# mop function
# > clean up function, called when program quit unexpectedly
def mop(start_dir, output_dir):
	global grep_temp

	# Attempt to remove any generated files
	os.chdir(start_dir)
	try:
		os.remove(grep_temp)
	except Exception:
		pass
	try:
		shutil.rmtree(output_dir)
	except Exception:
		pass
	print()
	
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

	# Remove temporary file
	os.remove(start_dir + "/" + grep_temp)

	# Return all sequence records from the original grep_string
	return all_recs

# parameter_check function
# > checks input parameters and imposes set limits
def parameter_check(seq_motif, group_count):
	global min_motif_len
	
	if group_count == 0:
		raise RuntimeError("group count set to zero!")
	elif len(seq_motif) <= min_motif_len:
		raise RuntimeError("motif too short")

	return 0

# frame_eval function
# > moves frame by frame to determine the most likely starting frame
# > returns the index of the starting base of the most likely starting frame, or -1 if no
# > sutiable frame could be found.
# > Additional: return of -2 indicates that the seq encountered a non-consistent frame array 
# > error. This indicates that two or more likely frames were found in which none had 
# > sister frames <= the max set mismatch
# > return of -3 indicates the multiple final frames were encountered. In this situation
# > the frame evaluation is unable to determine a correct starting frame, and places the sequence
# > in "nogroup"

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
	# Examine best match's sister frames
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
	# if triggered, none of the sister frames scored less or at the set max mismatch
	if final_frames == []:

		# Temporary bugfix, current algorithm discounts valid 1 rep sequences,
		# this is a workaround but should be replaced with something better as it
		# might cause other issues
		if len(best_mismatch_array) == 1:
			return(best_mismatch_array[0])
		elif len(best_mismatch_array) > 1:
			#print("Encountered non-consistent frame array error")
			return(-2)
			
		return(-1)

	# Should be only one final frame, if not, throw an error and send to nogroup
	# This is a case where frame evaluation fails
	if len(final_frames) > 1:
		#print("Multiple matches detected: Frame evaluation failed!")
		return(-3)

	# Else, should only be one viable frame, return this.
	return(final_frames[0])

# prime_seq_eval function
# > Evaluates sequences using frame shifts
def prime_seq_eval(seq_rec, seq_motif, max_mismatch):
	global frame_multi
	
	seg_len = len(seq_motif)
	seq_ptr = 0
	grouping = 0

	#tar_seg = seq_rec.seq[seq_ptr:seq_ptr+seg_len]
	#print(tar_seg)
	#temp = seq_rec.seq[seq_ptr+seg_len:seq_ptr+2*seg_len]
	#print(temp)

	# Extract initial segment for frame_eval
	init_seg = seq_rec.seq[seq_ptr:seq_ptr+(frame_multi*seg_len)]
	mismatch = 0
	
	# Determine best starting frame
	seq_ptr = frame_eval(init_seg, seq_rec.seq, seq_motif, max_mismatch)
	
	# Layer check triggered
	if seq_ptr == -1:
		return grouping
	# Non-consistent frame array error
	elif seq_ptr == -2:
		return grouping
	# Multiple matches error
	elif seq_ptr == -3:
		return grouping

	# Plausible starting frame found, begin counting repeats
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
# > calls sequence evaluation function(s), generates grouping matrix
# > element 0 of grouping mtx will always be nogroup, and element 1 is group 1, etc...
def group_eval(seq_recs, seq_motif, grp_count, mismatches):
	global start_grp # Currently unused

	# initialize grouping matrix	
	group_mtx = []
	for cnt in range(grp_count+1):
		group_mtx.append([])

	# Evaluate each sequence, return its grouping and insert into group mtx
	for seq_rec in seq_recs:
		#print(seq_rec.id)
		grouping = prime_seq_eval(seq_rec, seq_motif, mismatches)
		#print(grouping)

		# Disregard sequence if its groupings is higher than the set one
		if grouping <= grp_count:
			group_mtx[grouping].append(seq_rec)
		
	return group_mtx

# fasta_write function
# > writes out a group_mtx to the initialized output files in fasta format.
def fasta_write(group_mtx, start_dir, name_array):
	original_dir = os.getcwd()
	os.chdir(start_dir)

	#print(len(group_mtx))
	for group_num in range(len(group_mtx)):
		if group_mtx[group_num] != []:
			for entry in group_mtx[group_num]:
				with open(name_array[group_num], 'a') as fh:
					fh.write(">%s\n%s\n" % (entry.description, entry.seq))

	os.chdir(original_dir)

	return 0
# fastq_write function
# > writes out a group_mtx to the initialized output files in fastq format.
def fastq_write(group_mtx, start_dir, name_array):
	original_dir = os.getcwd()
	os.chdir(start_dir)

	for group_num in range(len(group_mtx)):
		if group_mtx[group_num] != []:
			for entry in group_mtx[group_num]:
				with open(name_array[group_num], 'a') as fh:
					#print(entry.format("fastq"))
					fh.write(entry.format("fastq"))

	os.chdir(original_dir)

	return 0

# main function
# > parses parameters and iterates over files
def main():
	# Set globals
	global default_mistmatches
	global default_groupings
	global default_output

	# Argument parser
	parser = argparse.ArgumentParser()
	parser.add_argument('fastq_directory', help='Path to directory containing fastq files to be juiced. Searches each file for set motif and groups them based on the number of times this motif repeats from the start of the sequence')
	parser.add_argument('motif', help='Nucleic motif to be used as the repeat base')
	parser.add_argument('-g', '--number-of-groupings', default=default_groupings, help='Specifies the number of groupings to be used (This also codes for the maximum repeat count)(Default = 30)', type=int)
	parser.add_argument('-m', '--mismatches', default=default_mismatches, help='Number of mismatches to be allowed in each instance of the motif (Default = 1)', type=int)
	parser.add_argument('-o', '--output-name', default=default_output, help='Name of the output directory for juicer (Default = juicer_out)')
	args = parser.parse_args()

	# Set parameters
	fastq_dir = args.fastq_directory
	seq_motif = args.motif
	group_count = args.number_of_groupings
	mismatches = args.mismatches
	output_name = args.output_name
	start_dir = os.getcwd()

	ret = parameter_check(seq_motif, group_count)

	# Root try statement
	try:
		parameter_check(seq_motif, group_count)
		name_array = init(group_count, output_name)
		os.chdir(fastq_dir)
	
		# File evaluation
		for seq_file in os.listdir(os.getcwd()):
			
			#OLD FUNCTION CALL
			#file_eval(seq_file, seq_motif, group_count, mismatches, start_dir, name_array)
			
			# Ignores files without fastq suffix
			if seq_file.split(".")[-1] == "fastq":

			# Grep out sequences containing motif
				print("Evaluating sequence file:", seq_file)

				cmd = "grep -B 1 -A 2 "+seq_motif+" "+seq_file
				ret = subprocess.run(["grep", "-B", "1", "-A", "2", seq_motif, seq_file], stdout=subprocess.PIPE, universal_newlines=True)
				if ret.returncode != 0:
					print("No matches found:", seq_file)
				else:
					# Parse grep output
					target_recs = grep_parse(ret.stdout, start_dir)

					# Generate Group matrix
					group_mtx = group_eval(target_recs, seq_motif, group_count, mismatches)
	
					# Write out file results
					fastq_write(group_mtx, start_dir, name_array) # LINE CHANGE

					#print(len(target_recs))
					#print(group_mtx)


		# Remove empty output files
		os.chdir(start_dir)
		for name in name_array:
			stats = os.stat(name)
			if stats.st_size == 0:
				os.remove(name)
			
	except KeyboardInterrupt:
		mop(start_dir, output_name)
		exit(1)
	except RuntimeError:
		exit(1)
	except Exception as err:
		mop(start_dir, output_name)
		crit_error(err)

	return 0

main()
