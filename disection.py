# File: disection.py
# Date: 06-08-18
# Author: Harrison Inocencio

# Description: Generates the color coded interface used when the disection option is toggled

import colorama

# disection function
# >
def dprint(start_idx, record, num_reps, motif_len):
	colorama.init()

	# Print id
	print("%s = Group %s" % (record.description, num_reps))
	#print(record.seq)
	
	shifted_rep_len = num_reps*motif_len+start_idx

	# Print prefix
	for idx in range(start_idx):
		print(record.seq[idx], end='')

	# Print motif repeats
	for idx in range(start_idx, shifted_rep_len):
		print(colorama.Back.GREEN + record.seq[idx], end='')
	print(colorama.Style.RESET_ALL, end='')
	
	# Print next frame (frame that failed) if sequence is long enough
	if shifted_rep_len+motif_len <= len(record.seq):
		for idx in range(shifted_rep_len, shifted_rep_len+motif_len):
			print(colorama.Back.RED + record.seq[idx], end='')
		print(colorama.Style.RESET_ALL, end='')

		for idx in range(shifted_rep_len+motif_len, len(record.seq)):
			print(record.seq[idx], end='')

	# Print the end of sequence if sequence not long enough
	else:
		for idx in range(shifted_rep_len, len(record.seq)):
			print(record.seq[idx], end='')
		
	print()

	input("Press a key to continue...")
