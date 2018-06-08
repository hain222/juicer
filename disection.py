# File: disection.py
# Date: 06-08-18
# Author: Harrison Inocencio

# Description: Generates the color coded interface used when the disection option is toggled

import colorama

# disection function
# >
def dprint(start_idx, record, num_reps, motif_len):
	colorama.init()

	print("%s = Group %s" % (record.description, num_reps))
	#print(record.seq)
	
	for idx in range(start_idx):
		print(record.seq[idx], end='')

	for idx in range(start_idx, num_reps*motif_len):
		print(colorama.Back.GREEN + record.seq[idx], end='')
	print(colorama.Style.RESET_ALL, end='')
	
	for idx in range(num_reps*motif_len, num_reps*motif_len+motif_len):
		print(colorama.Back.RED + record.seq[idx], end='')
	print(colorama.Style.RESET_ALL, end='')

	for idx in range(num_reps*motif_len+motif_len, len(record.seq)):
		print(record.seq[idx], end='')
	print()

	input("Press a key to continue...")
