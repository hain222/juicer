# File: group.py
# Date:
# Author: Harrison Inocencio

# Purpose: Grouping method 2 for juicer. Marks each feasible motif
#		   occurance, and then uses a set range to decide which gaps
#		   gaps are large enough to bridge and which are too large

import constants

#seq = "CTAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAAACCCAAACCCAAACCCAAACCCGAAGGGTTCCCAAGACGCCTAAACCCGAAGGGTTTAGGATATTATTTCGTTAGATCGGAAGAGAACACGTCTGAACTCCAGTCACCAGTCCCGAAATCGTATGCCGTCTTCTGCTTGAAAAAAAAAAA"
#motif = "CCCTAA"

# det_mismatch function
# > Given two strings of equal size
# > return the number of mismatches between the two strings
def det_mismatch(frame, motif):
	mismatch = 0
	for idx in range(len(frame)):
		if frame[idx] != motif[idx]:
			mismatch +=1

	return mismatch

# hit_print function
# > Print sequence and a second line containing vertical
# > bars that mark each hit in the sequence
def hit_print(seq, hits):
	print(seq)
	for idx in range(len(seq)):
		if idx in hits:
			print("|", end='')
		else:
			print(" ", end='')
	print()

# count function
# > count the number of repeats in the read using method 2
# > output the number of repeats found
def count(hits, m_array, motif):
	one_rep_min = constants.one_rep_min
	two_rep_min = constants.two_rep_min
	#three_rep_min = 15
	max_size = constants.max_size

	distances = []
	for idx in range(len(hits)-1):
		dis = hits[idx+1]-hits[idx]
		distances.append(dis)

	rep_count = 1
	for idx in range(len(distances)):
		if distances[idx] < one_rep_min:
			pass
		elif one_rep_min <= distances[idx] < two_rep_min:
			rep_count+=1
		elif two_rep_min <= distances[idx] < max_size:
			rep_count+=2
		#elif three_rep_min <= distances[idx] < max_size:
			#rep_count +=3
		else:
			break
			#return rep_count

	return rep_count

# group function
# > Give sequence, motif, and some threshold (1)
# > output the determined group for sequence
def group(seq, motif, threshold):
	#threshold = 1 # constant
	m_array = []
	hits = []
	for idx in range(len(seq)):
		frame = seq[idx:idx+len(motif)]
		if len(frame) == len(motif):
			m_array.append(det_mismatch(frame, motif))
	#print(m_array)
	
	
	for idx, mismatch in enumerate(m_array):
		if mismatch <= threshold:
			hits.append(idx)
	if hits[0] >= len(motif):
		return 0
	#print(hits)
	#hit_print(seq, hits)
	rep_count = count(hits, m_array, motif)
	#print(rep_count)

	return rep_count

#group(seq, motif)
