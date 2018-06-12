# File: group.py
# Date:
# Author: Harrison Inocencio

import constants

seq = "CTAACCCTAAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAAACCCAAACCCAAACCCAAACCCGAAGGGTTCCCAAGACGCCTAAACCCGAAGGGTTTAGGATATTATTTCGTTAGATCGGAAGAGAACACGTCTGAACTCCAGTCACCAGTCCCGAAATCGTATGCCGTCTTCTGCTTGAAAAAAAAAAA"
motif = "CCCTAA"

def det_mismatch(frame, motif):
	mismatch = 0
	for idx in range(len(frame)):
		if frame[idx] != motif[idx]:
			mismatch +=1

	return mismatch

#def start_frame(array):

def hit_print(seq, hits):
	print(seq)
	for idx in range(len(seq)):
		if idx in hits:
			print("|", end='')
		else:
			print(" ", end='')
	print()

def gap_eval(gap, m_array, motif):
	gap_seg = m_array[gap[0]:gap[1]]
	print(gap_seg)

	return 0

def count(hits, m_array, motif):
	#len_mod = 1 # constant
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
			return rep_count

		#if len(motif)-len_mod <= distances[idx] <= len(motif)+len_mod:
		#	rep_count += 1
		
		#elif len(motif) + len_mod < distances[idx] > len(motif)+len_mod:
		#	print(hits[idx:idx+2])
		#	rep_count += gap_eval(hits[idx:idx+2], m_array, motif)

	return rep_count

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
