### import modules
from argparse import ArgumentTypeError
import argparse
import os
from Bio import SeqIO
import re
import numpy as np
import random
import khmer
import warnings
import pandas as pd




import math
import subprocess
import screed
import bisect
import inspect
import copy
import gc
import statistics


### basic functions
def parsenumlist(k_sizes_str: str):
	"""
	Parses a string like 10-21-1 and turn it into a list like [10, 11, 12,...,21]
	:param k_sizes_str: the <start>-<end>-<increment> string
	:type k_sizes_str: str
	:return: list of k-mer sizes
	:rtype: list
	"""
	m = re.match(r'(\d+)(?:-(\d+))?(?:-(\d+))?$', k_sizes_str)
	# ^ (or use .split('-'). anyway you like.)
	if not m:
		raise ArgumentTypeError(
			"'" + k_sizes_str + "' is not a range of number. Expected forms like '1-5' or '2' or '10-15-2'.")
	start = int(m.group(1))
	end = int(m.group(2))
	if m.group(3):
		increment = int(m.group(3))
	else:
		increment = 1
	return list(range(start, end + 1, increment))

def build_kmer_dict_with_start_pos(seq, ksize: int):
	"""
	Build a dict with key=kmer, value=kmer start pos (include multiple hits)
	:param seq:
	:param k:
	:return: a dict
	"""
	temp=dict()
	for i in range(len(seq) - ksize + 1):
		kmer = seq[i:i + ksize]
		if kmer in temp:
			temp[kmer].append(i)
		else:
			temp[kmer] = [i]
	return temp

def dup_count_for_kmer_dict(input_dict, item_name):
	"""
	Count dup kmers in a given kmer dict: kmer:[pos]
	"""
	all_value_len = [len(x) for x in input_dict.values()]
	multi_hit = [x>1 for x in all_value_len]
	print("Dup count of %s \nTotal distinct: %s\nDuplicated items: %s \n\n" %(item_name, len(all_value_len), sum(multi_hit)))
	
def construct_strobemer(seq, n, l, wmin, wmax, prime=6229, x=3):
	"""
	Generate minstrobe for given seq
	:param seq:
	:param l: strobe length
	:param n: strobe number
	:param wmin: window gap min
	:param wmax: window gap max
	:param prime: a prime # to mod, don't need a large one
	#:param: x: partition parameter, use 3
	:return: dict with strobemer:[ [m1...mn], ...  ]
	"""
	minstrobe_dict = dict()
	randstrobe_dict = dict()
	hybridstrobe_dict = dict()
	# loop through index
	for i in range(len(seq) - n * l + 1):
		# max window size: only activate when approaching the tail of seq, o.w. is w_max
		wu = int(min(wmax, (len(seq) - i - l+1) / (n - 1)))  # in case of wu -wl < 1, use int
		# min window size: make sure no overlap strobes: lenth l < window min
		wl = max(wmin - (wmax - wu), l)
		# manually stop (may generate empty records in the j loop)
		if wu < wl:
			break
		elif wu == wl:
			# start from here, all remaining strobemers are regular kmer (may have wu - wl < 1 -> several kmers at end)
			fix_seq = seq[i:i + n * l]
			fix_pos = [i]
			for j in range(1, n):
				fix_pos.append(i + j * l)
			# add to min_dict
			if fix_seq in minstrobe_dict:
				minstrobe_dict[fix_seq].append(fix_pos)
			else:
				minstrobe_dict[fix_seq] = [fix_pos]
			# also add to rand_dict
			if fix_seq in randstrobe_dict:
				randstrobe_dict[fix_seq].append(fix_pos)
			else:
				randstrobe_dict[fix_seq] = [fix_pos]
			# end this loop here
			continue
		
		# minstrobe: 1st strobe
		m1 = seq[i:i + l]
		min_seq = m1
		min_pos = [i]
		# minstrobe loop
		for j in range(2, n + 1):
			# this will generate empty list when wu <= wl
			select_window = list(range(i + wl + (j - 2) * wu, i + (j - 1) * wu, 1))
			current_pos = 0
			current_hash = prime + 10  # init loop
			# find argmin in current window
			for pos in select_window:
				new_hash = khmer.hash_no_rc_murmur3(seq[pos:pos + l]) % prime
				if new_hash < current_hash:
					current_hash = new_hash
					current_pos = pos
			# record
			min_pos.append(current_pos)
			min_seq += seq[current_pos:current_pos + l]
		# add to the dict
		if min_seq in minstrobe_dict:
			# add strobe starting indices [x,x,x] into a list, so use append
			minstrobe_dict[min_seq].append(min_pos)
		else:
			minstrobe_dict[min_seq] = [min_pos]
		
		# ranstrobe: 1st strobe
		rand_seq = m1
		rand_pos = [i]
		last_strobe = m1
		# randstrobe loop
		for j in range(2, n + 1):
			# this will generate empty list when wu <= wl
			select_window = list(range(i + wl + (j - 2) * wu, i + (j - 1) * wu, 1))
			current_pos = 0
			current_hash = prime + 10  # init loop
			# find argmin in current window
			for pos in select_window:
				new_hash = (khmer.hash_no_rc_murmur3(seq[pos:pos + l]) + khmer.hash_no_rc_murmur3(last_strobe)) % prime
				if new_hash < current_hash:
					current_hash = new_hash
					current_pos = pos
			# update last strobe
			last_strobe = seq[current_pos:current_pos + l]
			# record
			rand_pos.append(current_pos)
			rand_seq += seq[current_pos:current_pos + l]
		# add to the dict
		if rand_seq in randstrobe_dict:
			randstrobe_dict[rand_seq].append(rand_pos)
		else:
			randstrobe_dict[rand_seq] = [rand_pos]
			
		# hybridstrobe: 1st strobe
		hy_seq = m1
		hy_pos = [i]
		last_strobe = m1
		# hybridstrobe loop
		for j in range(2, n + 1):
			# this will generate empty list when wu <= wl
			select_window = list(range(i + wl + (j - 2) * wu, i + (j - 1) * wu, 1))
			# additional step to subset selection window by x
			slice_index = khmer.hash_no_rc_murmur3(last_strobe) % x
			if len(select_window) < x:
				# no need to slice with few elements
				sub_window = select_window
			else:
				sub_window = list(np.array_split(select_window, x)[slice_index])
			# now we get the r-th subwindow to play
			current_pos = 0
			current_hash = prime + 10  # init loop
			# find argmin in current window
			for pos in sub_window:
				new_hash = (khmer.hash_no_rc_murmur3(seq[pos:pos + l]) + khmer.hash_no_rc_murmur3(last_strobe)) % prime
				if new_hash < current_hash:
					current_hash = new_hash
					current_pos = pos
			# update last strobe
			last_strobe = seq[current_pos:current_pos + l]
			# record
			hy_pos.append(current_pos)
			hy_seq += seq[current_pos:current_pos + l]
		# add to the dict
		if hy_seq in hybridstrobe_dict:
			hybridstrobe_dict[hy_seq].append(hy_pos)
		else:
			hybridstrobe_dict[hy_seq] = [hy_pos]
		
		
	
	# return results
	return minstrobe_dict, randstrobe_dict, hybridstrobe_dict




### workflow functions
def s1_generate_string_list_from_genome(input_genome):
	"""
	Transfer input genome to a string
	"""
	# store fasta to a list (in case multiple records line)
	fasta_sequences = SeqIO.parse(input_genome, 'fasta')
	seq_list = []
	for fasta in fasta_sequences:
		seq_list.append(str(fasta.seq))
	# merge every record to a string
	if len(seq_list) == 0:
		raise Exception("couldn't read fasta from input file, please double check!")
	else:
		print("There are " + str(len(seq_list)) + " record(s).")
		return seq_list

def s2_generate_mutated_strings(input_string, p):
	"""
	Generate a mutated string by simple mutation model
	"""
	# simple mutation event by Unif(0,1)
	string_length = len(input_string)
	# determine the location to be mutated
	random_unif1 = list(np.random.uniform(size=string_length))
	mutation_id = [x < p for x in random_unif1]  # use collections.Counter(mutation_id) to check number of elements
	# mutate string
	out_string = ""
	counter = 0
	while counter < string_length:
		if not mutation_id[counter]:
			# no mutation
			out_string += input_string[counter]
		else:
			# choose from insertion / deletion / substitution of 1/3 prob
			mut_type = random.choice(['ins', 'del', 'sub'])
			if mut_type == 'ins':
				# inserted letter
				out_string += random.choice('ACGT')
				# original letter
				out_string += input_string[counter]
			elif mut_type == 'del':
				out_string += "-"
			else:
				dest_space = "ACGT".replace(input_string[counter], "")
				out_string += random.choice(dest_space)
		counter += 1
	# return
	out_string = out_string.replace("-", "") # trail the deletion letter
	return out_string

def s3_generate_all_types_of_kmers(input_string, k_value, strobe_size=5, wmin=10, wmax=20):
	"""
	Generate regular k-mer, spaced k-mer, and 3 types of strobmer
	:return: a dict with keys ['kmer', 'spaced_sparse', 'spaced_dense', 'minstrobe', 'randstrobe', 'hybridstrobe'], values are dicts of kmers with starting point
	"""
	out_dict = dict()
	### regular kmer
	out_dict['kmer'] = build_kmer_dict_with_start_pos(input_string, k_value)
	
	### spaced_kmer sparse: L=3k, 1/3 occupied
	temp_dict = build_kmer_dict_with_start_pos(input_string, k_value*3)
	# this is a 3k-mer, then we modify it to a kmer with fixed pattern (uniform)
	for k3mer in list(temp_dict):
		# spaced_kmer sparse
		new_kmer = k3mer[::3]
		# refresh keys, this may introduce dup! So need to use "EXTEND", but in regular kmer we don't need it
		if new_kmer in temp_dict:
			temp_dict[new_kmer].extend(temp_dict.pop(k3mer))
		else:
			temp_dict[new_kmer] = temp_dict.pop(k3mer)
	out_dict['spaced_sparse'] = temp_dict
	del temp_dict, k3mer
	
	### spaced_kmer dense: L=1.5k, 2/3 occupied
	temp_dict = build_kmer_dict_with_start_pos(input_string, int(k_value * 1.5))
	for k1_5mer in list(temp_dict):
		# spaced_kmer sparse
		kmer_2_list = list(k1_5mer)
		del kmer_2_list[2::3] # delete every 3rd letter
		new_kmer = "".join(kmer_2_list)
		# refresh keys, this may introduce dup!
		if new_kmer in temp_dict:
			temp_dict[new_kmer].extend(temp_dict.pop(k1_5mer))
		else:
			temp_dict[new_kmer] = temp_dict.pop(k1_5mer)
	out_dict['spaced_dense'] = temp_dict
	del temp_dict, k1_5mer
	
	strobe_num = k_value // strobe_size
	if k_value / strobe_size != int(k_value / strobe_size):
		warnings.warn("Notice: k_value can't be divided by strobe size, the length may NOT match between kmer and strobemer")
		
	### minstrobe, randstrobe and hybridstrobe
	out_dict['minstrobe'], out_dict['randstrobe'], out_dict['hybridstrobe'] = \
		construct_strobemer(input_string, n=strobe_num , l=strobe_size, wmin=wmin, wmax=wmax, prime=6229, x=3)

	### dup count
	for item in out_dict:
		dup_count_for_kmer_dict(out_dict[item], item)
	
	### output
	return out_dict

def s4_get_eval_matrix_between_kmer_dict_from_2_strings(dict1, dict2):
	"""
	Get py-version SIM-R table for 1 pair, may repeat for many times
	:param dict1: kmer dict for str1
	:param dict2: kmer dict for mutated str1
	:return: a list [category, parameter, m, sc, mc, E]
	"""
	print("hello")
	




# test code
def local_test_variable():
	"""
	setup var for local test
	"""
	print("Loading local variables")
	genome_file = '/Users/shaopeng/Library/CloudStorage/OneDrive-ThePennsylvaniaStateUniversity/PSU_academic/Koslicki_lab/research/strobemer_related/input_genomes/GCA_000014345.1_ASM1434v1_genomic.fna'
	sample_size = 10
	rev_comp = True
	k_sizes = [15, 20, 25]
	p_values = [0.01, 0.05, 0.1]




### run analysis
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Strobemer simulation",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', help="Path to input genomes")
	parser.add_argument('-k', '--k_range', type=str, help="k-mer range", default="15-50-5")
	parser.add_argument('-n', '--sample_size', type=int, help="Number of technical replicates", default=10)
	parser.add_argument('-s', '--strobe_size', type=int, help="Length of single strobe", default=5)
	parser.add_argument('-p', '--p_range', type=str, help="mutation rate range",
	                    default="0.01,0.03,0.05,0.08,0.1")
	parser.add_argument('-r', '--rev_comp', type=str, help="Whether to keep the reverse complementary", default="False")
	parser.add_argument('-w', '--window_min', type=int, help="Min window size", default=10)
	parser.add_argument('-e', '--window_max', type=int, help="Max window size", default=20)
	
	########################################## read parameters
	args = parser.parse_args()
	
	if not os.path.exists(args.genome):
		raise Exception("Input file %s does not exist." % args.genome)
	genome_file = os.path.abspath(args.genome)
	print("Input genome file is " + genome_file)
	
	sample_size = args.sample_size
	print("Permutation sample size is " + str(sample_size))
	
	rev_comp = args.rev_comp
	rev_comp = rev_comp == 'True'
	if rev_comp:
		print("Use canonical kmers in analysis!!!")
	
	k_range = args.k_range
	k_sizes = parsenumlist(k_range)
	print("The range of k values are: ")
	print(k_sizes)
	
	p_range = args.p_range
	p_values = [float(x) for x in p_range.split(",")]
	print("The range of mutation rate p values are: ")
	print(p_values)
	
	
	
	
	########################################## Step1: transfer genome to a sequence list (as it may have multiple chrs)
	### transfer genome file to string list
	not_unknown = re.compile('[nN]')
	
	raw_seq = s1_generate_string_list_from_genome(genome_file)
	cleaned_seq_list=[]
	for fragment in raw_seq:
		fragment = fragment.upper()
		frag_split_onlyACTG = not_unknown.split(fragment)
		cleaned_seq_list.extend([x for x in frag_split_onlyACTG if len(x)>0])
	
	### for this test, we only need 1 string not the multiple records per chr, so manually concatenate it
	# paste seq
	out_seq = "".join(cleaned_seq_list)[:10000]
	print("The seq length is "+str(len(out_seq)))



	########################################## Step2~4: mutate, generate k-mers, compare statistics
	seq0 = out_seq
	seq1 = s2_generate_mutated_strings(seq0, 0.1)
	
	dict0 = s3_generate_all_types_of_kmers(seq0, k_value=30)
	dict1 = s3_generate_all_types_of_kmers(seq1, k_value=30)
	
	
	
	




