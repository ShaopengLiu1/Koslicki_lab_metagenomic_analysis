### import modules
import multiprocessing
from multiprocessing import Pool
import argparse
import os
import khmer
from itertools import *
import numpy as np
import pandas as pd
import re
from argparse import ArgumentTypeError
import subprocess
import screed
import bisect
import inspect
import copy
import gc
from Bio import SeqIO
import collections
import random


### regular functions
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

def estimate_genome_size(sequence, ksize):
	"""
	fast kmer cardinality estimate (copied from CMash, but didn't find it from MH)
	:param genome:
	:param ksize:
	:return: estimated kmer cardinality of a genomic file on a given k_size
	"""
	hll = khmer.HLLCounter(0.01, ksize)
	hll.consume_string(sequence)
	return hll.estimate_cardinality()

def containment_to_jaccard(C_est, cardi_A, cardi_B):
	"""
	transfer CI to JI by their definition
	:param C_est: CI(A, B) = |A \cap B|  /  |A|
	:return: JI(A, B)
	"""
	return C_est * cardi_A  / (cardi_A + cardi_B - C_est * cardi_A)

def kmers(seq, ksize):
	for i in range(len(seq) - ksize + 1):
		yield seq[i:i + ksize]


### Workflow related
def s1_1_generate_string_from_genome(input_genome, rev_comp):
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
		print("There are "+str(len(seq_list))+" record(s).")
		out_seq=""
		for i in seq_list:
			out_seq = out_seq+i
	# return
	return out_seq

def s2_1_generate_mutated_strings(input_string, p):
	"""
	Generate a mutated string by simple mutation model
	"""
	# simple mutation event by Unif(0,1)
	string_length = len(input_string)
	random_unif1 = list(np.random.uniform(size=string_length))
	mutation_id = [x<p for x in random_unif1]    # use collections.Counter(mutation_id) to check number of elements
	# simple mutation destination
	out_string=""
	counter=0
	while counter < string_length:
		if not mutation_id[counter]:
			# no mutation
			out_string += input_string[counter]
		else:
			dest_space="ACGT".replace(input_string[counter], "")
			out_string += random.choice(dest_space)
		counter += 1
	# return
	return out_string
	
def s2_2_check_mutated_ratio(str1, str2):
	"""
	Personal usage: validate the edit distance between 2 equal length strings
	"""
	if len(str1) != len(str2):
		print("2 strings are not equal length, stop!")
		return
	else:
		counter=0
		mismatch=0
		while counter < len(str1):
			if str1[counter] != str2[counter]:
				mismatch += 1
			counter += 1
	
	mis_ratio = mismatch / len(str1)
	return mis_ratio
	
def s2_3_simulation_for_given_p(seq, p, n):
	"""
	perform simulation for a given p and sample size n
	"""
	mutated_seq_list = []
	edit_dist_ratio_list = []
	for i in range(n):
		simu_seq = s2_1_generate_mutated_strings(seq, p)
		mutated_seq_list.append(simu_seq)
		# for validation purpose
		simu_mis_ratio = s2_2_check_mutated_ratio(seq, simu_seq)
		edit_dist_ratio_list.append(simu_mis_ratio)
	return mutated_seq_list, edit_dist_ratio_list
	
def s3_1_ci_ji_between_mutated_list(original_seq, mutated_list, kvalue):
	"""
	generate CI/JI array for a given k value
	"""
	original_seq = raw_seq
	mutated_list = p_seq_list
	kvalue = 20
	
	# intermediate item
	cardi_original = estimate_genome_size(original_seq, kvalue)
	cardi_mutated_list = []
	ref_sketch = [""] * 2000
	max_prime = 9999999999971. # similar to CMash, compress hash space
	_mins = [max_prime] * 2000
	
	# use lists to store results
	out_JI = []
	out_CI = []
	
	# generate MinHash sketch for original_seq, size = 2000
	for kmer in kmers(original_seq, kvalue):
		if rev_comp:
			kmer = min(kmer, khmer.reverse_complement(kmer))
		h = khmer.hash_no_rc_murmur3(kmer) % max_prime
		if h >= _mins[-1]:
			continue
		i = bisect.bisect_left(_mins, h)  # find index to insert h
		if _mins[i] == h:
			continue #already found
		else:
			_mins.insert(i, h)
			_mins.pop()
			ref_sketch.insert(i, kmer)
			ref_sketch.pop()
			
	# delete unused cells in sketch (in case of small genome)
	while _mins[-1] == max_prime:
		_mins.pop()
		ref_sketch.pop()
	
	# loop through mutated list and stream kmers
	
	
	
	
	
	







def local_test_variable():
	"""
	setup var for local test
	"""
	print("Loading local variables")
	genome_file='/Users/shaopeng/Desktop/Metagenomic_analysis/3_mini_analysis/20210705_point_mutation_ecoli_permutation_ANI_JI/GCA_002220215.1_ASM222021v1_genomic.fna'
	sample_size=100
	rev_comp=True
	k_sizes=[15, 20]
	p_values=[0.001, 0.005, 0.1]





### run analysis
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="ANI vs JI analysis",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', help="Path to input genomes")
	parser.add_argument('-k', '--k_range', type=str, help="k-mer range", default="15-60-5")
	parser.add_argument('-n', '--sample_size', type=int, help="Number of mutated strings.", default=10000)
	parser.add_argument('-p', '--p_range', type=str, help="mutation rate range", default="0.001,0.005,0.1,0.15,0.5")
	parser.add_argument('-r', '--rev_comp', type=str, help="Whether to keep the reverse complementary", default="True")


	########################################## read parameters
	args = parser.parse_args()

	if not os.path.exists(args.genome):
		raise Exception("Input file %s does not exist." % args.genome)
	genome_file=os.path.abspath(args.genome)

	sample_size=args.sample_size

	rev_comp = rev_comp == 'True'

	k_range=args.k_range
	k_sizes=parsenumlist(k_range)

	p_range=args.p_range
	p_values=[float(x) for x in p_range.split(",") ]



	########################################## Step1: generate n mutated strings
	### 1.1: transfer genome file to string
	raw_seq = s1_1_generate_string_from_genome(genome_file, rev_comp=True)
	raw_seq = raw_seq.upper()
	# remove non ACTG letter (just in case)
	notACTG = re.compile('[^ACTG]')
	seq_split_onlyACTG = notACTG.split(raw_seq)
	if len(seq_split_onlyACTG) > 1:
		raw_seq=""
		for short_seq in seq_split_onlyACTG:
			raw_seq += short_seq
	
	
	# for test purpose, trunc raw seq to 10000
	raw_seq = raw_seq[:10000]




	# Start from here:
	# we'll loop through all p values
	# those codes are for single p test
	# will merge all downstream steps into a function later
	# for mute_p in p_values:
	
	
	
	### Step2: generate mutated m sequences
	mute_p = 0.1 # delete this line in loop
	p_seq_list, p_mis_ratio_list = s2_3_simulation_for_given_p(seq=raw_seq, p=mute_p, n=sample_size)
	
	
	### Step3: generate 1 vs n (a list/array) of JIs and CIs by containment MinHash (requires smaller sample size)
	for k in k_sizes:
		print(k)
		
		
	


	





