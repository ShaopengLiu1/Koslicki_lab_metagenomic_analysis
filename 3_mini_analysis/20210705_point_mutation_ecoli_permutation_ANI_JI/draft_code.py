### import modules
import argparse
import os
import khmer
import re
from argparse import ArgumentTypeError
import subprocess
import screed
import bisect
import inspect
import copy
import gc
from Bio import SeqIO
import random
import pandas as pd
import numpy as np
import seaborn as sb
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import math
import statistics


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
	"""
	Transfer a seq into kmers
	"""
	for i in range(len(seq) - ksize + 1):
		yield seq[i:i + ksize]

def plot_list_density(input_list, out_name="temp.png"):
	"""
	generate density plot for a list of numbers
	just to double check the random number generator
	"""
	fig, axs = plt.subplots(1, 1, figsize=(12, 12))
	plt.rcParams.update({'font.size': 22})
	sb.histplot(input_list, kde=True, ax=axs)
	axs.set(xlabel="Mutation ratio")
	fig.suptitle("Histplot of mutation rate")
	fig.savefig(out_name, dpi=200)
	plt.close(fig)





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
	
	# intermediate item
	cardi_original = estimate_genome_size(original_seq, kvalue)
	cardi_mutated_list = []
	ref_sketch = [""] * 2000
	max_prime = 9999999999971.  # similar to CMash, compress hash space
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
			continue  # already found
		else:
			_mins.insert(i, h)
			_mins.pop()
			ref_sketch.insert(i, kmer)
			ref_sketch.pop()
	# delete unused cells in sketch (in case of small genome)
	while _mins[-1] == max_prime:
		_mins.pop()
		ref_sketch.pop()
	# transfer sketch to a set for faster operation, we use kmer match (not hash match for a hit)
	sketch_set = set(ref_sketch)
	
	# loop through mutated list -> stream all kmers for CMH calculation
	for mutated_string in mutated_list:
		temp_cardi = estimate_genome_size(mutated_string, kvalue)
		A_matches = dict()  # to record overlaped kmers
		for kmer in kmers(mutated_string, kvalue):
			if rev_comp:
				kmer = min(kmer, khmer.reverse_complement(kmer))
			if kmer in sketch_set:
				A_matches[kmer] = 1
		### C-MH
		C_est = np.sum(list(A_matches.values())) / len(sketch_set)
		J_est = containment_to_jaccard(C_est, cardi_A=cardi_original, cardi_B=temp_cardi)
		### store results
		out_CI.append(C_est)
		out_JI.append(J_est)
	
	# return
	return out_CI, out_JI

def s3_2_ci_ji_brute_force_validation(original_seq, mutated_list, kvalue):
	"""
	brute force validation for s3_1 calculation
	"""
	
	# get ref kmer collection
	set_original = set()
	for kmer in kmers(original_seq, kvalue):
		if rev_comp:
			kmer = min(kmer, khmer.reverse_complement(kmer))
		set_original.add(kmer)
	
	# list to store results
	out_JI = []
	out_CI = []
	
	# loop through mutated list
	for mutated_string in mutated_list:
		set_mutated_string = set()
		for kmer in kmers(mutated_string, kvalue):
			if rev_comp:
				kmer = min(kmer, khmer.reverse_complement(kmer))
			set_mutated_string.add(kmer)
		### brute force calculation of JI
		_inter = len(set_original.intersection(set_mutated_string))
		_union = len(set_original.union(set_mutated_string))
		bf_JI = _inter / _union
		out_JI.append(bf_JI)
		### brute force calculation for CI
		_size_original = len(set_original)
		bf_CI = _inter / _size_original
		out_CI.append(bf_CI)
	
	# return
	return out_CI, out_JI

def s3_3_find_mismatches_between_ci_ji(est_ci_list, bf_ci_list, est_ji_list, bf_ji_list, cutoff=0.1):
	"""
	Check the consistence between CMH and BF calculation
	"""
	if len(est_ci_list) != len(bf_ci_list) or len(est_ji_list) != len(bf_ji_list):
		raise Exception("Input lists don't have same length, please double check!")
	
	count_ci = 0
	count_ji = 0
	count_total = len(est_ci_list)
	
	for i in range(count_total):
		temp_dif1 = abs(est_ci_list[i] / bf_ci_list[i] - 1)
		if temp_dif1 > cutoff:
			count_ci += 1
		
		temp_dif2 = abs(est_ji_list[i] / bf_ji_list[i] - 1)
		if temp_dif2 > cutoff:
			count_ji += 1
	
	# summary
	ratio_CI = count_ci / count_total
	ratio_JI = count_ji / count_total
	print("Mismatch ratio for CI is " + str(ratio_CI))
	print("Mismatch ratio for JI is " + str(ratio_JI))
	return ratio_CI, ratio_JI

def s3_4_plot_ci_ji_matrix_under_given_p(df_est_ci, df_est_ji, mute_p):
	fig, axs = plt.subplots(2, 2, figsize=(20, 20))
	### ab: swarmplot
	sb.swarmplot(data=df_est_ci, ax=axs[0, 0], color="black", size=1, linewidth=0.3)
	axs[0, 0].set(ylabel="Est_CI")
	axs[0, 0].set_ylim([0, 1])
	sb.swarmplot(data=df_est_ji, ax=axs[0, 1], color="black", size=1, linewidth=0.3)
	axs[0, 1].set(ylabel="Est_JI")
	axs[0, 1].set_ylim([0, 1])
	### cd: boxplot
	sb.boxplot(data=df_est_ci, ax=axs[1, 0], color="black")
	axs[1, 0].set(ylabel="Est_CI")
	axs[1, 0].set_ylim([0, 1])
	sb.boxplot(data=df_est_ji, ax=axs[1, 1], color="black")
	axs[1, 1].set(ylabel="Est_JI")
	axs[1, 1].set_ylim([0, 1])
	### savefig
	fig.savefig("CI_JI_plot_mutation_rate_" + str(mute_p) + ".png", dpi=200)
	plt.close(fig)




# test code
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

def t1_validate_bf_ci_ji_codes():
	seq1 = "AAACCC"   # AAAC, AACC, ACCC
	seq_list = ["AAACCC", "AACCAC"]   # AACC, ACCA, CCAC
	ci_out, ji_out = s3_2_ci_ji_brute_force_validation(original_seq=seq1, mutated_list=seq_list, kvalue=4)
	assert ci_out == [1, 1/3]
	assert ji_out == [1, 0.2]

def t2_check_est_ci_ji_consistency_to_bf():
	raw_seq = s1_1_generate_string_from_genome(genome_file, rev_comp=True)
	raw_seq = raw_seq.upper()
	raw_seq = raw_seq[:10000]  # faster
	mutated_list, mismatch_ratio = s2_3_simulation_for_given_p(seq=raw_seq, p=0.08, n=100)  # 0.1 will have a very small JI
	test_est_ci, test_est_ji = s3_1_ci_ji_between_mutated_list(raw_seq, mutated_list, 25)
	test_bf_ci, test_bf_ji = s3_2_ci_ji_brute_force_validation(raw_seq, mutated_list, 25)
	mis_ratio_ci, mis_ratio_ji = s3_3_find_mismatches_between_ci_ji(test_est_ci, test_bf_ci, test_est_ji, test_bf_ji, 0.1)





### run analysis
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="ANI vs JI analysis",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', help="Path to input genomes")
	parser.add_argument('-k', '--k_range', type=str, help="k-mer range", default="15-60-5")
	parser.add_argument('-n', '--sample_size', type=int, help="Number of mutated strings.", default=500)
	parser.add_argument('-p', '--p_range', type=str, help="mutation rate range", default="0.001,0.005,0.01,0.05,0.1,0.15,0.2,0.5")
	parser.add_argument('-r', '--rev_comp', type=str, help="Whether to keep the reverse complementary", default="True")


	########################################## read parameters
	args = parser.parse_args()

	if not os.path.exists(args.genome):
		raise Exception("Input file %s does not exist." % args.genome)
	genome_file=os.path.abspath(args.genome)
	print("Input genome file is "+genome_file)

	sample_size=args.sample_size
	print("Permutation sample size is "+str(sample_size))

	rev_comp=args.rev_comp
	rev_comp = rev_comp == 'True'
	if rev_comp:
		print("Use canonical kmers in analysis!!!")

	k_range=args.k_range
	k_sizes=parsenumlist(k_range)
	print("The range of k values are: ")
	print(k_sizes)

	p_range=args.p_range
	p_values=[float(x) for x in p_range.split(",") ]
	print("The range of mutation rate p values are: ")
	print(p_values)



	########################################## Step1: raw seq data processing
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
	
	
	########################################## Step2: loop through p values and generate plot
	mean_ci = pd.DataFrame()
	mean_ji = pd.DataFrame()
	
	for mute_p in p_values:
		print("Current mutation rate is "+str(mute_p))
		# generate mutated list
		p_seq_list, p_mis_ratio_list = s2_3_simulation_for_given_p(seq=raw_seq, p=mute_p, n=sample_size)
		# confirm the edit distance for samples
		plot_list_density(p_mis_ratio_list, out_name="Mutation_ratio_" + str(mute_p) + "_size_" + str(sample_size) + ".png")
		
		# get est_CI, est_JI by C-MH method
		df_est_ci = pd.DataFrame()
		df_est_ji = pd.DataFrame()
		for kvalue in k_sizes:
			print("k"+ str(kvalue)+ " under mutation ratio of "+str(mute_p))
			df_est_ci[str(kvalue)], df_est_ji[str(kvalue)] = s3_1_ci_ji_between_mutated_list(raw_seq, p_seq_list, kvalue)
			
		# generate plot from the dict
		s3_4_plot_ci_ji_matrix_under_given_p(df_est_ci, df_est_ji, mute_p)

		# store mean value
		temp_mean_ci = []
		temp_mean_ji = []
		for kvalue in k_sizes:
			temp_mean_ci.append(round(df_est_ci[str(kvalue)].mean(), 4))
			temp_mean_ji.append(round(df_est_ji[str(kvalue)].mean(), 4))
			
		mean_ci[str(mute_p)] = temp_mean_ci
		mean_ji[str(mute_p)] = temp_mean_ji
		
	### store the mean df
	mean_ci.set_axis(["k"+str(x) for x in k_sizes], axis=0, inplace=True) #change rowname to kvalue
	mean_ji.set_axis(["k"+str(x) for x in k_sizes], axis=0, inplace=True)
	mean_ci.to_csv("Mean_est_CI_by_mute_p_kvalue.csv", index=True, encoding='utf-8')
	mean_ji.to_csv("Mean_est_JI_by_mute_p_kvalue.csv", index=True, encoding='utf-8')

	





