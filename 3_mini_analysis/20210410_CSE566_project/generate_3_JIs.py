################################## import modules
from multiprocessing import Pool
import multiprocessing
import argparse
from argparse import ArgumentTypeError
import os
import khmer
from itertools import *
import numpy as np
import pandas as pd
import re
import subprocess
import screed
import bisect





################################## define functions
# Read all file paths from an input file and check if they exist
def check_files(input_file):
	print("Reading paths from %s to a list." %input_file)
	out_list = list()
	with open(input_file, 'r') as temp:
		for line in temp.readlines():
			line = line.strip()
			if not os.path.exists(line):
				raise Exception("Input file %s does not exist." % line)
			out_list.append(os.path.abspath(line))
		out_list = sorted(out_list, key=os.path.basename)
		return out_list

# fast kmer cardinality estimate by hyperloglog
def estimate_genome_size(genome, input_ksize):
	hll = khmer.HLLCounter(0.01, input_ksize)
	hll.consume_seqfile(genome)
	return hll.estimate_cardinality()

# transfer C(A, B) to J(A, B) by definiation for 2 CE objects
def containment_to_jaccard(C_est, CE1, CE2):
	A = CE1.cardinality
	B = CE2.cardinality
	return C_est * A * 1.0 / (A + B - C_est * A)

# canonical form of kmer
def canonical_kmer(kmer):
	return min(kmer, khmer.reverse_complement(kmer))

###### several wrapper function for parallel running of est_J, gt_J, and J_from_C based on CE objects
# 1 vs list for a vector of est_JIs
def wrapper_est_jaccard_vector(arg):
	"""
	wrapper function for CCE object: est_jaccard function, 1 vs a list
	:param arg:
	:return:
	"""
	return arg[0].est_jaccard(arg[1])

# list vs list for a matrix of est_JIs
def wrapper_est_jaccard_matrix(list1, list2, out_file_name, threads=8):
	"""
	wrapper function to get a pairwise est_JI matrix for CCE object list1 vs list2
	:param threads:
	:param list1:
	:param list2:
	:param out_file_name:
	:param thread:
	:return: return nothing but an output file with the matrix
	"""
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].input_file_name) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].input_file_name)
		row_name.append(name)
		Y = list1[i].ji_vector(list2, threads=threads)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return

# 1 vs list for a vector of J_Cs
def wrapper_streaming_method_ji_from_ci_vector(arg):
	"""
	wrapper function for streaming method CI -> JI in CCE obj to handle 1 vs a list of objects
	:param arg:
	:return:
	"""
	return arg[0].stream_cmash_for_ji(arg[1])

# list vs list for a matrix of J_Cs
def wrapper_streaming_method_ji_from_ci_matrix(list1, list2, out_file_name, threads=8):
	"""
	wrapper functino to generate a JI matrix for 2 lists of CCE objects
	:param list1:
	:param list2:
	:param out_file_name:
	:param threads:
	:return:
	"""
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].input_file_name) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].input_file_name)
		row_name.append(name)
		Y = list1[i].cmash_stream_ji_vector(list2, threads=threads)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return df

# 1 vs list for a vector of gt_JIs
def wrapper_ground_truth_ji_vector(arg):
	"""
	wrapper function for ground truth JI calculation for CCE objects 1 vs a list
	:param arg:
	:return:
	"""
	return arg[0].ground_truth_jaccard(arg[1])

# list vs list for a matrix of gt_JIs
def wrapper_ground_truth_ji_matrix(list1, list2, out_file_name, threads=8):
	"""
	wrapper functino to generate a GT_JI matrix for 2 lists of CCE objects
	:param list1:
	:param list2:
	:param out_file_name:
	:param threads:
	:return:
	"""
	lst = []
	row_name = []
	col_name = [os.path.basename(list2[i].input_file_name) for i in range(len(list2))]
	for i in range(len(list1)):
		name = os.path.basename(list1[i].input_file_name)
		row_name.append(name)
		Y = list1[i].ground_truth_ji_vector(list2, threads=threads)
		lst.append(Y)
	df = pd.DataFrame(lst)
	df.columns = col_name
	df.index = row_name
	df.to_csv(out_file_name, index=True, encoding='utf-8')
	return





################################## define the CountEstimator object (CE object)
# CE object is the major part of this script, it will store all important infor of a given file
# using CE objects would need extra running time but it's important for algorithm illustration purpose

### dependent function
# check for potential breakup symbol
notACTG = re.compile('[^ACTG]')

# yield all k-mers of len ksize from seq.
def kmers(seq, ksize):
	for i in range(len(seq) - ksize + 1):
		yield seq[i:i + ksize]

### object definition
class CMash_CountEstimator:
	
	# initialization
	def __init__(self, n=1000, max_prime=9999999999971., ksize=31, input_file_name=None, rev_comp=True,
	             full_kmer=False, threads=8):
		"""
		:param n: sketch size
		:param max_prime: use a large prime to rehash the values so they are more evenly distributed
		:param ksize: kmer length
		:param input_file_name:
		:param rev_comp: use canonical kmers
		:param full_kmer: store all kmers instead of sketches only, required for gt_JI calculation
		:param threads: parallel usage
		"""
		
		# read parameters
		self.ksize = ksize  # max k value for the truncation
		self.rev_comp = rev_comp  # binary indicator: use canonical k-mer, default Yes
		self.full_kmer = full_kmer  # binary indicator: store full k-mer sets besides the n-sketch for CMash validation
		self.threads = threads
		self.p = max_prime
		self.input_file_name = input_file_name
		
		# initialize kmer/hash vector
		self._mins = [self.p] * n  # hash sketch
		self._kmers = [''] * n  # kmer sketch
		self._all_kmer = set()  # use set for the project purpose, use dict for production code (using weights)
		
		# build bottom sketch
		if self.input_file_name:
			self.cardinality = estimate_genome_size(input_file_name, self.ksize)
			for record in screed.open(self.input_file_name):
				self.add_sequence(record.sequence, rev_comp=self.rev_comp)
	
	# add sequence
	def add_sequence(self, seq, rev_comp):
		"""
		 Sanitize and add a sequence to the sketch.
		"""
		seq = seq.upper()  # otherwise their hashes will be different
		seq_split_onlyACTG = notACTG.split(seq)
		if len(seq_split_onlyACTG) == 1:
			for kmer in kmers(seq, self.ksize):
				self.add(kmer, rev_comp=rev_comp)
		else:
			for sub_seq in seq_split_onlyACTG:
				if sub_seq:
					self.add_sequence(sub_seq, rev_comp=rev_comp)
	
	# add (continued)
	def add(self, kmer, rev_comp):
		"""
		Add kmer into sketch, keeping sketch sorted
		Also add to full kmer set if needed
		"""
		_mins = self._mins
		_kmers = self._kmers
		
		# canonical form if required
		if rev_comp:
			kmer = canonical_kmer(kmer)
		
		# get hash value
		h = khmer.hash_no_rc_murmur3(kmer)
		h = h % self.p
		
		# insert into full kmer set if needed
		if self.full_kmer:
			self._all_kmer.add(kmer)
			
		# insert into sketch
		if h >= _mins[-1]:
			return
		else:
			i = bisect.bisect_left(_mins, h)  # find index to insert h
			_mins.insert(i, h)
			_mins.pop()
			_kmers.insert(i, kmer)
			_kmers.pop()
			return
	
	# estimate JI between 2 CE objects, slightly dif from CMash (same idea)
	def est_jaccard(self, other):
		if self.ksize != other.ksize:
			raise Exception("different k-mer sizes - cannot compare")
		if self.p != other.p:
			raise Exception("different primes - cannot compare")
		i = 0
		j = 0
		overlap = 0
		processed = 0
		mins1 = self._mins
		mins2 = other._mins
		random_sample_size = min(len(mins1), len(mins2))
		while processed < random_sample_size:  # loop stop when min(|A|, |B|) sample size reached
			processed += 1
			if mins1[i] < mins2[j]:
				i += 1
			elif mins1[i] > mins2[j]:
				j += 1
			elif mins1[i] == mins2[j]:
				# a match from A overlap B
				overlap += 1
				i += 1
				j += 1
		est_ji = overlap * 1.0 / processed
		print("est_JI: " + str(est_ji))
		return est_ji
	
	# generate a JI vector for self vs a list of CCE objects, "unwrap_jaccard_vector" from MH
	def ji_vector(self, other_list, threads=8):
		pool = Pool(processes=threads)
		Y = np.array(pool.map(wrapper_est_jaccard_vector, zip([self] * len(other_list), other_list)))
		pool.terminate()
		return Y
	
	# streaming method to est JI from CI
	def stream_cmash_for_ji(self, other):
		"""
		Est JI from CI by streaming method
		We can only get truncated_JI in this way. It's not correct to directly truncate the hash lists from 2 CCE objects (wrong sampling method)
		:param other:
		:return:
		"""
		if self.ksize != other.ksize:
			raise Exception("different k-mer sizes - cannot compare")
		if self.p != other.p:
			raise Exception("different primes - cannot compare")
		if not other.full_kmer:
			raise Exception("Need to enable full_kmer feature for the genome being streamed")
		
		A_kmers = self._kmers
		A_matches = 0
		for kmer in A_kmers:
			if kmer in other._all_kmer:
				A_matches += 1
		
		# return results
		C_est = A_matches*1.0 / len(A_kmers)
		J_est = containment_to_jaccard(C_est, self, other)
		print("Containment_JI: " + str(J_est))
		return J_est
	
	# generate a trunc_JI vector (from CI by streaming method) for self vs a list of CCE objects
	def cmash_stream_ji_vector(self, other_list, threads=8):
		pool = Pool(processes=threads)
		Y = np.array(pool.map(wrapper_streaming_method_ji_from_ci_vector, zip([self] * len(other_list), other_list)))
		pool.terminate()
		return Y
	
	# ground truth JI: discard kmc which takes toooooo many storage spaces
	def ground_truth_jaccard(self, other):
		if self.ksize != other.ksize:
			raise Exception("different k-mer sizes - cannot compare")
		if self.p != other.p:
			raise Exception("different primes - cannot compare")
		if not self.full_kmer or not other.full_kmer:
			raise Exception("full kmer not enabled for the CE object")
		
		_inter = self._all_kmer.intersection(other._all_kmer)
		_union = self._all_kmer.union(other._all_kmer)
		gt_ji = len(_inter) * 1.0 / len(_union)
		print("GT_JI: "+str(gt_ji))
		return gt_ji
	
	# get a GT_JI vector
	def ground_truth_ji_vector(self, other_list, threads=8):
		pool = Pool(processes=threads)
		Y = np.array(pool.map(wrapper_ground_truth_ji_vector, zip([self] * len(other_list), other_list)))
		pool.terminate()
		return Y






################################## Script running functions
# workflow related wrappers

# a wrapper function to create CE object
def make_ce_object(genome, max_h, ksize, rev_comp, full_kmer, threads):
	CE = CMash_CountEstimator(n=max_h, ksize=ksize, input_file_name=genome, rev_comp=rev_comp, threads=threads, full_kmer=full_kmer)
	return CE

# wrapper to create CE objects in parallel by pool
def wrapper_make_ce_star(arg):
	return make_ce_object(*arg)

# create a list of CE objects
def get_ce_lists(input_file, max_h, input_k, rev_comp, full_kmer, threads):
	para = Pool(processes=threads)
	genome_sketches = para.map(wrapper_make_ce_star, zip(input_file, repeat(max_h), repeat(input_k), repeat(rev_comp), repeat(full_kmer), repeat(threads)))
	para.close()
	return genome_sketches







################################## test the script on known data
# test GT_JI
def f1_test_gt_ji():
	subprocess.run(f"echo '>test_seq\nTTAATTAA' > test1.fq", shell=True)
	subprocess.run(f"echo '>test_seq\nAAATTTAA' > test2.fq", shell=True)
	obj1 = make_ce_object(genome="test1.fq", max_h=10, ksize=4, rev_comp=True, full_kmer=True, threads=1)
	obj2 = make_ce_object(genome="test2.fq", max_h=10, ksize=4, rev_comp=True, full_kmer=True, threads=1)
	os.remove("test1.fq")
	os.remove("test2.fq")
	#print("This is for canonical kmers: |A|=3, |B|=4, overlap=2")
	temp_gt_ji = obj1.ground_truth_jaccard(obj2)
	assert temp_gt_ji == 2 / 5.0

# test est_JI
def f2_test_est_ji():
	E1 = CMash_CountEstimator(n=0, ksize=5)
	E2 = CMash_CountEstimator(n=0, ksize=5)
	E1._mins = [1, 2, 3, 4, 5]
	E2._mins = [1, 2, 3, 4]
	assert E1.est_jaccard(E2) == 4 / 4.0    #sample until min(|A|, |B|), which is 4 here
	assert E2.est_jaccard(E1) == 4 / 4.0
	E2._mins = [1, 3, 4, 5, 6]
	assert E1.est_jaccard(E2) == 4 / 5.0









################################## local test
# local test at personal laptop
def local_test():
	input_file = 'file_path.txt'
	ksize = 30
	threads = 10
	
	query_list = check_files(input_file)
	ref_list = check_files(input_file)
	file1 = query_list[0]
	file2 = query_list[1]
	
	# single file test
	E1 = make_ce_object(genome=file1, max_h=2000, ksize=ksize, rev_comp=True, full_kmer=True, threads=threads)
	E2 = make_ce_object(genome=file2, max_h=2000, ksize=ksize, rev_comp=True, full_kmer=True, threads=threads)
	print(E1.cardinality)
	print(E2.cardinality)
	E1.est_jaccard(E2)
	E1.ground_truth_jaccard(E2)
	E1.stream_cmash_for_ji(E2)
	
	# list vs list for JI matrix
	ce_list1 = get_ce_lists(query_list, max_h=2000, input_k=ksize, rev_comp=True, full_kmer=True, threads=5)
	ce_list2 = get_ce_lists(ref_list, max_h=2000, input_k=ksize, rev_comp=True, full_kmer=True, threads=5)
	###
	out_name = "test_est_JI.csv"
	wrapper_est_jaccard_matrix(ce_list1, ce_list2, out_name, threads=5)
	### non-parallel running mode is faster: possible explaination: https://stackoverflow.com/questions/38217449/python-multiprocessing-slower-than-single-thread
	for i in ce_list1:
		for j in ce_list2:
			i.est_jaccard(j)
			
	
	






################################## Pipeline script!
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="CSE566 project analysis",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-q', '--query', help="Path to the query file.")
	parser.add_argument('-r', '--ref', help="Path to the ref file.")
	parser.add_argument('-c', '--rev_comp', type=str, help="Whether to keep the reverse complementary", default="True")
	parser.add_argument('-n', '--num_hashes', type=int, help="Number of hashes to use.", default=2000)
	parser.add_argument('-k', '--k_size', type=int, help="k-mer size", default=30)
	parser.add_argument('-t', '--threads', type=int, help="Number of threads to use",
	                    default=min(12, int(multiprocessing.cpu_count() / 2)))
	parser.add_argument('-p', '--prime', help='Prime (for modding hashes)', default=9999999999971)
	# for partial results
	parser.add_argument('-x', '--skip_est', type=str, help="Skip Est JI calculation", default="False")
	parser.add_argument('-y', '--skip_containment', type=str, help="Skip Containment JI calculation", default="False")
	parser.add_argument('-z', '--skip_gt', type=str, help="Skip GT JI calculation", default="False")
	
	# read parameters
	args = parser.parse_args()
	num_threads = args.threads
	prime = args.prime  # taking hashes mod this prime
	ksize = args.k_size  # max k used in the analysis
	max_h = args.num_hashes  # number of hashes to use
	query_file = args.query
	ref_file = args.ref
	rev_comp = args.rev_comp
	skip_gt = args.skip_gt
	skip_est = args.skip_est
	skip_containment = args.skip_containment
	
	
	
	# setup
	rev_comp = rev_comp == 'True'
	if rev_comp:
		print("Using true for rev_comp!!!!!!")
	
	skip_gt = skip_gt == 'True'
	if skip_gt:
		print("Skipping all GT JI calculation steps")
	
	skip_est = skip_est == 'True'
	if skip_est:
		print("Skipping all Est JI calculation steps")
		
	skip_containment = skip_containment == 'True'
	if skip_containment:
		print("Skipping all Containment JI calculation steps")
	
	query_file_names = os.path.abspath(query_file)
	if not os.path.exists(query_file_names):
		raise Exception("Input file %s does not exist." % query_file_names)
	
	ref_file_names = os.path.abspath(ref_file)
	if not os.path.exists(ref_file_names):
		raise Exception("Input file %s does not exist." % query_file_names)
	
	
	
	
	# Pipe start:
	### read files
	query_list = check_files(query_file)
	ref_list = check_files(ref_file)
	
	### build CE objects
	if skip_containment and skip_gt:
		full_kmer = False
	else:
		full_kmer = True  #slightly faster if only est_JI needed
	ce_list1 = get_ce_lists(query_list, max_h=max_h, input_k=ksize, rev_comp=rev_comp, full_kmer=full_kmer, threads=num_threads)
	ce_list2 = get_ce_lists(ref_list, max_h=max_h, input_k=ksize, rev_comp=rev_comp, full_kmer=full_kmer, threads=num_threads)
	
	### generate results csv
	if not skip_est:
		print("Calculating est_JI......")
		out_name = "est_JI_matrix_k" + str(ksize) + ".csv"
		wrapper_est_jaccard_matrix(ce_list1, ce_list2, out_name, threads=num_threads)
	if not skip_containment:
		print("Calculating Containment_JI......")
		out_name = "containment_JI_matrix_k" + str(ksize) + ".csv"
		wrapper_streaming_method_ji_from_ci_matrix(ce_list1, ce_list2, out_name, threads=num_threads)
	if not skip_gt:
		print("Calculating GT_JI......")
		out_name = "GT_JI_matrix_k" + str(ksize) + ".csv"
		wrapper_ground_truth_ji_matrix(ce_list1, ce_list2, out_name, threads=num_threads)
			
	print("Script done")
	
	

