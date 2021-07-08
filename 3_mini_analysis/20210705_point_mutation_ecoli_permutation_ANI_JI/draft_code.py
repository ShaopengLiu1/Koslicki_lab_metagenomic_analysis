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



### function definition
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



def generate_string_from_genome(input_genome):
	"""
	Transfer input genome to a string
	"""
	out_string = input_genome
	return out_string



def generate_mutated_strings(input_string, p, n):
	"""
	Generate n mutated strings by simple mutation model
	"""
	out_str_list=[]
	return out_str_list



### run analysis
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="ANI vs JI analysis",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', help="Path to input genomes")
	parser.add_argument('-k', '--k_range', type=str, help="k-mer range", default="15-60-5")
	parser.add_argument('-n', '--sample_size', type=int, help="Number of mutated strings.", default=1000)
	parser.add_argument('-p', '--p_range', type=str, help="mutation rate range", default="0.05-0.9-0.05")
	parser.add_argument('-r', '--rev_comp', type=str, help="Whether to keep the reverse complementary", default="True")


	# read parameters
	args = parser.parse_args()

	if not os.path.exists(args.genome):
		raise Exception("Input file %s does not exist." % args.genome)
	genome_file=os.path.abspath(args.genome)

	sample_size=args.sample_size

	rev_comp = rev_comp == 'True'

	k_range=args.k_range
	k_sizes=parsenumlist(input_range)

	p_range=args.p_range
	p_values=parsenumlist(p_range)



	# run analysis







