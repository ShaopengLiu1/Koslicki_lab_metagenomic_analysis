import argparse
import os
import glob
import pandas as pd
import numpy as np
import sys
import seaborn as sb
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import re

### def functions
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
		raise argparse.ArgumentTypeError(
			"'" + k_sizes_str + "' is not a range of number. Expected forms like '1-5' or '2' or '10-15-2'.")
	start = int(m.group(1))
	end = int(m.group(2))
	if m.group(3):
		increment = int(m.group(3))
	else:
		increment = 1
	return list(range(start, end + 1, increment))

# clean matrix for symmatric input
def upper_triangle(input_df):
	if lower_tril:
		print("Ploting by the lower triangle of the input matrix!!!")
		out_df = input_df.where(np.tril(np.ones(input_df.shape), -1).astype(np.bool))
	else:
		out_df = input_df.where(np.triu(np.ones(input_df.shape), 1).astype(np.bool))
	return out_df

# use dict to store all dfs for each section
def read_all_k_df(sect_name, k_range, symmatric=False, maxk_for_bf=None):
	"""
	read dfs for all k values of each section
	:param sect_name: est, gt, trunc, bias
	:param k_range: [15, 20, 25]
	:param symmatric: if the matrix is symmatric, only upper half will be kept
	:param maxk_for_bf: maxk is needed for bias factor usage
	:return: a dict storing all dfs
	"""
	
	if sect_name == "est":
		prefix = "est_JI_k"
	elif sect_name == "gt":
		prefix = "GroundTruth_JI_k"
	elif sect_name == "trunc":
		prefix = "trunc_JI_k"
	elif sect_name == "bias":
		prefix = "bias_factor_k"
		if maxk_for_bf == None:
			raise Exception("Please specify maxk when reading BF dfs")
	else:
		raise Exception("Sect_name must be one of ['est', 'gt', 'trunc', 'bias']")
	
	out_dict=dict()
	
	for k in k_range:
		if sect_name == "bias":
			f_name = prefix + str(k) + "_to_k" + str(maxk_for_bf) + ".csv"
		else:
			f_name = prefix + str(k) + ".csv"
		# read df for a single k value
		temp_df = pd.read_csv(f_name, header=0, index_col=0)
		if symmatric:
			temp_df = upper_triangle(temp_df)
		out_dict[str(k)] = temp_df
		
	return out_dict

# plot single section of dfs
def plot_single_section_all_k(input_dict, k_range, plot_name, ylab="JI value"):

	# merge the values into a big DF
	p_out = pd.DataFrame()
	for k in k_range:
		temp_colname = 'k=' + str(k)
		all_values_1k = input_dict[str(k)].stack().to_numpy()
		p_out[temp_colname] = all_values_1k
	
	# generate plot
	fig, axs = plt.subplots()
	sb.boxplot(data=p_out, ax=axs, color="cyan")
	fig.suptitle(plot_name)
	axs.set(ylabel=ylab)
	fig.savefig(plot_name + ".png", dpi=200)
	plt.close(fig)
	
	# output dfs if needed
	return p_out

# plot abs/rela error of 2 dfs
def plot_abs_rela_error_all_k(dict1, dict2, k_range, plot_name):
	
	# merge out values into a big DF
	abs_out = pd.DataFrame()
	rela_out = pd.DataFrame()
	for k in k_range:
		temp_colname = 'k=' + str(k)
		# abs error
		dif_df = dict1[str(k)] - dict2[str(k)]
		all_abs_values = dif_df.stack().to_numpy()
		abs_out[temp_colname] = all_abs_values
		# rela error
		rela_df = dif_df / (dict2[str(k)]+0.00001)
		all_rela_values = rela_df.stack().to_numpy()
		rela_out[temp_colname] = all_rela_values
		
	# generate plot
	fig, axs = plt.subplots(1, 2, figsize=(16, 12))
	sb.boxplot(data=abs_out, ax=axs[0], color="black")
	axs[0].set(ylabel="Truncation - Ground truth")
	axs[0].title.set_text("Value change: Truncation-based JI vs Ground truth")
	sb.boxplot(data=rela_out, ax=axs[1], color="black")
	axs[1].set(ylabel="Difference as a ratio of Ground truth")
	axs[1].title.set_text("Ratio change: difference over the value of Ground truth")
	fig.suptitle(plot_name)
	fig.savefig("Trunc_JI_vs_GT_" + plot_name + ".png", dpi=200)
	plt.close(fig)
	
	# return
	return [abs_out, rela_out]
		
	
input_range="15-60-5"
k_sizes = parsenumlist(input_range)
lower_tril = True
maxk=60

plot_est = read_all_k_df('est', k_sizes, True)
all_est_jis = plot_single_section_all_k(plot_est, k_sizes, "Boxplot_of_est_JI")




