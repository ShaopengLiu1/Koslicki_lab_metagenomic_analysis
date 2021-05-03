import os
import pandas as pd
import numpy as np
import seaborn as sb
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import math
import statistics


####################### define functions
# read data
### clean matrix for symmatric input: pairwise JI matrix is symmatric so we only need half of the matrix
def upper_triangle(input_df, lower_tril=True):
	if lower_tril:
		print("Ploting by the lower triangle of the input matrix!!!")
		out_df = input_df.where(np.tril(np.ones(input_df.shape), -1).astype(bool))
	else:
		out_df = input_df.where(np.triu(np.ones(input_df.shape), 1).astype(bool))
	return out_df


### use a dict to store all csv for gt_CI / trunc_CI
def read_all_k_df(sect_name, k_range, symmatric=True):
	"""
	read dfs for all k values of each section and store in a dict
	:param sect_name: est, gt, trunc, bias
	:param k_range: a list of k values, e.g. [15, 20, 25]
	:param symmatric: if the matrix is symmatric, only lower half will be kept
	:param maxk_for_bf: maxk is needed for bias factor usage
	:return: a dict storing all dfs
	"""
	
	if sect_name == "est":
		prefix = "est_JI_matrix_k"
	elif sect_name == "gt":
		prefix = "GT_JI_matrix_k"
	elif sect_name == "cj":
		prefix = "containment_JI_matrix_k"
	else:
		raise Exception("Sect_name must be one of ['est', 'gt', 'cj']")
	
	out_dict = dict()
	
	for k in k_range:
		f_name = prefix + str(k) + "_out.csv"
		# read df for a single k value
		temp_df = pd.read_csv(f_name, header=0, index_col=0)
		if symmatric:
			temp_df = upper_triangle(temp_df)
		out_dict[str(k)] = temp_df
	
	return out_dict


### get dif matrix
def get_dif_matrix(dict1, d_gt, k_range):
	# merge out values into a big DF
	abs_out = pd.DataFrame()
	rela_out = pd.DataFrame()
	for k in k_range:
		temp_colname = 'k=' + str(k)
		# abs error
		dif_df = dict1[str(k)] - d_gt[str(k)]
		all_abs_values = dif_df.stack().to_numpy()
		abs_out[temp_colname] = all_abs_values
		# rela error
		rela_df = dif_df / (d_gt[str(k)] + 0.00001)
		all_rela_values = rela_df.stack().to_numpy()
		rela_out[temp_colname] = all_rela_values
		
	return rela_out


### generate single boxplot for rela error
def single_box_of_df(input_df, out_name):
	fig, axs = plt.subplots(1, 1, figsize=(10, 10))
	plt.rcParams.update({'font.size': 10})
	sb.boxplot(data=input_df, ax=axs, color="black")
	axs.set(ylabel=out_name + " - Ground truth")
	axs.set_ylim([-0.3, 0.3])
	fig.savefig("rela_error_" + out_name + ".png", dpi=200)
	plt.close(fig)


### generate boxplot for the unbalanced case
def multibox_ub_rela_error(input_df, out_name):
	fig, axs = plt.subplots(1, 1, figsize=(10, 10))
	plt.rcParams.update({'font.size': 10})
	sb.boxplot(data=input_df, ax=axs, color="black")
	axs.set(ylabel="Relative error to GT")
	axs.set_ylim([-0.5, 0.5])
	fig.savefig("UB_rela_error_" + out_name + ".png", dpi=200)
	plt.close(fig)



# local test
def local_test():
	work_dir = '/Users/shaopeng/Desktop/local_CMash_dev'
	f1 = work_dir + '/566_data/out1_simi/'
	f2 = work_dir + '/566_data/out2_ub/'



# for pipe purpose
work_dir = os.getcwd()
f1 = work_dir + '/out1_similar_size_output/'
f2 = work_dir + '/out2_unbalance_setting'




# get fig1
os.chdir(f1)
os.listdir()
k_range=list(range(15, 61, 5))

dict_est = read_all_k_df('est', k_range=k_range, symmatric=True)
dict_gt = read_all_k_df('gt', k_range=k_range, symmatric=True)
dict_cj = read_all_k_df('cj', k_range=k_range, symmatric=True)

dif_df_est = get_dif_matrix(dict_est, dict_gt, k_range=k_range)
dif_df_cj = get_dif_matrix(dict_cj, dict_gt, k_range=k_range)
# plot dif_est
single_box_of_df(dif_df_est, "Est_JI")
single_box_of_df(dif_df_cj, "Con_JI")


# plot the unbalance situation
os.chdir(f2)
### read data:
df_gt = pd.read_csv("GT_JI_matrix_k30_out.csv", header=0, index_col=0)
est_m5_rela = pd.DataFrame()
cj_m5_rela = pd.DataFrame()

df_m5_est = pd.read_csv("est_JI_matrix_k30_m500.csv", header=0, index_col=0)
df_m5_cj = pd.read_csv("containment_JI_matrix_k30_m2000.csv", header=0, index_col=0)
df_rela_m5_est = (df_m5_est - df_gt) / df_gt
df_rela_m5_cj = (df_m5_cj - df_gt) / df_gt
df_rela_m5_est.columns.values
df_rela_m5_est.columns = ["1:3", "1:8", "1:22", "1:1.5"]
df_rela_m5_cj.columns = ["1:3", "1:8", "1:22", "1:1.5"]
# make box plot
multibox_ub_rela_error(df_rela_m5_est.iloc[:,0:3], "Est_ji")
multibox_ub_rela_error(df_rela_m5_cj.iloc[:,0:3], "Con_ji")
# box plot of GT value
fig, axs = plt.subplots(1, 1, figsize=(10, 10))
plt.rcParams.update({'font.size': 10})
sb.boxplot(data=df_gt.iloc[:,0:3], ax=axs, color="black")
axs.set(ylabel="GT_JI distribution")
axs.set_ylim([0, 1])
fig.savefig("GT_JI_distribution.png", dpi=200)
plt.close(fig)
# compare rela error when increasing random sample
df_increase_sample_size = pd.DataFrame()
for sample_size in ['m500', 'm1000', 'm2000', 'm5000', 'm10000']:
	temp_df = pd.read_csv("est_JI_matrix_k30_" + sample_size + ".csv", header=0, index_col=0)
	temp_rela = (temp_df - df_gt) / df_gt
	df_increase_sample_size[sample_size] = temp_rela['merged_f14.fa.gz']
# generate plot of increasing sample size
fig, axs = plt.subplots(1, 1, figsize=(10, 10))
plt.rcParams.update({'font.size': 10})
sb.boxplot(data=df_increase_sample_size, ax=axs, color="black")
axs.set(ylabel="Relative error")
axs.set_ylim([-0.5, 0.5])
fig.savefig("Increase_sample_size_rela_error.png", dpi=200)
plt.close(fig)




