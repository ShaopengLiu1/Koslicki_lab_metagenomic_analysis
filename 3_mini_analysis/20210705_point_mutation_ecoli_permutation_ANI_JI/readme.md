### Analysis plan for JI vs ANI in E.coli permutation model

---

### Local files: 

```
/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210705_point_mutation_ecoli_permutation_ANI_JI
```

(Processing, will be summarized to Onedrive)



#### Aim: 

Correlates average nucleotide identity (ANI) and Jaccard index (JI) by k-mers in range of k values.



#### Parameter:

1. Input genome (default: e.coli for now), the baseline
2. k range (default: 15-60 by 5), length for k-mer
3. p range (default: 0.001, 0.005, 0.01,0.05,0.1, 0.15, 0.5), simple point mutation
4. simulation size m (default: 10000)



#### More details to confirm:

1. independent & simple point mutation: 
   1. generate n numbers ~ Unif(0,1) and use p as the cutoff for 0/1
   2. mutation target: equally to the rest 3 letters
2. Expectation of JI/CI: that's empirical mean or median?



#### Analysis plan:

1. data preparation:
   1. transfer input genome -> a long string S
2. loop through mutation rate p in p range, for each p value:
   1. generate m mutated strings
   2. for every mutated string S', get the following between S and S':
      1. ANI: close to p, for QC purpose only in case of small n
      2. JI for all k values
      3. CI for all k values
   3. plot the distribution of JI, CI by k value. Get the **mean/median** to compare with the theoretical expectation.



#### Output:

1. many line/box plots of ANI vs JI: x axis is k, y axis is JI.
2. many tables records the expectation of JI/CI 







---

### Analysis records

#### Step1, get a E. coli data

1. There are 260 E.coli data uploaded in 2021, just randomly pick 1:

```
Infor: GCA_002220215.1, strain=E62, complete genome & full
link: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/220/215/GCA_002220215.1_ASM222021v1
```



2. download data

```
file="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/220/215/GCA_002220215.1_ASM222021v1"
suffix=$(echo ${file##*/})
echo "downloading $suffix"
wget ${file}/${suffix}_genomic.fna.gz 2>/dev/null  
gunzip GCA_002220215.1_ASM222021v1_genomic.fna.gz
```



3. paste fq into a single string

```
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
```



#### Step2, generate m mutated sequences

1. simple mutation model by Unif(0,1)

```
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
```



2. repeat for n times, the mismatch ratio is supposed to be normally distributed with small variance. Start from here, all downstream steps will be merged into a loop through all mutation rates.

```
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
```



#### Step3, generate JI/CI for a k value by Contaiment MH (sketch size=2000)

Shadowed from CMash, under coding

```
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
	
	
```



### Overall process:

```bash
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
```

