### Analysis plan for JI vs ANI in E.coli permutation model

---

#### Aim: 

Correlates average nucleotide identity (ANI) and Jaccard index (JI) by k-mers in range of k values.



#### Parameter:

1. Input genome (default: e.coli for now), the baseline
2. k range (default: 15-60 by 5), length for k-mer
3. p range (default: 0.001, 0.005, 0.1, 0.15, 0.5), simple point mutation
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
	
	
```

