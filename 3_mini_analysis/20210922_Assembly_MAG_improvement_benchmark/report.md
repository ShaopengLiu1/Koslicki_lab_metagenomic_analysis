## MAG improvement benchmark analysis

### Information

1. Local location:

   ```bash
   /data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark
   ```

2. Raw data:

   1. Data information: 
      1. CAMI2 challenge (simulated): extremely high strain diversity dataset (strain madness)
      2. CAMI DB: https://data.cami-challenge.org/participate (make an account)
         1. strain madness: https://www.biorxiv.org/content/10.1101/2021.07.12.451567v1.full.pdf
      3. FTP location: https://portal.nersc.gov/dna/RD/Adv-Seq/StrainMadness/
   2. Assembly results from MetaHipmer (copied, "final_assembly.fasta.gz"): `/data/temp/For_Shaopeng`

3. Data information:

   1. reads: 100 raw reads fq files
   2. final_assembly.fasta.gz: pool of multi-assembly. The depth file is abundance in 100 read files.
   3. assembly file: ****
      1. pooled gold standard assembly: 35318 (records), long and complete, should be the contig from input genomes **(confirmed: this is equivalent to the truth, raw data are not open now)**
      2. pooled megahit assembly: 1404890 (records), short contigs, the numbers make sense. This should be assembies from pooled data. **(confirmed: this is the multi-assembly output)**
      3. final assembly.fa: 40603 pooled contigs, with average length ~3k. Looks like all contigs have been aggregated to reach such numbers? Possibly lose of resolution due to aggregation. **(confirmed: co-assembly output from Metahipmer2)**
   4. Note: our lab server is not appropriate to run Metahipmer2 (needs large MEM).

4. Refs:

   1. Metahipmer2 help: https://hpcadvisorycouncil.atlassian.net/wiki/spaces/HPCWORKS/pages/1827307543/MetaHipMer+2.0+for+ISC21#Build-and-Install-mhm2
   2. Install (No Github or Conda): https://bitbucket.org/berkeleylab/mhm2/src/master/
   3. https://bitbucket.org/berkeleylab/mhm2/src/master/docs/mhm_guide.md (check here)
   
4. Questions:

   - [ ] Wha't the number in MHM2 output: ">Contig9178 48.770420", values range from 0 to 144.95. Depth? 



---

### Analysis plan

#### Hypothesis: 

1. new metagenomic sequence data can improve current assembly results
2. read-derived strain bins is better than species genome bins



#### Data:

1. Randomly assign the 100 sequencing data to "**current**" and "**new**" group by ratio 6:4. 
   1. "Current" group is used to build the current database
   2. "New" groups servers as the new data that can potentially increase assembly results
2. Gold standard: the provided pooled gold standard seesmbly file 



#### Procedure:

<u>For hypothesis 1:</u> 

1. build the current database $D_0$ as a starting point:
   1. **Option1**: run MetaHipMer2 for all data in the "current" group to get co-assembly results, then cluster them into genome bins
   2. **Option2**: 
      1. run MetaHit for all data in the "current" group to get multi-assembly results
      2. pool all the multi-assembly results together, deduplication, identify and remove sequence errors (the number of contigs is overwhemling here)
      3. merge all unitigs in the pooled results
      4. cluster all contigs into genome bins

2. for every data $i$ in the "new" group:
   1. run **MetaHit** for de novo assembly $A_i$
   2. assign all contigs from $A_i$ to current genome bin by some similarity cutoff (e.g. 95% ANI by **FastANI**)
   3. for remaining contigs, cluster them and create new genome bins
   3. **SGB dedup (TBD)** 
   4. update current database from $D_{i-1}$ to $D_i$
   
3. Select some $D_i$ to compare the quality of assembly.



<u>For hypothesis 2:</u>

1. Build **species genome bins (SGBs)**, same as $D_0$ in hypothesis 1 plan, but use full 100 data (which is the co-assembly output)
2. Map raw reads from **ALL** datasets to SGBs (BBTools), then repeat the assembly for each SGB (**STRONG**). We then get **read-derived-strain-bins**.
3. Compare SGBs vs read-derived-strain-bins



---

### Benchmark interactive code

```bash
#Local file
cd /data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark



# env preparation (use "metagenomic_py37" for now, will build separate one in full analysis)
work_dir=$PWD
data_dir=${work_dir}/data
src_dir=${work_dir}/src
out_dir=${work_dir}/output
conda activate metagenomic_py37



# run megahit v1.2.9 (new version) for all data
mkdir -p ${out_dir}/megahit_single_assembly

for single_data in ls ${data_dir}/raw_reads_100; do
        echo "Running ${single_data} now......"
        name=$(echo ${single_data%.fq.gz})
        megahit -r ${data_dir}/raw_reads_100/${single_data} -o ${out_dir}/megahit_single_assembly/megahit_out_${name}
        unset single_data name
done



# exploration:
### task1: clustering alg (based on fastANI?)

### task3: map raw reads from **ALL** datasets to SGBs (BBTools)
```



### Unit test (notes)

```bash
# fastANI
### code
$ fastANI -q genome1.fa -r genome2.fa -o output.txt
$ fastANI -q genome1.fa --rl genome_list.txt -o output.txt
### note
1. needs to separate fa file, each file represent a genome. 
2. For our analysis, we need to cluster first to get SGBs then assign every new assembled fragments to the SGB
3. extremely fast
4. we may need to cluster assembled fragments too before running


# cluster of fragments:
### by ANI
1. needs to disect to single files


### by cosine similarity of TNF (VAE paper)


```

