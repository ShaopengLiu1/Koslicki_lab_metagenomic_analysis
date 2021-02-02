# Overview

Date: 2021.1.31

1. Aim: simulate error-free and regular (with error) metagenomic reads by CAMISIM

2. Input: a list of strains

3. Output: metagenomic reads stored in fastq file

4. Result location (this is one of the 10 replicates)

```
/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210130_CAMISIM_errorfree_simulation/output/tailored_output

### file infor:
1. distribution_*_sample_0.txt: abundance of the simulated data
2. new_genome_to_id.tsv: genome id and the corresponding file
3. *.fq.gz: simulated fastq file
```

5. Reproducibility: all interactive scripts can be found [here](https://github.com/ShaopengLiu1/Koslicki_lab_metagenomic_analysis/tree/cami_simu/3_mini_analysis/20210130_CAMISIM_errorfree_simulation)



### Step1: find genomes

---

Input strains can be found at [here](https://github.com/ShaopengLiu1/Koslicki_lab_metagenomic_analysis/blob/cami_simu/3_mini_analysis/20210130_CAMISIM_errorfree_simulation/input_species.txt)

- all 22/22 species were found, however, only 8 of matched strains were found.
- 4/22 records were not found in the downloaded GenBank annotation file: 2 of them were recorded by another name, 2 of them were found in the webportal (I guess the annotation file is not regularly updated)

```
1. Actinomyces odontolyticus ATCC 17982: equivalent "Schaalia odontolytica", GCA_000154225.1
2. Candida albicans ATCC MY-2876, GCA_003454735.1
3. Methanobrevibacter smithii ATCC 35061, GCA_000824705.1
4. Propionibacterium acnes DSM 16379, equivalent "Cutibacterium acnes KPA171202", GCA_000008345.1
```

- the genomic data can be found at:

```
/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210130_CAMISIM_errorfree_simulation/input_data/genome_data
```



### Step2: prepare CAMISIM environment

---

- Followed the [manual of CAMISIM](https://github.com/CAMI-challenge/CAMISIM/wiki/User-manual#File-Formats) (Installation from source, not Docker), I created a py27 conda environment with few modifications:
  1. use Samtools v1.9 instead of v1.0
  2. change size from default 0.1 to 0.5 (size of simulated sample in GB)
  3. change seed to 1234567
- Error profile
  1. by default, regular simulation (with sequence errors) is performed by [ART](https://academic.oup.com/bioinformatics/article/28/4/593/213322) which relies on well trained error profiles. That's the default setting in CAMISIM and no changes needed. See [here](https://github.com/ShaopengLiu1/Koslicki_lab_metagenomic_analysis/blob/cami_simu/3_mini_analysis/20210130_CAMISIM_errorfree_simulation/CAMISIM_config_file/regular_run/new_mini_config.ini) for the input for the regular run.
  2. to make an error-free simulation, CAMISIM relies on [wgsim](https://github.com/lh3/wgsim), another simulator wrote by Heng Li (the author of BWA). Unfortunately, CAMISIM is not well tuned for this purpose. But I found an official "guidance" in [issue53](https://github.com/CAMI-challenge/CAMISIM/issues/53) and followed it to modify the input settings. See [here](https://github.com/ShaopengLiu1/Koslicki_lab_metagenomic_analysis/blob/cd8af908036a8818c6edccced87d04d154480323/3_mini_analysis/20210130_CAMISIM_errorfree_simulation/CAMISIM_config_file/error_free_run/new_mini_config.ini#L14-L19) for the error-free run. 



### Step3: run CAMISIM

---

- Manually put config files into CAMISIM folder and use the de novo running:

```
python metagenomesimulation.py defaults/personal/new_mini_config.ini
```

- There were 10 independent communities (same species but different abundances) , output files can be found at. 

```
/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210130_CAMISIM_errorfree_simulation/output
```





### Note: unequal size

---

Though the size of the two fq.gz files are different, they did contain comparable sequencing depth. The size difference is due to compression (error-free one is more compressable because the quality line is just IIIIIIIIIIIIII).
