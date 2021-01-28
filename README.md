# Koslicki_lab_metagenomic_analysis
For metagenomic-related analysis pipeline and resources  
Shaopeng Liu (sml6467@psu.edu)  

## **Table of Contents**
- [1. File introduction](#1_file_intro)
- [2. Analysis and progress](#2_analysis)
  - [2.1 CMash](#21_CMash_repro)
  - [2.2 ANI-based similarity aggregation](#22_ani)
- [3. Reproducibility](#3_reproducible)
  - [3.1 Install dependency (only once)](#31_depen)
  - [3.2 Activate working environment](#32_activate)
  - [3.3 Reproduce 2.1_CMash results](#33_CMash)
  - [3.4 Reproduce 2.2_similarity_tree results](#34_simi)
- [4. Change log and note](#4_log)
  
## 1. File introduction <a name="1_file_intro"></a>
- [src](https://github.com/ShaopengLiu1/Koslicki_lab_metagenomic_analysis/tree/main/src): a source folder to store necessary tools and scripts
   - install_dependency.sh: a bash script to download tools and to create a conda envrionment for downstream usage
- [CMash-reproducible](): reproducibility for all CMash-related analysis
- [ANI-tree](): reproducibility for all ANI-based similarity tree analysis

## 2. Analysis <a name="2_analysis"></a>
### 2.1 CMash <a name="21_CMash_repro"></a>
Implement truncation-based algorithm to accelerate k-mer based estimations.  
- [ ] generate a random sample of species under phylum/class/order/family/genus (though we only need 3, but not sure which three works best) and download the corresponding complete & (representative || random pick 1 if no representative genomes defined)
- [ ] create a py wrapper for CMash based on the previous CMash-reproducible repo to implement previous analysis
- [ ] pipeline previous steps for automatic run
- [ ] finish manuscript

### 2.2 Similarity tree <a name="22_ani"></a>
Try ANI-based aggregation tree for species clustering.  
- [ ] collect random samples from different taxonomic ranks
- [ ] use JI/CI/weighted overlap?(as MetaPlatte) to estimate the ANI
- [ ] cluster species based on ANI matrix
- [ ] find some methods to validate the performance of the clustering 

## 3. Reproducibility <a name="3_reproducible"></a>
Repeat analysis results by bash/py pipelines.
### 3.1 Install dependency (only need to be done once) <a name="31_depen"></a>
(1) clone this repo:
```
git clone https://github.com/ShaopengLiu1/Koslicki_lab_metagenomic_analysis.git
```
(2) go to src folder and install necessary tools
```
cd src
bash install_dependency.sh
```

### 3.2 Activate working environment (every time before running the analysis) <a name="32_activate"></a>

### 3.3 Reproducibles for CMash <a name="33_CMash"></a>

### 3.4 Reproducibles for ANI-similarity-based tree <a name="34_simi"></a>


## 4. Change log and note <a name="4_log"></a>

