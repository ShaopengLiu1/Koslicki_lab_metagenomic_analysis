# Koslicki_lab_metagenomic_analysis
For metagenomic-related analysis pipeline and resources  
Shaopeng Liu (sml6467@psu.edu)  

## **Table of Contents**
- [1. File introduction](#1_file_intro)
- [2. Analysis and progress](#2_analysis)
  - [1. CMash](#21_CMash_repro)
  - [2. ANI-based similarity aggregation tree](#22_ani)
  - [3. Small_analysis](#23_small)
- [3. Reproducibility](#3_reproducible)
  - [3.1 Install dependency (only once)](#31_depen)
  - [3.2 Reproduce 2.1_CMash results](#32_CMash)
  - [3.3 Reproduce 2.2_similarity_tree results](#33_simi)
- [4. Change log and note](#4_log)
  
## 1. File introduction <a name="1_file_intro"></a>
- [src](https://github.com/ShaopengLiu1/Koslicki_lab_metagenomic_analysis/tree/main/src): a source folder to store necessary tools and scripts
   - install_dependency.sh: a bash script to download tools and to create a conda envrionment for downstream usage
- [CMash-reproducible](): reproducibility for all CMash-related analysis
- [ANI-tree](): reproducibility for all ANI-based similarity tree analysis

## 2. Analysis <a name="2_analysis"></a>
### 1. CMash <a name="21_CMash_repro"></a>
Implement truncation-based algorithm to accelerate k-mer based estimations.  
- [x] generate a random sample of species under phylum/class/order/family/genus (though we only need 3, but not sure which three works best) and download the corresponding complete & (representative || random pick 1 if no representative genomes defined)
- [ ] (ongoing) **create a py wrapper for CMash based on the previous CMash-reproducible repo to implement previous analysis**
- [x] pipeline wrapper for previous steps for automatic run
- [ ] finish manuscript

### 2. Similarity tree <a name="22_ani"></a>
Try ANI-based aggregation tree for species clustering.  
- [ ] collect random samples from different taxonomic ranks
- [ ] use JI/CI/weighted overlap?(as MetaPlatte) to estimate the ANI
- [ ] cluster species based on ANI matrix
- [ ] find some methods to validate the performance of the clustering 

### 3. Small analysis <a name="23_small"></a>
Mini and usually independent analysis.
- [ ] Error free simulation by CAMISIM

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
bash 0.install_all_required_dependency_run_once.sh
```

### 3.2 Reproducibles for CMash <a name="32_CMash"></a>
(1) prepare input genomes
```
```


### 3.3 Reproducibles for ANI-similarity-based tree <a name="33_simi"></a>


## 4. Change log and note <a name="4_log"></a>

