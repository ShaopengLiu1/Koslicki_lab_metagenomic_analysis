# CSE566 Project: Light-version containment MinHash analysis
For reproducibility analysis of of the CSE566 projects  
Note: this is the demo python version, maybe I'll have Rust version in the future......ummmmm, maybe  


## Install dependencies (only run once)
1. Install conda or [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) (version >= 4.6)
2. Clone this repo
```
git clone https://github.com/ShaopengLiu1/Koslicki_lab_metagenomic_analysis.git
```
3. create a Conda env for required tools
```
cd Koslicki_lab_metagenomic_analysis/src
bash CMash_prepare_running_env_conda_py37.sh
```



## Reproducibility: repeat the results in the reports
```
cd Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210410_CSE566_project
nohup bash project_reproducibility.sh &  #it may take some time to run
```




## Extended usage: for customized containment JI matrix
```
cd Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210410_CSE566_project
# activate the conda env
conda activate ../../src/conda_env/CMash_env_py37/
# store the abs paths of your query/ref files
python generate_3_JIs.py -q <path_to_query> -r <path_to_ref> -x True -z True 
```




