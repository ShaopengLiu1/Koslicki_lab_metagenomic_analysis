# To reproduce CMash results
# Usage: bash <pipe>



pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"



# generate random sample of species: phylum/family/genus x 3 replicates
bash ${pipe_path}/CMash_random_sample_generation.sh



# some bash wrapper to setup next step analysis
cd ${pipe_path}/../1_CMash/
time_tag=$(date +"%Y%m%d_%H")
mkdir -p CMash_out/result_${time_tag}
work_dir=${PWD}/CMash_out/result_${time_tag}
### prepare input file paths
cd taxon_records
for folder in $(ls -d sample_*); do
	realpath ${folder}/*.gz > full_path_${folder}.txt
done
mv full_path*.txt ${work_dir}
cd ${work_dir}



# call py script to do pairwise comparison
### activate conda env
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
        . ${conda_path}/etc/profile.d/conda.sh
else
        echo "ERROR: conda path can't be corrected identified!!!"
        exit 1
fi
conda activate ${pipe_path}/conda_env/CMash_env_py37

for i in 1 2 3; do
	echo "Processing sample $i"
	for input_file in $(ls full_path_sample_${i}_*txt); do
		python CMash_py_wrapper.py -q ${input_file} -r ${input_file}
		python CMash_plot_results.py 
		mkdir out_${input_file}
		mv *csv out_${input_file}
		mv *png out_${input_file}
	done
done



date
echo "pipe done"




