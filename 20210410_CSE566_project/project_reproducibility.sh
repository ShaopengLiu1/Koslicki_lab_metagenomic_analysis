# to reproduce all results (except figures) for the class project

# Usage: bash <pipe>
threads=$1


date
### local variables
[ -z $threads ] && threads=1
ltime="/usr/bin/time -av -o temp_runLog"
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd ${pipe_path}



### activate conda
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
        . ${conda_path}/etc/profile.d/conda.sh
else
        echo "ERROR: conda path can't be corrected identified!!!"
        exit 1
fi
unset conda_path
conda activate ../../src/conda_env/CMash_env_py37/




### 1. download data
mkdir -p  genome_data
cd genome_data
mkdir -p Brucella_30_genomes
cd Brucella_30_genomes
for file in $(cut -f 3 ../../Brucella_genus_30_genomes.txt); do
	suffix=$(echo ${file##*/})
	wget -q ${file}/${suffix}_genomic.fna.gz 2>/dev/null
done
readlink -f GCA*_genomic.fna.gz > ../../fullpath_Brucella_30.txt

cd ${pipe_path}
cd genome_data
mkdir -p random_30_genomes
cd random_30_genomes
for file in $(cut -f 3 ../../top30_of_random_1k_genomes.txt); do
	suffix=$(echo ${file##*/})
	wget -q ${file}/${suffix}_genomic.fna.gz 2>/dev/null
done
readlink -f GCA*_genomic.fna.gz > ../../fullpath_random_1000.txt

cd ${pipe_path}





### 2. create artificial "merged" data for unbalanced settings
for file in $(head -2 fullpath_Brucella_30.txt); do
	cat $file >> f1.fa.gz
done

for file in $(head -3 fullpath_random_1000.txt); do
	cat $file >> f2.fa.gz
done

for file in $(head -8 fullpath_random_1000.txt); do
	cat $file >> f3.fa.gz
done

for file in $(head -18 fullpath_random_1000.txt); do
	cat $file >> f4.fa.gz
done

mkdir -p genome_data/merged_data
mv f*.fa.gz genome_data/merged_data
cd genome_data/merged_data
cat f1.fa.gz f2.fa.gz > merged_f12.fa.gz
cat f1.fa.gz f3.fa.gz > merged_f13.fa.gz
cat f1.fa.gz f4.fa.gz > merged_f14.fa.gz
readlink -f merged_f*fa.gz > ../../fullpath_ub_input.txt
cd ${pipe_path}





### 3. run pairwise JI for Genus Brucella: compare JI for |A| = |B|
mkdir -p out1_similar_size_output
cp fullpath_Brucella_30.txt ./out1_similar_size_output
cd out1_similar_size_output
# run JI for multiple k values
for k in $(seq 15 5 60); do
	echo "Running for k ${k}"
	python ../non_parallel_generate_3_JIs.py -q fullpath_Brucella_30.txt -r fullpath_Brucella_30.txt -t $threads -k ${k}
done
# then there are est_JI, GT_JI, containment_JI for each k values




### 4. run 2 JI methods in unbalanced settings
cd ${pipe_path}
mkdir -p out2_unbalance_setting/
cp fullpath_Brucella_30.txt ./out2_unbalance_setting/
cp fullpath_ub_input.txt ./out2_unbalance_setting/
cd out2_unbalance_setting/
# generate GT_JI (only need to do once) for different sketch size choices
python ../non_parallel_generate_3_JIs.py -q fullpath_Brucella_30.txt -r fullpath_ub_input.txt -x True -y True  
# generate est_JI, containment_JI for various m choices
for size in 500 1000 2000 5000 10000; do
	echo "Using size ${size} for sketch"
	python ../non_parallel_generate_3_JIs.py -q fullpath_Brucella_30.txt -r fullpath_ub_input.txt -z True -n ${size} -l m${size}
done
cd ${pipe_path}




### 5. generate temp plots
# however, those are primary plots so the code are very temporary
python plot_cse566_project_results.py





echo "pipe done"
date



































