# select sample: ls -ltr | sed -n '1,30p' | awk '{print $9}'

date

echo "pipe start"

# active conda
pipe_path=/data/sml6467/github/Koslicki_lab_metagenomic_analysis/src

# activate conda
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
        . ${conda_path}/etc/profile.d/conda.sh
else
        echo "ERROR: conda path can't be corrected identified!!!"
        exit 1
fi
conda activate /data/sml6467/github/Koslicki_lab_metagenomic_analysis/src/conda_env/CMash_env_py37


get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}


work_dir=$PWD


for file in `ls species_under_*.txt`
do
	#download files
	mkdir output_${file}
	name=${file}
	cut -f 3 ${file} | shuf --random-source=<(get_seeded_random 1234567) | head -30 > output_${file}/accession.txt
	cd output_${file}
	for file in $(cat accession.txt);
	do
		suffix=$(echo ${file##*/})
		echo "downloading $suffix"
		wget -q ${file}/${suffix}_genomic.fna.gz 2>/dev/null
	done
	rm accession.txt
	realpath *gz > input_path.txt

	#run est_JI
	python ${pipe_path}/CMash_python_wrapper_of_CMash_github_to_perform_prefix_tree_bias_analysis.py -q input_path.txt -r input_path.txt -z True -y True -o True
	python /data/sml6467/github/Koslicki_lab_metagenomic_analysis/1_CMash/taxon_records/est_ji_scan_genus/temp_py_plot_est_ji.py 
	cp Boxplot_of_est_JI.png ../Boxplot_est_JI_${name}.png
	cd ${work_dir}
done


echo "pipe done"
date
	


