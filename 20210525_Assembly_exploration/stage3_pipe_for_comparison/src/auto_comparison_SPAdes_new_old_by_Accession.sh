# This is an automatic pipe to compare SPAdes3.15 vs 3.7 on simulated reads from arbitrary input genome (recognized by RefSeq Accession)
# Usage: bash <pipe> -a <Accession> (optional parameter)


### Parameters
work_dir="/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210525_Assembly_exploration/stage3_pipe_for_comparison"
src_dir="/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210525_Assembly_exploration/stage3_pipe_for_comparison/src"
cur_dir="$PWD"


while getopts a:d:l:m:s:p:n:h  opts
do case "$opts" in
a) accession="$OPTARG";;    # Accession ID
d) coverage="$OPTARG";;    # default 100
l) length="$OPTARG";;    # default 150
m) fragment="$OPTARG";;    # default 500
s) sysSeq="$OPTARG";;    # default HS25
p) paired="$OPTARG";;    # default PE
n) out_name="$OPTARG";;    # default to use accession ID
h) echo "
Compare SPAdes 3.7 vs 3.15 measured by Quast
usage:  <pipe> -a <Accession_ID> (e.g. GCF_900128725.1)
Optional parameter: -d <coverage> -l <length> -m <fragment_size> -s <sysSeq> -p PE/SE
Options:    -d      coverage of seq data, default 100
            -l      read length, default 150
            -m      fragment size, default 500
            -s      default HS25
                    The name of Illumina sequencing system of the built-in profile used for simulation
                    GA1 - GenomeAnalyzer I (36bp,44bp), GA2 - GenomeAnalyzer II (50bp, 75bp)
                    HS10 - HiSeq 1000 (100bp),          HS20 - HiSeq 2000 (100bp),      HS25 - HiSeq 2500 (125bp, 150bp)
                    HSXn - HiSeqX PCR free (150bp),     HSXt - HiSeqX TruSeq (150bp),   MinS - MiniSeq TruSeq (50bp)
                    MSv1 - MiSeq v1 (250bp),            MSv3 - MiSeq v3 (250bp),        NS50 - NextSeq500 v2 (75bp)
            -p      only supports PE now, needs to tune SPAdes code for SE.
	    -n      Output folder name, default to use the accession ID
"
exit;;
esac
done


date
echo "Pipe start......"

if [ -z "$accession" ]
then
    echo "Please specifying the Accession ID of input genome!!!"
    echo "Aborted......"
    exit 1
else
    ftp=$(grep ${accession} ${work_dir}/simplified_refseq_record.txt | cut -f 7)
    if [ -z "$ftp" ]; then
        echo "Didn't find the FTP link of input Accession ID, please double check!!!"
        exit 1
    fi
    echo -e "\n\n"
    echo "Processing $accession"
    echo "The corresponding FTP link is $ftp"
    echo -e "\n\n"
fi

if [ -z "$coverage" ];then
    coverage=100
    echo -e "Using coverage of $coverage \n"
fi

if [ -z "$length" ];then
    length=150
    echo -e "Using length of $length \n"
fi

if [ -z "$fragment" ];then
    fragment=500
    echo -e "Using fragment of $fragment \n"
fi

if [ -z "$sysSeq" ];then
    sysSeq=HS25
    echo -e "Using sysSeq of $sysSeq \n"
fi

if [ -z "$paired" ];then
    paired=PE
    echo -e "Using paired of $paired"
fi
if [[ $paired == PE ]]; then
    art_pair="-p"
    echo "Art parameter pair is ${art_pair}"
elif [[ $paired == SE ]]; then
    art_pair=""
    echo "Don't support SE now"
    exit 1
else
    echo "Wrong -p input, only PE or SE accepted"
    exit 1
fi

if [ -z "$out_name" ];then
    out_name=$accession
    echo -e "\nUsing out_name of $out_name \n"
fi


# use time cmd for some details
ltime="/usr/bin/time -av -o temp_runLog"

# to activate conda in script
temp=$(which conda)
conda_path=$(echo ${temp%/*bin/conda})
. ${conda_path}/etc/profile.d/conda.sh



### folder prep
time_tag=`date +"%y_%m_%d_%H-%M"`
mkdir -p ${work_dir}/output/${out_name}_depth_${coverage}_length_${length}_sysSeq_${sysSeq}_Time_${time_tag}
out_dir="${work_dir}/output/${out_name}_depth_${coverage}_length_${length}_sysSeq_${sysSeq}_Time_${time_tag}"
cd $out_dir
echo "Now working under $PWD......"
echo ""



### download data
echo -e "Downloading data...... \n"
suffix=$(echo ${ftp##*/})
wget ${ftp}/${suffix}_genomic.fna.gz
# make name consistent with folder anme
mv ${suffix}_genomic.fna.gz  ${out_name}_genomic.fna.gz
gunzip ${out_name}_genomic.fna.gz




### simulation
echo -e "Start simulation by ART...... \n"
conda activate metagenomic_py37
${ltime} art_illumina -ss ${sysSeq} -sam -i ${out_name}_genomic.fna ${art_pair} -l ${length} -f ${coverage} -m ${fragment} -s 10 -o simu_${out_name}
mv temp_runLog art_illumina.log
conda deactivate




### Run 2 versions of SPAdes
echo -e "Running SPAdes v3.15...... \n"
conda activate new_assmb_py37
${ltime} spades.py --isolate -1 simu_${out_name}1.fq -2 simu_${out_name}2.fq -o SPAdes_New_${out_name}
mv temp_runLog SPAdes_3.15.log
conda deactivate

echo -e "Running SPAdes v3.07...... \n"
conda activate old_assmb_py35
${ltime} spades.py -1 simu_${out_name}1.fq -2 simu_${out_name}2.fq -o SPAdes_Old_${out_name}
mv temp_runLog SPAdes_3.07.log
conda deactivate




### Run QUAST
conda activate metagenomic_py37
echo -e "Running QUAST to collect results...... \n"

${ltime} quast -o ./Quast_SPAdes_New_${out_name} -m 250 --circos --glimmer --rna-finding -1 simu_${out_name}1.fq -2 simu_${out_name}2.fq -r ${out_name}_genomic.fna  SPAdes_New_${out_name}/contigs.fasta
mv temp_runLog quast_SPAdes_New.log


${ltime} quast -o ./Quast_SPAdes_Old_${out_name} -m 250 --circos --glimmer --rna-finding -1 simu_${out_name}1.fq -2 simu_${out_name}2.fq -r ${out_name}_genomic.fna  SPAdes_Old_${out_name}/contigs.fasta
mv temp_runLog quast_SPAdes_Old.log

conda deactivate

#merge results
find . -name "report.tsv" | head -1 | xargs cut -f 1 > merged_Quast_report.tsv
for folder in $(ls -d Quast_SPAdes*/); do
	echo $folder
	name=$(echo ${folder%.fna/})
	cat <(echo ${name})  <(cut -f 2 ${folder}/report.tsv | sed '1d') > _temp_data.txt
	paste merged_Quast_report.tsv _temp_data.txt > _temp_merge && rm _temp_data.txt
	mv _temp_merge merged_Quast_report.tsv
done
mv merged_Quast_report.tsv merged_Quast_report_${out_name}.tsv


cd ${cur_dir}
echo "Pipe done!"
date


