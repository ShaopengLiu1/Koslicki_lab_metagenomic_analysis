# This is INTERACTIVE code that only works for lab server
# dependencies have been installed, this is to record my codes for future review
date
echo "This is interactive scripts!!!"
exit 1


### setup var
work_dir="/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210525_Assembly_exploration/stage3_pipe_for_comparison"
src_dir="/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210525_Assembly_exploration/stage3_pipe_for_comparison/src"


### RefSeq download link (use latest records)
cd ${work_dir}
wget -O ref_seq_bac_metadata.txt  https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
awk -F '\t' -v OFS='\t'  '{if($12=="Complete Genome" && $11=="latest") print $1,$5,$6,$7,$8,$15,$20}' ref_seq_bac_metadata.txt  | cat <(echo -e "accession\tref_category\ttaxid\tspecies_taxid\tname\trel_date\tftp_path") - > simplified_refseq_record.txt



### conda env for tools
conda create -y -n old_assmb_py35 python=3.5
conda activate old_assmb_py35
conda install -y -c bioconda spades=3.7.0
conda deactivate

conda create -y -n new_assmb_py37 python=3.7
conda activate new_assmb_py37
conda install -y -c bioconda spades=3.15.3
conda deactivate

conda activate metagenomic_py37 #my own metagenomic ENV
conda install -y -c bioconda quast
conda install -y  -c bioconda art
conda deactivate



