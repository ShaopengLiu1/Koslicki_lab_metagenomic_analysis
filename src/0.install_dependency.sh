# To install all necessary dependencies for reproducibility analysis
# Usage: bash <pipe>


# src folder
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $pipe_path
mkdir -p conda_env
mkdir -p git_repo


# install CAMISIM env
bash ${pipe_path}/CAMISIM_prepare_environment.sh


# install CMash env
bash ${pipe_path}/CMash_prepare_environment.sh


# download NCBI GenBank bacteria database
mkdir -p GenBank_db
cd GenBank_db
wget -O NCBI_GenBank_bacteria_assembly_summary.txt  ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
### note: this filter doesn't distinguish strain-level infor
### when multiple hits found, I usually use the 1st one (which is newer)
awk -F '\t'  '{if($12=="Complete Genome" && $11=="latest") print $6"\t"$8"\t"$20}' NCBI_GenBank_bacteria_assembly_summary.txt > NCBI_GenBank_download_link.txt



