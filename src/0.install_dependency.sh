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




