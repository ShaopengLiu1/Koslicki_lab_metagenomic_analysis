# Aim: generate a random sample of species for CMash analysis
# Input: a seed
# Output: group of species under same taxon rank (from phylum to genus), group size ~30, replicates 3. Genomic files would be downloaded for future usage.


# need to activate conda env in the main bash script

# dependent variables
pipe_path="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ltime="/usr/bin/time -av -o temp_runLog"
taxon="taxonkit --data-dir ${pipe_path}/conda_env/CMash_env_py37/taxonkit_ncbi"
Genbank_ref="${pipe_path}/GenBank_db/NCBI_GenBank_download_link.txt"



# read parameter
seed=$1
size=$2
[ -z $seed ] && seed=1234567
[ -z $size ] && size=30

echo "pipe start......"
echo "The seed for random sampling is $seed"
echo "The size of random sample is $size"


# pipe start
cd ${pipe_path}/../1_CMash
work_dir=${PWD}


# taxon infor
mkdir -p taxon_records
cd taxon_records

### list all categories
${taxon} list --ids 2 --show-name --show-rank > all_bacteria.txt

### modify Genbank file for exact id match
awk -F"\t"  '{print "_"$1"_""\t"$2"\t"$3}' ${Genbank_ref} > _temp_genbank_ref.txt

### get all records
grep '\[phylum\]' all_bacteria.txt | sed 's/^[ \t]*//g' > all_phylum_record.txt
grep '\[family\]' all_bacteria.txt | sed 's/^[ \t]*//g' > all_family_record.txt
grep '\[genus\]' all_bacteria.txt | sed 's/^[ \t]*//g' > all_genus_record.txt

### for all taxon, check number of species with complete genomes in GenBank
check_species_num()
{
  input_file="$1"
  size="$2"
  echo -e "taxon_id\tfound_species" > species_num_${input_file}  #this is a record file
  mkdir output_${input_file}
  while read p; do
	  taxid=$(echo $p | cut -d" " -f 1)
	  ${taxon} list --ids ${taxid} --show-name --show-rank | grep '\[species\]'  > temp_all_species.txt
	  awk '{print "_"$1"_"}' temp_all_species.txt > temp_ids.txt
	  echo -n > temp_found.txt
	  grep -w -f temp_ids.txt _temp_genbank_ref.txt | sort -k1,1V -u >> temp_found.txt
	  hits_num=$(cat temp_found.txt | wc -l)
	  echo -e "${taxid}\t${hits_num}" >> species_num_${input_file}
	  echo "${taxid} has ${hits_num} species found......"
	  if [ $hits_num -gt $size ]; then
		  mv temp_found.txt ././output_${input_file}/species_under_${taxid}.txt
	  fi
	  rm temp_*.txt
  done < $input_file
}

### species counting
check_species_num all_phylum_record.txt ${size}
check_species_num all_family_record.txt ${size}
check_species_num all_genus_record.txt ${size}

### filter count file and then do a random selection
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

### random select 3 taxon for each level
generate_random_sample()
{
  input_file="$1"
  out_name="$2"  #phylum / family / genus
  # random sample of 3 taxon
  sed '1d' ${input_file}  | awk -v size=$size '$2>size' | \
	  shuf --random-source=<(get_seeded_random ${seed}) | \
	  head -3 > random_${out_name}_siez${size}_seed${seed}.txt
  # download samples
  for i in 1 2 3; do
	  taxon_id=$(sed -n "${i}p" random_${out_name}_siez${size}_seed${seed}.txt | cut -f 1)
	  mkdir sample_${i}_${out_name}_${taxon_id}
	  cat ./output_all_${out_name}_record.txt/species_under_${taxon_id}.txt | \
		  shuf --random-source=<(get_seeded_random ${seed}) | \
		  head -${size} | cut -f 3 > ./sample_${i}_${out_name}_${taxon_id}/download_links.txt
	  cd sample_${i}_${out_name}_${taxon_id}
	  for file in $(cat download_links.txt); do
		suffix=$(echo ${file##*/})
		echo "downloading $suffix"
		wget -q ${file}/${suffix}_genomic.fna.gz 2>/dev/null
	  done
	  rm download_links.txt
	  cd ..
  done
}

generate_random_sample species_num_all_phylum_record.txt phylum
generate_random_sample species_num_all_family_record.txt family
generate_random_sample species_num_all_genus_record.txt genus

### may need some cleaning
date
echo "Pipe done!"












