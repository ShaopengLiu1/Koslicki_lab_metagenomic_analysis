# interactive code to download genomes files for species in the input file
# in "input_data" folder

### download the NCBI GenBank bacteria database
wget -q -O NCBI_GenBank_bacteria_assembly_summary.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
### keep records for "complete"
awk -F '\t'  '{if($12=="Complete Genome") print $0}' NCBI_GenBank_bacteria_assembly_summary.txt > filtered_NCBI_assembly.txt
### search genus+species then strain (manual inspection needed)
echo -e "species-wise records\n\n" > check_genomes.txt
while read p; do
	genus=$(echo $p | cut -d" " -f 1)
	species=$(echo $p | cut -d" " -f 2)
	strain=$(echo $p | cut -d" " -f 4)
	echo $p >> check_genomes.txt
	grep $genus filtered_NCBI_assembly.txt | grep $species | grep $strain >> check_genomes.txt
	echo -e "\n\n" >> check_genomes.txt	
done < input_species.txt
### 6 of 22 have been found



### download criteria: if strain found, then 1st record of the matched strain (usually newer); if no strain match, just download the 1st appearance of the full genome (could be another strain though)
echo "Records for missing data" > missing_record.txt
while read p; do
        genus=$(echo $p | cut -d" " -f 1)
        species=$(echo $p | cut -d" " -f 2)
        strain=$(echo $p | cut -d" " -f 4)
	strain_hit=$(grep $genus filtered_NCBI_assembly.txt | grep $species | grep -c -w  $strain)
	if [ "$strain_hit" -gt 0 ]; then
		download_link=$(grep $genus filtered_NCBI_assembly.txt | grep $species | grep -w -m 1 $strain | cut -f 20)
		gca_id=$(echo ${download_link##*/})
		wget -q -O ${genus}_${species}_${strain}_${gca_id}_genomic.fna.gz ${download_link}/${gca_id}_genomic.fna.gz
	else
		#no strain match
		species_hit=$(grep $genus filtered_NCBI_assembly.txt | grep -c $species)
		if [ "$species_hit" -gt 0 ]; then
			download_link=$(grep $genus filtered_NCBI_assembly.txt | grep -m 1 $species | cut -f 20)
			gca_id=$(echo ${download_link##*/})
			wget -q -O ${genus}_${species}_${strain}_${gca_id}_genomic.fna.gz ${download_link}/${gca_id}_genomic.fna.gz
		else
			echo $p >> missing_record.txt
		fi
	fi
done < input_species.txt



### results: 6/22 strain match, 18/22 species match. 
### for the missing 4:
# Actinomyces odontolyticus ATCC 17982: equivalent "Schaalia odontolytica", GCA_000154225.1
# Candida albicans ATCC MY-2876, GCA_003454735.1
# Methanobrevibacter smithii ATCC 35061, GCA_000824705.1
# Propionibacterium acnes DSM 16379, equivalent "Cutibacterium acnes KPA171202", GCA_000008345.1








