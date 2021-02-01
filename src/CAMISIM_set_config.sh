# to setup config file for a run based on selected data
# input file: TSV file w/o header and 1st column = NCBI_ID, 2nd column = absolute path
# output: new_genome_to_id.tsv, new_metadata.tsv. Two required input file for CAMISIM tailored by input information
# Usage: bash <pipe> <input_file> <cami_config_template> <yes/na for error free>


input_file=$1
config_temp=$2
[ -z $input_file ] && echo "Missing input file!!!" && exit 1
[ -z $config_temp ] && echo "Missing input file!!!" && exit 1

# genome_2_id file
cut -f 2 $input_file | awk '{print "genome"NR"\t"$1}' > new_genome_to_id.tsv


# metadata
### OTU is only related to a drawning limit (default 1 per OTU, changed to 100 per OTU as we really don't care it)
echo -e "genome_ID\tOTU\tNCBI_ID\tnovelty_category" > new_metadata.tsv
cut -f 1 $input_file | awk '{print "genome"NR"\t"NR"\t"$1"\t""Known_strain"}' > temp.txt
cat new_metadata.tsv temp.txt > fss.txt && mv fss.txt new_metadata.tsv && rm temp.txt


# edit config file
# assume those files would be put in the same folder as the config file
cp $config_temp ./new_mini_config.ini
### relative_path points to the position folder in CAMISIM repo where we will put those config files
relative_path=$(echo ${config_temp#*CAMISIM/})
relative_path=$(echo ${relative_path%/*})
### edit config file
sed -i "s#defaults/metadata.tsv#${relative_path}/new_metadata.tsv#g" new_mini_config.ini
sed -i "s#defaults/genome_to_id.tsv#${relative_path}/new_genome_to_id.tsv#g" new_mini_config.ini
total_genome=$(cat new_genome_to_id.tsv | wc -l)
sed -i "s/genomes_total=24/genomes_total=${total_genome}/g" new_mini_config.ini
sed -i "s/genomes_real=24/genomes_real=${total_genome}/g" new_mini_config.ini



# error free mode
error_free=$3
if [[ ! -z $error_free ]]; then
	if [ $error_free == "yes" ]; then
		echo "Setting config into error-free mode!!!!!!"
		sed -i "s#readsim=tools/art_illumina-2.3.6/art_illumina#readsim=tools/wgsim/wgsim#g" new_mini_config.ini
		sed -i "s#error_profiles=tools/art_illumina-2.3.6/profiles#error_profiles=.#g" new_mini_config.ini
		sed -i "s#profile=mbarc#profiles=0.0#g" new_mini_config.ini
		sed -i "s#type=art#type=wgsim#g" new_mini_config.ini
	fi
fi




# sum up
echo "please edit seed, size if necessary"
echo "please move generated files to parent folder of $config_temp"
echo "pipe done"


