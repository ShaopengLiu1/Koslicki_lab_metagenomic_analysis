date

# for data preparation
### download the 100 raw reads
out_path="/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark/data/raw_reads_100"

for i in {0..99}; do
	echo "downlaoding seq data $i"
	wget -nv -P ${out_path} https://portal.nersc.gov/dna/RD/Adv-Seq/StrainMadness/reads/strmgCAMI2_long_read_sample_${i}_reads.fq.gz
done







