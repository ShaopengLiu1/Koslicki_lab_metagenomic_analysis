date

# for data preparation
### download the 100 raw reads
out_path="/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark/data/raw_reads_100"

for i in {0..99}; do
	echo "downlaoding seq data $i"
	wget -nv -P ${out_path} https://portal.nersc.gov/dna/RD/Adv-Seq/StrainMadness/reads/strmgCAMI2_long_read_sample_${i}_reads.fq.gz
done


### download other data
mkdir -p /data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark/data/raw_reads_100
mkdir -p /data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark/data/assembly


### assembly output file
wget -nv -P /data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark/data https://portal.nersc.gov/dna/RD/Adv-Seq/StrainMadness/final_assembly.fasta.gz
wget -nv -P /data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark/data  https://portal.nersc.gov/dna/RD/Adv-Seq/StrainMadness/final_assembly_depths.txt.gz


### gold standard?
wget -nv -P /data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark/data/assembly  https://portal.nersc.gov/dna/RD/Adv-Seq/StrainMadness/assembly/strmgCAMI2_short_read_pooled_gold_standard_assembly.fasta.gz
wget -nv -P /data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark/data/assembly  https://portal.nersc.gov/dna/RD/Adv-Seq/StrainMadness/assembly/strmgCAMI2_short_read_pooled_megahit_assembly.fasta.gz




