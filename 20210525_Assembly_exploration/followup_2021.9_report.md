## Infor

Continuation of initial exploration analysis

More details check slides in "Personal_meeting" folder

Date:  09/06/2021

Local files (not in Git)

```bash
/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210525_Assembly_exploration/stage2_followup
```

---

### Analysis plan

1. Take some RefSeq genome and simulate reads from it (one of the species that were in your analysis you sent earlier)
2. run the oldest version of SPADES on it and the newest version
3. generate the QUAST report including the RefSeq genome as the “ground truth”
4. generate a histogram of the number of WGS isolate raw sequencing data from the SRA versus date of submission (or sample collection) (Perhaps there is some master metadata file somewhere in the SRA)

---

### Contents

- [1. Data and SPADES prep](#1_data)
- [2. Simulation](#2_simu)
- [3. QUAST report with groud truth](#3_quast)
- [4. SRA data infor](#4_sra)
- [5. Summary](#5_sum)



---

### Step1: Data and SPADES prep <a name="1_data"></a> 

1.1 RefSeq metadata

```
# download metadata
wget -O ref_seq_bac_metadata.txt  https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt

# check records
sed '1,2d' ref_seq_bac_metadata.txt | cut -f 11 | sort | uniq -c
### all latest: 223746 latest
sed '1,2d' ref_seq_bac_metadata.txt | cut -f 12 | sort | uniq -c
### keep Complete genomes only
   3792 Chromosome
  23159 Complete Genome
 118544 Contig
  78251 Scaffold

### simplify the list
awk -F '\t' -v OFS='\t'  '{if($12=="Complete Genome" && $11=="latest") print $1,$5,$6,$7,$8,$15,$20}' ref_seq_bac_metadata.txt  | cat <(echo -e "accession\tref_category\ttaxid\tspecies_taxid\tname\trel_date\tftp_path") - > simplified_refseq_record.txt
```



1.2 Input genome

```bash
species="Acinetobacter baumannii"

# number of records
grep -w "${species}" simplified_refseq_record.txt | awk '$3==$4' | wc -l 

# get the lastest
grep -w "${species}" simplified_refseq_record.txt | awk '$3==$4' | sort -t$'\t' -k6,6rV | head

```

| Species (strain)                          | Full genome | Records number | Shortcut in iniital report | Latest_Accession | Rel_date   |
| ----------------------------------------- | ----------- | -------------- | -------------------------- | ---------------- | ---------- |
| Acinetobacter baumannii MRSN15038         | species     | 254            | S1                         | GCF_019458485.1  | 2021/08/08 |
| Escherichia coli C305                     | species     | 1418           | S2                         | GCF_019843995.1  | 2021/09/01 |
| Staphylococcus aureus 17                  | species     | 564            | S3                         | GCF_019738975.1  | 2021/08/24 |
| Campylobacter jejuni PNUSAC009032         | species     | 165            | S4                         | GCF_019754215.1  | 2021/08/25 |
| Salmonella enterica PNUSAS029866          | species     | 87             | S6                         | GCF_019710375.1  | 2021/08/21 |
| Vibrio cholerae O1 Inaba El Tor           | species     | 72             | S7                         | GCF_019504425.1  | 2021/08/10 |
| Listeria monocytogenes FDA1077798-012-001 | species     | 189            | S8                         | GCF_019332005.1  | 2021/07/23 |



1.3 download all genome data

```bash
accession=(GCF_019458485.1 GCF_019843995.1 GCF_019738975.1 GCF_019754215.1 GCF_019710375.1 GCF_019504425.1 GCF_019332005.1)

for id in "${accession[@]}"; do 
	echo $id;
	ftp=$(grep ${id} simplified_refseq_record.txt | cut -f 7);
	suffix=$(echo ${ftp##*/});
	wget ${ftp}/${suffix}_genomic.fna.gz 
done

# rename file for easier recognition
cp GCF_019332005.1_ASM1933200v1_genomic.fna.gz S8.fna.gz
cp GCF_019458485.1_ASM1945848v1_genomic.fna.gz S1.fna.gz
cp GCF_019504425.1_ASM1950442v1_genomic.fna.gz S7.fna.gz
cp GCF_019710375.1_ASM1971037v1_genomic.fna.gz S6.fna.gz
cp GCF_019738975.1_ASM1973897v1_genomic.fna.gz S3.fna.gz
cp GCF_019754215.1_ASM1975421v1_genomic.fna.gz S4.fna.gz
cp GCF_019843995.1_ASM1984399v1_genomic.fna.gz S2.fna.gz
gunzip S*.fna.gz
```



1.4 clean folder

```bash
mkdir raw_data
mv GCF_*_genomic.fna.gz ./raw_data
mv *txt ./raw_data

mkdir input_fasta
mv S*.fna ./input_fasta

mkdir simulation
```



1.5 SPADES prep

| Software      | Version        | Category | Source                  | Local executive      |
| ------------- | -------------- | -------- | ----------------------- | -------------------- |
| SPAdes (2012) | 3.15.3 @2021.7 | New      | Conda (see init report) | Env "new_assmb_py37" |
| SPAdes (2012) | 3.7.0 @2016.3  | Old      | Conda (see init report) | Env "old_assmb_py35" |

Note:

1. Though published in 2012, the oldest release on Github is 3.8 (2016.6)
2. There are 2 earlier version in conda: https://anaconda.org/bioconda/spades/files?version=3.7.0
   1. 3.5 was upload in 2019 (no idea why) for py2.7. I didn't use this considering the time and py version
   2. 3.6.2 was uploaded in 2016.1, pretty close to 3.7. But it requires earlier version. That's why I use v3.7 in the initial report. 





### Step2: Simulation <a name="2_simu"></a> 

Don't need CAMISIM as this is not metagenomic condition.

Use "ART" for Illumina simulation.

Ref:

1. https://www.nature.com/articles/nrg.2016.57

Links:

1. Art conda: https://anaconda.org/bioconda/art
2. NIH list: https://www.niehs.nih.gov/research/resources/software/biostatistics/art/
3. Docu (better to use art -h): http://manpages.ubuntu.com/manpages/bionic/man1/art_illumina.1.html



Run simulation

```bash
conda activate metagenomic_py37  #that's my personal work env

cd simulation/
for file in $(ls ../input_fasta/); do
	art_illumina -ss HS25 -sam -i ../input_fasta/${file} -p -l 150 -f 200 -m 500 -s 10 -o simu_${file}
done

cd ..

# -ss HS25: HiSeq 2500
# -p: paired
# -l: length
# -f: coverage
# -m: fragment size
# -s: sd of fragment size
```





### Step3: SPAdes + QUAST <a name="3_quast"></a> 

Assembly by 2 versions of SPAdes

```bash
cd assembly_by_SPAdes/

# v3.15.3
conda activate new_assmb_py37
### add --idolate per recommendation, not available for v3.7
for file in $(ls ../input_fasta/); do
	echo $file
	spades.py --isolate -1 ../simulation/simu_${file}1.fq -2 ../simulation/simu_${file}2.fq -o SPAdes_New_${file}
done

conda deactivate



# v3.7.0
conda activate old_assmb_py35

for file in $(ls ../input_fasta/); do
	echo $file
	spades.py -1 ../simulation/simu_${file}1.fq -2 ../simulation/simu_${file}2.fq -o SPAdes_Old_${file}
done

conda deactivate

cd ..
```



Run Quast for all the results (add GT genome)

```bash
mkdir -p quast_out && cd quast_out

conda activate metagenomic_py37

#ignore: --ref-sam ../simulation/simu_${file}.sam
#seems skip BWA for ref_genome compared with -r, didn't find the document for the difference
#only add -r ref_genome
for file in $(ls ../input_fasta/); do
  echo "Processing $file"
  for version in New Old; do
    echo "Processing SPAdes_${version}_${file}"
    quast -o ./Quast_SPAdes_${version}_${file} -m 250 --circos --glimmer --rna-finding -1   ../simulation/simu_${file}1.fq -2 ../simulation/simu_${file}2.fq  -r ../input_fasta/${file}   ../assembly_by_SPAdes/SPAdes_${version}_${file}/contigs.fasta		
  done
done

conda deactivate
```



Collect all Quast results

```bash
# still in the quast_out folder

# row names:
find . -name "report.tsv" | head -1 | xargs cut -f 1 > merged_Quast_report.tsv

# need test
for folder in $(ls -d */); do
	echo $folder
	name=$(echo ${folder%.fna/})
	cat <(echo ${name})  <(cut -f 2 ${folder}/report.tsv | sed '1d') > _temp_data.txt
	paste merged_Quast_report.tsv _temp_data.txt > _temp_merge && rm _temp_data.txt
	mv _temp_merge merged_Quast_report.tsv
done
```





### Step4: SRA data information <a name="4_sra"></a> 

ref: 

1. An example: https://github.com/linsalrob/SRA_Metadata
2. SRA knowledge base: https://www.ncbi.nlm.nih.gov/books/NBK569234/#search.each_srx_entry_in_the_entrez_sr
   1. SRX: experiment
   2. SRR/ERR/DRR: 1 SRX may have 1 or many runs
   3. SRS: sample accession, imported from **BioSample**



```
cd /data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210525_Assembly_exploration/SRA_data_stats

# download SRA metadata
# 10GB, 53M records, may take several mins
wget ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab

file="SRA_Accessions.tab"
sed '1d' ${file} | cut -f 18 | sort | uniq -c
# total records: 53,536,353 
# unique biosamples: 12,592,974
sed '1d' ${file}  | cut -f 18 | sort | uniq -c | wc -l

### status distri
sed '1d' ${file} | cut -f 3 | sort | uniq -c

# 43813002 live
# 4451174 suppressed
# 5272176 unpublished


### data clean
# only keep {Accession, Received, BioSample}, and remove "-" records in BioSample
sed '1d' ${file}  | awk -F "\t" '$3=="live" && $18!="-" {print $1"\t"$6"\t"$18}'  > cleaned_SRA_accessions.txt
wc -l cleaned_SRA_accessions.txt
# 40599079 cleaned_SRA_accessions.txt
# dedup (keep oldest date)
sort -t$'\t' -k3,3V -k2,2V  cleaned_SRA_accessions.txt > temp && mv temp cleaned_SRA_accessions.txt
awk -F "\t" '!a[$3]++' cleaned_SRA_accessions.txt  > deduped_BioSample_list.txt
wc -l deduped_BioSample_list.txt
# 12592973
# this is a reasonable number based on ref1


### distribution collection
cut -f 2 deduped_BioSample_list.txt | cut -d"-" -f 1 | sort | uniq -c | awk '{print $2"\t"$1}' > distri_all_SRA.txt



### get microbiome only distri
### the GenBank/RefSeq metadata has no raw data infor, need to check with SRA metadata


sed '1,2d' ../stage2_followup/raw_data/ref_seq_bac_metadata.txt | cut -f 3 | sort -k1,1V | uniq -c | awk '{print $2}' > RefSeq_BioSample_records.txt
grep -w -f RefSeq_BioSample_records.txt  deduped_BioSample_list.txt > overlap_RefSeq_SRA.txt
wc -l overlap_RefSeq_SRA.txt
# 115225/223315
cut -f 2 overlap_RefSeq_SRA.txt  | cut -d"-" -f 1 | sort | uniq -c | awk '{print $2"\t"$1}' > distri_Ref_Seq.txt



sed '1,2d' ../NCBI_GenBank_bacteria_assembly_summary.txt | cut -f 3 | sort -k1,1V | uniq -c | awk '{print $2}' | sed '1d' >  GenBank_BioSample_records.txt
grep -w -f GenBank_BioSample_records.txt  deduped_BioSample_list.txt > overlap_GenBank_SRA.txt
wc -l overlap_GenBank_SRA.txt
# 773051/965953
cut -f 2 overlap_GenBank_SRA.txt  | cut -d"-" -f 1 | sort | uniq -c | awk '{print $2"\t"$1}' > distri_GenBank.txt

```



Note:

1. Time:

   1. Update time: col4
   2. Pub time: col5
   3. Received time (submission): col6
   4. No collection data, can be found (not guaranteed) inside the information page of the biosample

2. Identifier:

   1. Accession: col1, unique ID
   2. Submission: col2, 1 biosample may have MULTIPLE runs to submit, I'll ignore all dups. 1 for 1 biosample

3. Status: col3

   1. keep only "live"

4. Dup biosample:

   1. Keep the earlies date

5. BioSample: col18

   1. Unique identifier of the seq material, but may be corresponding to multiple SRR record

6. Adjustment to the table

   1. Some records do NOT include biosample information (a bad habit), but they usually belong to a series of which the biosample will be mentioned at least once. E.g: "DRA000006". These "-" records will be removed 
   2. Some records use numbers in 18th column (ummmmm), e.g. "DRR315880". No idea how they name their data, while I didn't find the corresponding BioSample, I'll count those numbers as BioSample records (keep original)
   3. Many microbiome BioSample can NOT be found in the SRA metadata (which is very strange)
      1. e.g. "SAMD00018573"
      2. Overlap rate:
         1. RefSeq: 115225/223315 = 51.6%
         2. GenBank: 773051/965953= 80.0%
      3. Based on the distri data: seems miss recent years (slower than other data submission). Haven't verified though.

   



Plot data in R

```R
setwd("./assembly_followup/")

dist_all <- read.csv("distri_all_SRA.txt", header=FALSE, sep="\t", col.names = c("Year", "Number_BS"))
dist_genbank <- read.csv("distri_GenBank.txt", header=FALSE, sep="\t", col.names = c("Year", "Number_BS"))
dist_refseq <- read.csv("distri_Ref_Seq.txt", header=FALSE, sep="\t", col.names = c("Year", "Number_BS"))

# separate box plot
single_bar <- function(input_df, name_key="temp") {
  png(paste0("Single_boxplot_", name_key, ".png"), width = 800, height = 600)
  barplot(input_df$Number_BS, main=name_key, xlab="Year", ylab="Number of BioSample", names.arg = input_df$Year)
  dev.off()
}

single_bar(dist_all, "All_SRA")
single_bar(dist_genbank, "GenBank")
single_bar(dist_refseq, "RefSeq")
```





### Step5: Summary <a name="5_sum"></a> 

Temp compare by R code

```R
quast <- read.csv("SPAdes_simu_compare_merged_Quast_report.tsv", sep="\t")
row.names(quast) <- quast[,1]
quast <- quast[,-1]
col_quast <- colnames(quast)

# rearrange pairs together
new_order=c()
for (i in seq(1,7)) {
  new_order = c(new_order, i, i+7)
}
pair_col <- col_quast[new_order]
pair_quast <- quast[, pair_col]
write.table(pair_quast, file="paired_SPAdes_simu_Quast_report.tsv", sep="\t", quote=FALSE, col.names=NA)


### Not used, just for backup
# def subset df function
get_subdf <- function(key_word, out_name="temp.tsv") {
  out_logic <- grepl(key_word, col_quast)
  out_df <- quast[,out_logic]
  write.table(out_df, file=out_name, sep="\t", quote=FALSE, col.names=NA)
  return(out_df)
}
```



Result overview:

BioSamples in SRA:

1. increasing trend, especially after 2018
2. RefSeq is more evenly distributed, possibly due to manually curation of input data
3. GenBank has lower missing records rate (in SRA database), possibly due to delayed update in 2021 (May file)
4. Most of data are later than 2015~2016



Assembly comparison:

1. N50: similar (slightly worse?) except S3
2. N75: similar except S2, S3
3. L50/L75: similar
4. Gene: similar
5. Conclusion: didn't see too much difference. Maybe due to good data quality? (Depth 200, previous 47~257, mean 121)
