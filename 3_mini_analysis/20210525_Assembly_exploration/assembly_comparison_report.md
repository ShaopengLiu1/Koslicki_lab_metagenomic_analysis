# Assembly comparison for old and new data

Shaopeng Liu (sml6467@psu.edu)

05/25/2021

Local files: 

```bash
/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210525_Assembly_exploration
```

---

### Table of contents

- [1. Data filter](#1_filter)
- [2. Search for pairs](#2_pair)
- [3. Demo: locate data for a selected pair](#3_data)
- [4. Pipe: comparison](#4_pipe)





---

### Step1: data filter <a name="1_filter"></a> 

1.1 Download ref file

```bash
# metadata
wget -q -O NCBI_GenBank_bacteria_assembly_summary.txt  ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt
# documents for metadata
wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt
```



1.2 Quality filter

- Col11: "Version status": "latest" (the most recent of all the versions for this assembly chain)
- Col14: "genome_rep": "Full" (the data used to generate the assembly was obtained from the whole genome, as in Whole Genome Shotgun (WGS) assemblies for example. There may still be gaps in the assembly.)

```bash
# interestingly, the metadata ONLY has "latest" records......
#while "latest" refers to an assembly chain
# if you get sth here, that means it's a new assembly chain for the same species
cut -f 11 NCBI_GenBank_bacteria_assembly_summary.txt | sort | uniq -c
################################
 967224 latest
      1 #   See ftp://ftp.ncbi.nlm.nih.gov/genomes/README_assembly_summary.txt for a description of the columns in this file.
      1 version_status
################################

awk -F '\t'  '{if($14=="Full" && $11=="latest") print $0}' NCBI_GenBank_bacteria_assembly_summary.txt > filtered_metadata.txt
```



1.3 Use "contig" only because scaffolding is another (big) computational challenge

- Col12: "assembly_level": "Contig"

```bash
# get some distribution infor
cut -f 12 filtered_metadata.txt | sort | uniq -c

################################
   4881 Chromosome
  25551 Complete Genome
 798307 Contig
 138278 Scaffold
################################

# keep only contig data
awk -F '\t'  '{if($12=="Contig") print $0}' filtered_metadata.txt > contig_only_filtered_metadata.txt  
```





---

### Step2.1: search for species-level old/new pairs <a name="2_pair"></a> 

2.1.1 Species-level pairs (just in case)

- Col7: "species_taxid"

```bash
cut -f 7 contig_only_filtered_metadata.txt | sort | uniq -c | awk '$1>1' > species_level_hits.txt
wc -l species_level_hits.txt   # 5023 -> because there are less species
```



2.1.2 Check time-span

- Col 15: "seq_rel_date"

```bash
# originally set for 2 purposes, but then found that I can't direct search strain in this way
function get_date_range() {
  search_key=$1 #strain name OR species id for 2 levels
  match_col=$2 #strain name: col9, taxid: col7
  #echo "searching $search_key"
  awk -F "\t" -v input_key="$search_key" -v input_col="$match_col" '$input_col==input_key' contig_only_filtered_metadata.txt | cut -f 15 > _temp.txt
  # write into sec since 1970-1-1
  echo -n ""  > _time.txt
  for time_tag in $(cat _temp.txt); do
    date -d $time_tag "+%s" >> _time.txt
  done
  # get time range
  sort -k1,1n _time.txt  > _range.txt
  start=$(head -1 _range.txt)
  end=$(tail -1 _range.txt)
  range=$(echo "scale=1;($end - $start)/3600/24" | bc)  # round to day
  # return range value
  rm _time.txt _range.txt _temp.txt
  echo $range
}

# test function
awk '$7==100371' contig_only_filtered_metadata.txt | cut -f 15  # 16 records, date range: 2019.3.1-2020.7.20
get_date_range 100371 7  # 507 days
aa=$(get_date_range 100371 7)
echo $aa  # 507

# get date range for species data
echo "" > date_range_species.txt
for taxid in $(awk '{print $2}' species_level_hits.txt); do
  echo processing $taxid
  date_range=$(get_date_range "$taxid" 7)
  echo -e "$taxid\t$date_range" >> date_range_species.txt
done
```



### Step2.2: search for strain-level old/new pairs

2.2.1 Re-format the metadata for better search

- Col9: Infraspecific name: the strain, breed, cultivar or ecotype of the organism from which the sequences in the genome assembly were derived.
- Because 1 strain name can be shared by many species, e.g. strain=1. Need to paste species name and strain together to remove (many) false positive
- col1 is GCA identifier, col7 is taxid, col8 is sci name, col9 is strain name, col15 is release date

```bash
awk -F"\t" '$9' contig_only_filtered_metadata.txt | awk -F"\t" '{print $1"\t"$7"\t"$8" "$9"\t"$15}' > merge_strain_name_contig_only_filtered_metadata.txt
```



2.2.2 Strain-level pairs (ideal situation)

- Count unique on merged name+strain

```bash

cut -f 3 merge_strain_name_contig_only_filtered_metadata.txt | sort | uniq -c | awk '$1>1' | sed '1d' > strain_level_hits.txt
wc -l strain_level_hits.txt    # 3741 -> pretty much
```

- get date range

```bash
# tune the grep function
function get_strain_date_range() {
  search_key=$1 #strain name OR species id for 2 levels
  awk -F "\t" -v input_key="$search_key" '$3==input_key' merge_strain_name_contig_only_filtered_metadata.txt | cut -f 4 > _temp.txt
  # write into sec since 1970-1-1
  echo -n ""  > _time.txt
  for time_tag in $(cat _temp.txt); do
    date -d $time_tag "+%s" >> _time.txt
  done
  # get time range
  sort -k1,1n _time.txt  > _range.txt
  start=$(head -1 _range.txt)
  end=$(tail -1 _range.txt)
  range=$(echo "scale=1;($end - $start)/3600/24" | bc)  # round to day
  # return range value
  rm _time.txt _range.txt _temp.txt
  echo $range
}

# test function
grep -w "Acidovorax caeni strain=R-24608" merge_strain_name_contig_only_filtered_metadata.txt
get_strain_date_range "Acidovorax caeni strain=R-24608"  #415


### due to some space in the strain names, need to print all but $1, e.g. "ecotype=Quinolone resistant". This can't be properly handled by awk '{for(i=2,i<=NF,i++)}' OR cut -d" " -f 8-. 
### Notice all first 8 characters are not needed, so sed is the easiest solution
sed 's/^........//' strain_level_hits.txt > _temp_strain_name.txt
echo "" > date_range_strain.txt
while read strain_name; do
  echo processing $strain_name
  date_range=$(get_strain_date_range  "$strain_name")
  echo -e "$strain_name\t$date_range" >> date_range_strain.txt
done  < _temp_strain_name.txt
```



### Step2.3 Summary stats for available data

2.3.1 clean results

```bash
rm _temp_strain_name.txt 

sed '1d' date_range_species.txt | sort -k2,2rn  > _temp && mv _temp date_range_species.txt
sed '1d' date_range_strain.txt  | sort -t$'\t' -k2,2rn  > _temp && mv _temp date_range_strain.txt
```



2.3.2 summary stats

```python
# py script
import os
import pandas as pd
import seaborn as sb
from matplotlib import pyplot as plt
import math

#local dir
#os.chdir("/Users/shaopeng/Desktop/assembly_mini")

# read table and round up (cause bash is just trunc)
f_strain = pd.read_table("date_range_strain.txt", header=0)
f_strain.columns=['name', 'range']
f_species = pd.read_table("date_range_species.txt", header=None)
f_species.columns=['name', 'range']
f_species['range'] = f_species['range'].fillna(0)

# get date ranges
l1=list(f_strain['range'])
l_strain = [int(math.ceil(x)) for x in l1]
l2=list(f_species['range'])
l_species = [int(math.ceil(x)) for x in l2]

# plot
def plot_range_distribution(input_list, out_name="test"):
	min_v = min(input_list)
	max_v = max(input_list)
	print("Min is "+str(min_v))
	print("Max is "+str(max_v))
	plot_bins = list(range(min_v, max_v, 365))
	fig, axs = plt.subplots(1, 1, figsize=(16, 16))
	plt.rcParams.update({'font.size': 30})
	sb.histplot(input_list, ax=axs, bins=plot_bins, fill=True)
	axs.title.set_text("Histplot of date range by year")
	#axs.set_ylim([0, 1000])
	axs.set(ylabel="Count", xlabel="Date range")
	fig.savefig("Date_range_distri_" + out_name + ".png", dpi=200)
	plt.close(fig)
	
plot_range_distribution(l_species, "species")
plot_range_distribution(l_strain, "strain")
```



#### Date range distribution for strain-level pairs (bin size: 1 year).

- a significant portion of strains have data entries that span more than 1 year

<img src="Date_range_distri_strain.png" style="zoom: 15%">





---

### Step3. Locate data for a given pair <a name="3_data"></a> 

3.1 locate data

```bash
target="Staphylococcus aureus strain=A69"  # use this as an example
grep "$target" merge_strain_name_contig_only_filtered_metadata.txt | sort -t$'\t' -k4,4rV  
# Output:
# GCA_903804755.1	1280	Staphylococcus aureus strain=A69	2020/06/15
# GCA_000698045.1	1280	Staphylococcus aureus strain=A69	2014/06/04
```



3.2 Use the GCA identifier to find raw data

Search from the [NCBI Assembly database](https://www.ncbi.nlm.nih.gov/assembly/).



3.2.1 NCBI FTP only contains the assemblied data

The [FTP link](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/903/804/755/GCA_903804755.1_21673_5_85) contains processed data based on the Readme file. 

We can find the assembly stats here, but not sure if we need it.



3.2.2 Find the raw data

Raw data

- Search GCA: https://www.ncbi.nlm.nih.gov/assembly/GCA_010429105.1/
- Go to Biosample: https://www.ncbi.nlm.nih.gov/biosample/SAMN01823741/
- SRA: https://www.ncbi.nlm.nih.gov/sra/?term=SRS510247

Assembler: skip those without assmbler infor (not in NCBI or can't find in publication)

```
# target="Salmonella enterica subsp. enterica serovar Rissen str. 150 strain=150"
# GCA_010429105.1	28901	Salmonella enterica subsp. enterica serovar Rissen str. 150 strain=150	2020/02/11
# GCA_000335835.1	28901	Salmonella enterica subsp. enterica serovar Rissen str. 150 strain=150	2013/02/01

# new: SKESA v. 2.2 (bingo)
# old: Newbler v. 2.3   (https://www.biostars.org/p/260053/) not live anymore, need to contact Roche for a copy
```



3.3.3 Notes before the analysis:

- time span: crossing different years.

- fix old/new assemblers? Or be consistent with data?
  - old assembly on old data: download the uploaded files  (assembled fna.gz)
  - pick some state-of-art (currently the best)
  - Reach to Paul (Shiki?)  http://rayan.chikhi.name/



---

### Step4. Exploration plan (will build an auto-pipe)

4.1 Pick several data with time span within 1~7 years:

| Strain                                       | Last update | Range/yr | Newest data                                                  | Assembler         | Oldest data                                                  | Assembler          |
| -------------------------------------------- | ----------- | -------- | ------------------------------------------------------------ | ----------------- | ------------------------------------------------------------ | ------------------ |
| Escherichia coli strain=G5                   | 2021/04/16  | 6.5      | [GCA_018044965.1](https://www.ncbi.nlm.nih.gov/assembly/GCA_018044965.1/) | SKESA v. 2.2      | [GCA_000768485.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000768485.1/) | Velvet v. 1.2.10   |
| Aliivibrio fischeri strain=KB2B1             | 2021/04/13  | 4.9      | [GCA_017921855.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_017921855.1/) | HGAP v. 2019      | [GCA_001640305.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_001640305.1/) | CLC NGS Cell v. 6  |
| Staphylococcus haemolyticus strain=IIF2SW-P5 | 2020/08/22  | 3.9      | [GCA_014266945.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_014266945.1/) | SPAdes v. v3.11.1 | [GCA_001743425.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_001743425.1/) | A5 v. 20150522     |
| Klebsiella pneumoniae strain=70              | 2020/09/11  | 2.9      | [GCA_014526145.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_014526145.1/) | SPAdes v. 3.9.0   | [GCA_002411925.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_002411925.1/) | SPAdes v. 3.9.1    |
| Listeria monocytogenes strain=LM11           | 2019/12/29  | 2.2      | [GCA_009807375.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_009807375.1/) | SPAdes v. 3.10    | [GCA_002776275.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_002776275.1/) | SPAdes v. v3.9.0   |
| Escherichia coli strain=C127                 | 2020/05/11  | 1.5      | [GCA_013030815.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_013030815.1/) | SPAdes v. 3.12    | [GCA_003721935.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_003721935.1/) | SPAdes v. 3.12.0   |
| Clostridium perfringens strain=2C45          | 2020/12/21  | 0.5      | [GCA_016236225.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_016236225.1/) | SPAdes v. 3.9.0   | [GCA_013305215.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_013305215.1/) | unicycler v. 0.4.8 |



```
### code example:
target="Escherichia coli strain=A3"
grep -w "$target" merge_strain_name_contig_only_filtered_metadata.txt | sort -t$'\t' -k4,4rV  

```



4.2 Comparisons:

1. Method improvement: New assembler on old data vs old assembler on old data (directly download)
2. Data improvement: New assembler on new data (directly download) vs new assembler on old data
3. ***Might be better to process both side instead of downloading them for fair comparison



4.3 Assemblers to use:

May need some followup

| New                                   | Old               |
| ------------------------------------- | ----------------- |
| SKESA v. 2.2 (2018)                   | Velvet v. 1.2.10  |
| HGAP v. 2019                          | CLC NGS Cell v. 6 |
| SPAdes v. 3.9~3.12  (which is better) | SPAdes            |
|                                       |                   |
|                                       |                   |



4.4 Pipeline:

1. Take an arbitrary input genome and process by 5+5 assemblers 
2. Bridged to the comparison platform



More infor for ref:

1. https://assemblathon.org/

2. http://bioinf.spbau.ru/quast  

   
