## Infor

Automatic pipeline to compare SPAdes v3.15 and v3.07 based on arbitrary genome inputs

Date: 09/23/2021

[Code](https://github.com/ShaopengLiu1/Koslicki_lab_metagenomic_analysis/blob/main/3_mini_analysis/20210525_Assembly_exploration/stage3_pipe_for_comparison/src/auto_comparison_SPAdes_new_old_by_Accession.sh)

Local file:

```bash
/data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210525_Assembly_exploration/stage3_pipe_for_comparison
```



---

### Pipeline process:

1. Parameters (default):
   - coverage=100
   - length=150
   - fragment size=500
   - sysSeq=HS25
2. Download RefSeq genome by Accession ID
3. Simulate Illumina reads by "ART"

```bash
art_illumina -ss ${sysSeq} -sam -i ${out_name}_genomic.fna ${art_pair} -l ${length} -f ${coverage} -m ${fragment} -s 10 -o simu_${out_name}
```

4. Run SPAdes on simulated reads

```bash
# v3.15
conda activate new_assmb_py37
spades.py --isolate -1 simu_${out_name}1.fq -2 simu_${out_name}2.fq -o SPAdes_New_${out_name}

# v3.07
conda activate old_assmb_py35
spades.py -1 simu_${out_name}1.fq -2 simu_${out_name}2.fq -o SPAdes_Old_${out_name}
```

5. Quast

```bash
quast -o ./Quast_SPAdes_New_${out_name} -m 250 --circos --glimmer --rna-finding -1 simu_${out_name}1.fq -2 simu_${out_name}2.fq -r ${out_name}_genomic.fna  SPAdes_New_${out_name}/contigs.fasta

quast -o ./Quast_SPAdes_Old_${out_name} -m 250 --circos --glimmer --rna-finding -1 simu_${out_name}1.fq -2 simu_${out_name}2.fq -r ${out_name}_genomic.fna  SPAdes_Old_${out_name}/contigs.fasta
```



---

### Validation of code

(Files are in local folder)

1. analysis consistency:
   1. run S8 with coverage=200 (previous report), and the results are similar
2. parameter check: -isolate for SPAdes v3.15
   1. run S1 without -isolate parameter, with coverage=20, and the results are similar to current 



---

### New results

(local file)

1. previous S3 outperformance (N50) drop back
2. It looks like the SPAdes v3.07 slightly outperforms v3.15
   1. fewer contigs, higher N50
   2. similar coverage and recovery
   3. similar identified genes
3. Still over quality data???