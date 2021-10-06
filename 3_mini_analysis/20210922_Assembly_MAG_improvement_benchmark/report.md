## MAG improvement benchmark analysis

### Information

1. Local location:

   ```bash
   /data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark
   ```

2. Raw data:

   1. Data information: 
      1. CAMI2 challenge (simulated): extremely high strain diversity dataset (strain madness)
      2. CAMI DB: https://data.cami-challenge.org/participate (make an account)
         1. strain madness: https://www.biorxiv.org/content/10.1101/2021.07.12.451567v1.full.pdf
      3. FTP location: https://portal.nersc.gov/dna/RD/Adv-Seq/StrainMadness/
   2. Assembly results from MetaHipmer (copied, "final_assembly.fasta.gz"): `/data/temp/For_Shaopeng`

3. Data information:

   1. reads: 100 raw reads fq files
   2. final_assembly.fasta.gz: pool of multi-assembly. The depth file is abundance in 100 read files.
   3. assembly file: **needs confirm**
      1. pooled gold standard assembly: 35318 (records), long and complete, should be the contig from input genomes
      2. pooled megahit assembly: 1404890 (records), short contigs, the numbers make sense. This should be assembies from pooled data.
      3. final assembly.fa: 40603 pooled contigs, with average length ~3k. Looks like all contigs have been aggregated to reach such numbers? Possibly lose of resolution due to aggregation.

4. Refs:

   1. Metahipmer2 help: https://hpcadvisorycouncil.atlassian.net/wiki/spaces/HPCWORKS/pages/1827307543/MetaHipMer+2.0+for+ISC21#Build-and-Install-mhm2
   2. Install (No Github or Conda): https://bitbucket.org/berkeleylab/mhm2/src/master/
   3. https://bitbucket.org/berkeleylab/mhm2/src/master/docs/mhm_guide.md (check here)



---

### Analysis plan (to be piped)

1. Assemble metagenome into **metagenomic-assembled-genomes (MAGs)** by MetaHipmer from each datasets (original, new)
2. Aggregate MAGs from the datasets (pooled MAGs) and cluster them into **species genome bins (SGBs)** using 95% ANI cutoff (**FastANI**)
3. Map raw reads from **ALL** datasets to SGBs (BBTools), then repeat the assembly for each SGB (**STRONG**). We then get **read-derived-strain-bins**.
4. Merge the unitigs from MetaHipmer for each SGB, then re-bin them by MetaBAT2 to get **unitig-derived-strain-bins**, try 2 things:
   1. within original SGB
   2. ignore SGB and re-bin all contigs (shouldn't be too far from1)
5. Compare the quality of the following by AMBER (CAMI tool to access performance of binning):
   1. SGBs
   2. <u>read-derived-strain-bins</u>
   3. unitig-derived-strain-bins



### Datasets

1. CAMI2 strain-madness simulated dataset, some of the samples will be treated as "original", and the rest as "new"
2. Real-world datasets?



### To-be confirmed:

- [ ] Why need original vs new?
- [ ] Co-assembly (Rob Egan) vs multi-assembly (If so, MegaHit2 -> dedup/de-variation -> assembly)

---

### Benchmark interactive code

```bash
#Local file
cd /data/sml6467/github/Koslicki_lab_metagenomic_analysis/3_mini_analysis/20210922_Assembly_MAG_improvement_benchmark

# env preparation

```

