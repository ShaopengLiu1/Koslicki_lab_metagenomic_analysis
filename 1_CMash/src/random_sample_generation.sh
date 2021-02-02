# Note: this pipe has a dependent: taxonkit (not covered in install_dependency.sh yet)

# Aim: generate a random sample of species for CMash analysis
# Input: a seed
# Output: group of species under same taxon rank (from phylum to genus), group size ~30, replicates 3. Genomic files would be downloaded for future usage.

seed=$1
[ -z $seed ] && seed=1234567











