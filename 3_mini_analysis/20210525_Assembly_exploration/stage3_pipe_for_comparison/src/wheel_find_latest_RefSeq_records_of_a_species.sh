### This is an example:
species="Acinetobacter baumannii"

# number of records
grep -w "${species}" simplified_refseq_record.txt | awk '$3==$4' | wc -l

# get the lastest
grep -w "${species}" simplified_refseq_record.txt | awk '$3==$4' | sort -t$'\t' -k6,6rV | head

