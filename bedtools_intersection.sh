#!/bin/sh -f

# this calculate the overlapping nucleotides between each sRNA locus and genome feature
for i in sRNA_loci/*.bed; do
    #echo $i
    for j in genome_features/*.bed; do
        #echo $j
        bedtools intersect -a $i -b $j | awk '{sum+=$3-$2} END {print sum}' -
    done
done

# this calculate the total nucleotides coverage of each type of sRNA locus 
for i in sRNA_loci/*.bed; do
    awk '{sum+=$3-$2} END {print sum}' $i
done

# this calculate the total nucleotides coverage of each type of genome_features
for i in genome_features/*.bed; do
    awk '{sum+=$3-$2} END {print sum}' $i
done