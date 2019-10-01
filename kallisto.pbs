#!/bin/sh -f
# qsub kallisto.pbs

#PBS -l nodes=1:ppn=16
##PBS -l pmem=1gb
#PBS -A open
#PBS -l walltime=23:50:00
#PBS -j oe
#PBS -o kallisto.out
#PBS -N kallisto
#PBS -M juc326@psu.edu
#PBS -m abe
#######################################
# kallisto v44.0 installed with conda
module load gcc/5.3.1 
cd $PBS_O_WORKDIR

# build index based on assembled transcriptome from StringTie
kallisto index stringtie_merged.fa

mkdir kallisto_output
kallisto quant -t 20 -i stringtie_merged.idx -o kallisto_output/E1pc --single -l 200 -s 20 E1pc/data.fa
kallisto quant -t 20 -i stringtie_merged.idx -o kallisto_output/E2pc --single -l 200 -s 20 E2pc/data.fa
kallisto quant -t 20 -i stringtie_merged.idx -o kallisto_output/E3pc --single -l 200 -s 20 E3pc/data.fa

kallisto quant -t 20 -i stringtie_merged.idx -o kallisto_output/W1pc --single -l 200 -s 20 W1pc/data.fa
kallisto quant -t 20 -i stringtie_merged.idx -o kallisto_output/W2pc --single -l 200 -s 20 W2pc/data.fa
kallisto quant -t 20 -i stringtie_merged.idx -o kallisto_output/W3pc --single -l 200 -s 20 W3pc/data.fa

cd kallisto_output
# combine tpm count from all samples
paste E1pc/abundance.tsv E2pc/abundance.tsv E3pc/abundance.tsv W1pc/abundance.tsv W2pc/abundance.tsv W3pc/abundance.tsv temp.txt
cut -f 1, 5, 10, 15, 20, 25, 30 temp.txt > kallisto_count.txt

#######################################
echo "Job Started"
echo "Job Ended"
