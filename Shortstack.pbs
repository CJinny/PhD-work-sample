#!/bin/sh -f
# qsub Shortstack.pbs

#PBS -l nodes=1:ppn=4
##PBS -l pmem=1gb
#PBS -A open
#PBS -l walltime=23:50:00
#PBS -j oe
#PBS -o shortstack.out
#PBS -N shortstack
#PBS -M juc326@psu.edu
#PBS -m abe
#######################################

export PATH=~/work/shortstack:$PATH
cd $PBS_O_WORKDIR

# align_only to save merged_alignments.bam so that different custom locifiles can be used
ShortStack --readfile U-E_1.fa U-E_2.fa U-E_3.fa U-S_1.fa U-S_2.fa U-S_3.fa WT_1.fa WT_2.fa WT_3.fa --genomefile ~/work/Zea_mays.AGPv4.fa -bowtie_cores 4 --align_only

# de novo identification of siRNA and miRNA
ShortStack --bamfile merged_alignments.bam --genomefile ~/work/Zea_mays.AGPv4.fa -bowtie_cores 4

# custom locifiles
ShortStack --bamfile merged_alignments.bam --locifile class_I_TE.txt --genomefile ~/work/Zea_mays.AGPv4.fa -bowtie_cores 4
ShortStack --bamfile merged_alignments.bam --locifile class_II_TE.txt --genomefile ~/work/Zea_mays.AGPv4.fa -bowtie_cores 4
ShortStack --bamfile merged_alignments.bam --locifile lincRNA.txt --genomefile ~/work/Zea_mays.AGPv4.fa -bowtie_cores 4
ShortStack --bamfile merged_alignments.bam --locifile up2k_genes.txt --genomefile ~/work/Zea_mays.AGPv4.fa -bowtie_cores 4
ShortStack --bamfile merged_alignments.bam --locifile dn2k_genes.txt --genomefile ~/work/Zea_mays.AGPv4.fa -bowtie_cores 4


#######################################
echo "Job Started"
echo "Job Ended"
