#!/bin/sh -f
# qsub pacbio_SV.pbs

#PBS -l nodes=1:ppn=20
##PBS -l pmem=1gb
#PBS -A open
#PBS -l walltime=23:50:00
#PBS -j oe
#PBS -o pacbio.out
#PBS -N pacbio
#PBS -M juc326@psu.edu
#PBS -m abe
#######################################
#ufo.fastq contains the PacBio long-read sequencing data

module load gcc/5.3.1 samtools/0.1.19 
export PATH=~/work/ngmlr-0.2.7:$PATH
export PATH=~/work/sniffles:$PATH
cd $PBS_O_WORKDIR


ngmlr -t 20 -r ~/work/Zea_mays.AGPv4.fa -q ufo.fastq -o ufo.bam
samtools sort ufo.bam ufo.sort.bam
samtools index ufo.sort.bam

sniffles -m ufo.sort.bam -v ufo.SV.vcf

#######################################
echo "Job Started"
echo "Job Ended"
