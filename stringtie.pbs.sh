#!/bin/sh -f
# qsub stringtie.pbs

#PBS -l nodes=1:ppn=20
##PBS -l pmem=1gb
#PBS -A open
#PBS -l walltime=23:50:00
#PBS -j oe
#PBS -o stringtie.out
#PBS -N stringtie
#PBS -M juc326@psu.edu
#PBS -m abe
#######################################
# StringTie v1.3.4 was installed with conda
module load gcc/5.3.1 
cd $PBS_O_WORKDIR

for DIR in */ ; do
    echo "$DIR".gtf >> mergelist.txt
    cd "$DIR"
    stringtie -p 20 -G ~/work/Zea_mays.AGPv4.32.gtf -o "$DIR".gtf Aligned.sortedByCoord.out.bam

stringtie --merge -p 20 -G ~/work/Zea_mays.AGPv4.32.gtf -o stringtie_merged.gtf mergelist.txt


#######################################
echo "Job Started"
echo "Job Ended"
