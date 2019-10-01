#!/bin/sh -f
# qsub STAR_2pass.pbs

#PBS -l nodes=1:ppn=20
##PBS -l pmem=1gb
#PBS -A open
#PBS -l walltime=23:50:00
#PBS -j oe
#PBS -o STAR_2pass.out
#PBS -N STAR_2pass
#PBS -M juc326@psu.edu
#PBS -m abe
#######################################

module load gcc/5.3.1 samtools/0.1.19 bedtools/2.26.0
export PATH=~/work/STAR-2.5/bin/Linux_x86_64:$PATH
export PATH=~/work/subread-1.5.3-Linux-x86_64/bin:$PATH
cd $PBS_O_WORKDIR

# STAR make genome index
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ~/work/STAR.idx --genomeFastaFiles ~/work/Zea_mays.AGPv4.fa --sjdbGTFfile ~/work/Zea_mays.AGPv4.32.gtf --sjdbOverhang 49


for DIR in */ ; do
    cd "$DIR"
    # first alignment
    STAR --runThreadN 20 --genomeDir ~/work/STAR.idx --readFilesIn data.fa  --outSAMtype BAM SortedByCoordinate
    mkdir idx
    # generate new genome index using splice junction information
    STAR --runThreadN 20 --runMode genomeGenerate --genomeDir idx/STAR-2pass.idx --genomeFastaFiles ~/work/Zea_mays.AGPv4.fa --sjdbFileChrStarEnd SJ.out.tab --sideOverhang 49
    # 2nd alignment
    STAR --runThreadN 20 --genomeDir idx/STAR-2pass.idx --readFilesIn data.fa --outSAMtype BAM SortedByCoordinate
    cd ../
done

featureCounts -T 20 -t exon -g gene_id -a ../Zea_mays.AGPv4.32.gtf -o Counts_exon_geneESW.txt */Aligned.sortedByCoord.out.bam


#######################################
echo "Job Started"
echo "Job Ended"