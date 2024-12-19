#!/bin/bash

#3_count.sh

#cjm 2024/12/05

#job submission script for making count tables from reads for scRNAseq analysis. 

#define some variables to make this more reuseable

#transcriptome index path
txome="/work/cjm124/scRNAanalysis/Lvar3"

#data directory path
fastqs="/work/cjm124/scRNAanalysis/reads"

cd $fastqs

for samp in `cat sample.list`; do
echo "submitting job for $samp via wrap"
echo "command sent: module load Cell-Ranger/7.2.0 && cellranger count --id=$samp --transcriptome=$txome --fastqs=$fastqs --expect-cells=3000 --sample=$samp  --localcores=6 --localmem=15"
sbatch -N 1 -n 6 --mem=15G -p scavenger -J fq2count.$reads --mail-user=c.m@duke.edu --mail-type=FAIL --wrap="module load Cell-Ranger/7.2.0 && cellranger count --id=$samp --transcriptome=$txome --fastqs=$fastqs --expect-cells=3000 --sample=$samp --localcores=6 --localmem=15"
done
