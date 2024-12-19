#! /bin/bash 
#SBATCH -n 8
#SBATCH --mem 32G


module load Cell-Ranger/7.2.0

Scaffolds="./ref/GCF_018143015.1_Lvar_3.0_genomic.fna"
Annotations="./ref/GCF_018143015.1_Lvar_3.0_genomic.gtf"

cellranger mkref --genome=Lvar3 --fasta=$Scaffolds --genes=$Annotations \
	--nthreads=8 \
	--memgb=32 \
	--localcores 8 \
	--localmem 32 
