#!/bin/bash
#SBATCH --ntasks=32
#SBATCH --mem=32G
#SBATCH --partition=RT
#SBATCH --job-name="TEcandidates"
#SBATCH --comment="TEcadidates_urchin"
#SBATCH --mail-user="panyushev@nextmail.ru"
#SBATCH --mail-type=ALL

module load java/13 
module load python/3.7

TEcandidates.sh -t=32 -r=32 -c=1 -l=20 -te=Spur_repeats.gff3 -g=Spur_5.0_genome.fasta -fq=. -m=PE -N=10000 1> TEcandidates_urchin.log 2> TEcandidates_urchin.err

