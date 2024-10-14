#!/bin/bash

#SBATCH --job-name=bwa
#SBATCH --partition=amd
#SBATCH --time=5-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

# TODO: Tweak bwa_mem memory parameter. 8 GB was too little, 32 GB worked, and I couldn't be bothered with testing further at the time.

module load bwa/0.7.17
module load any/samtools
module load nextflow

ref_gen_prefix=Homo_sapiens.GRCh38.dna.primary_assembly  # Change this to your reference genome's filename but strip '.fa' from the end

nextflow -log logs/.nextflow.log run bwa.nf -profile tartu_hpc \
    --samples data/samples.tsv \
    --ref_gen bwa_idx \
    --ref_gen_prefix $ref_gen_prefix \
    -resume
