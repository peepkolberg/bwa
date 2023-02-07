#!/bin/bash

#SBATCH --job-name=bwa_idx
#SBATCH --partition=amd
#SBATCH --time=1-00:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
# If runs out of memory, add more.

module load bwa/0.7.17

ref_gen_prefix=Homo_sapiens.GRCh38.dna.primary_assembly  # Change this to your reference genome's filename but strip '.fa' from the end

cd bwa_idx
bwa index -a bwtsw -p ${ref_gen_prefix} ${ref_gen_prefix}.fa
