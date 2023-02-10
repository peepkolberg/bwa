## bwa
*Nextflow workflow for aligning fastq reads to a reference genome*

To run:

In case you already have the bwa index, skip to step 5.

1. Copy your reference genome and its index (`.fai`) to directory `bwa_idx`.
2. In file `run_bwa_index.sh` line 13 change *ref_gen_prefix*'s value to your reference genome's filename but strip "`.fa`" from the end.
3. Make sure your current working directory is `bwa/`.
4. Submit `run_bwa_index.sh` to HPC. This will create the bwa index in the `bwa_idx` folder.

- - - - - - -

5. In file `run_bwa.sh`, line 15, change *ref_gen_prefix*'s value to the same prefix you used to create the bwa index (step 2).
6. In file `run_bwa.sh`, make sure --ref_gen points to the bwa index **directory** and --samples points to a `.tsv` file with your samples (example in `data/samples.tsv`).
7. Submit `run_bwa.sh` to HPC.

Results will appear in `results/`
