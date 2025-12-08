# ASVpaper Code

## Introduction

Here is the code underlying a draft paper on identification and quantification of bacterial microbiomes using ribosomal RNA amplicons processed by an Oxford Nanopore Techology's (ONT) sequencer. That study had real sequenced data from two published papers, and simulated nanopore-sequenced reads for a bacterial microbiome based on a published indicative human gut microbiome of healthy people.

## Outline of Key computer code and processing pipeline.

Steps 1 to 4 below are only relevant for the simulated reads mock microbiomes. They generate a store of multifasta files for each operon in the mock microbiome and construct a single multifastq file that corresponds to the simulated library from a sequencing run on the microbiome.  Each multifasta file has multiple fasta records that are indpendent noisy simulated reads of a single operon. The number of reads is determined to be at least as large as the number that will be required for the intended relative abundances of the generated microbiome.

The 3 items below are executed through make_mock_store.R.

  1. Select the operons of the strains that are in the mock microbiome and generate individual fasta files for each one.
  2. Determine the number of simulated reads to be generated for each operon, based on an approximate number of reads that the library of the pre-filtered mock microbiome is to have for the operons of the most abunmdant strain.
  3. For each operon submit a batch job that generates a moderate excess of simulated reads over the number actually required for the mock microbiome library.  Note that the reads will be randomly selected from the store without replacement.

Item 4 is executed throgh make_mock_denoise.R, Part 1.

  4. On completion of simulated read generation of all operons, randomly select the required number of reads for each operon and store as a single, multifastq file.  This is the dataset of the designed mock microbiome.

The remaining processing pipeline applies to any single fastq file that corresponds to a ribosomal RNA 16S , 23S or 16S-ITS-23S amplicon library for some microbiome.  

  5. Submit the mock microbiome fastq library for quality and length filtering and subsequent denoising using teh Robust Amplicon Denoising algorithm of Kumar et al.. These steps are exewcuted in Julia code.  The Julia  script used here is a simple elaboration of of that provided by Ben Murrell (cite(Murrell_github, and private communication).  Item 5 can be executed through part 2 of  make_mock_denoise.R .
