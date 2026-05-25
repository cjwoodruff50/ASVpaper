# ASVpaper Code

## Introduction

Here is the code underlying a draft paper on identification and quantification of bacterial microbiomes using ribosomal RNA amplicons processed by an Oxford Nanopore Techologies' (ONT) sequencer. That study had real sequenced data from two published papers, and simulated nanopore-sequenced reads for a bacterial microbiome based on a published indicative human gut microbiome of healthy people.


## Outline of Key computer code and processing pipeline.

Steps 1 to 4 below are only relevant for the simulated reads mock microbiomes. They generate a store of multifasta files for each operon in the mock microbiome and construct a single multifastq file that corresponds to the simulated library from a sequencing run on the microbiome.  Each multifasta file has multiple fasta records that are independent noisy simulated reads of a single operon. The number of reads is determined to be at least as large as the number that will be required for the intended relative abundances of the generated microbiome.

The 3 items below are executed through make_mock_store.R.

  1. Select the operons of the strains that are in the mock microbiome and generate individual fasta files for each one.
  2. Determine the number of simulated reads to be generated for each operon, based on an approximate number of reads that the library of the pre-filtered mock microbiome is to have for the operons of the most abunmdant strain.
  3. For each operon submit a batch job that generates a moderate excess of simulated reads over the number actually required for the mock microbiome library.  Note that the reads will be randomly selected from the store without replacement.

Item 4 is executed through make_mock_denoise.R, Part 1.

  4. On completion of simulated read generation of all operons, randomly select the required number of reads for each operon and store as a single, multifastq file.  This is the dataset of the designed mock microbiome.

The remaining processing pipeline applies to any single fastq file that corresponds to a ribosomal RNA 16S , 23S or 16S-ITS-23S amplicon library for some microbiome.  

  5. Submit the mock microbiome fastq library for quality and length filtering and subsequent denoising using the Robust Amplicon Denoising algorithm of Kumar et al.. These steps are executed in Julia code.  The Julia  script used here is a simple elaboration of that provided by Ben Murrell (cite(Murrell_github, and private communication).  Item 5 can be executed through part 2 of  make_mock_denoise.R, in which a shell script, juliaRADrun.sh, is executed which sets up for calling julia and a julia script, RADrun.jl, implements the filtering, generation of 2 plots, and calling of RAD (the Robust Amplicon Denoiser). Tailoring of the call is done by having the shell script called with a set of parameters - described elsewhere in this document - which are then passed to the julia script.

Having denoised the mock microbiome library, and hence generated a set of amplicon sequence variants (ASVs), the next part of the processing pipeline uses these ASVs to identify what strains, species or genera are present, and their relative abundances. That analysis, elaborated in items 6 to 12 below, is all included in the code  blastn_call_ident_quant.R .  

  6. Set up characterising data on the mock microbiome, such as strain names, number of operons for each strain, expected relative abundances.  Note that, because the microbiomes considered in this analysis are mock microbiomes, there are expected relative abundances available.  These may be determined based on amounts of DNA for each strain that was prepared in wet lab operations for sequencing, or as designed if a mock microbiome is generated from simulated reads.
  7. Read the key text files generated as output by RAD, specifying the sequence of each ASV and the reads associated with each ASV.
  8. Set up a call of blastn to align each ASV against the reference strain database.  We call blastn via R's system2() function.  Two such calls are made to generate different output sets.
  9. Parse the blastn output to determine the best match of each ASV sequence to the operons in the database, breaking ties for best alignment by selecting the first listed database entry.  Each ASV then has a strain label.
  10. Save the read counts for each labelled ASV.  Using the number of operons for the relevant strain convert operon counts to cellular counts for each strain, and then return the relative cellular abundance of each strain identified, distinguishing strains that were in the designed mock microbiome and those not in the mock microbiome.
  11. Generate 2D UMAP projection plots that illustrate the identification performance being achieved with the analysis.
  12. Generate scatterplots of observed proportions vs. designed (or expected) proportions - or the logarithm of these quantitites - showing the quantification performance of the analysis.

## Running the Code
The R code, shell script, and julia scripts are all assumed to be in the directory specified by **basepath**.  This directory has an associated sub-directory structure shown at  
![dirStruct](https://github.com/user-attachments/files/23666822/ASVcodeLayout.pdf). In that figure the directory  ASVcode  corresponds to **basepath**.

The code has been written to run any of 21 datasets, and could be modified relatively straightforwardly to process other datasets, whether from real sequenced microbiomes or simulated reads of microbiomes.  In the following we detail how to run one of the simulated read datasets and one of the real sequenced datasets.  To run the code it is necessary to load a number of files from Figshare, including the two blastn databases that are used (here identified as blastdb1 and blastdb2), and the 16S, 23S and rrn sequences of all operons of each strain in the strain-only database, blastndb1. Details of the Figshare material and how to use this material to set up the database environment is given below.

### Figshare repository description and content use.

#### Description
Figshare holds the set of files related to the blast reference databases 
and the 18 datasets of real nanopore-sequenced data.  These are accessed at 
  https://figshare.com/articles/dataset/Sereika_10_datasets_Srinivas_8_datasets/31052200 
and  
  https://figshare.com/articles/dataset/DB1_DB2_multifasta_specR_KGDB_files/31041886    
The blast-related set of files are in the item  DB1_DB2_multifasta_specR_KGDB_files.  The datasets are in two parts - that from Srinivas M et al. (2025) and that from extracting 16S and 23S rRNA gene sequences from whole genome sequencing of microbiomes in the study of Sereika M et al. (2022).  The 
Srinivas datasets are individually available as gzipped files, while the Sereika datasets, 
being much smaller, were tarred before gzipping.

#### How to Use the Content
It is assumed that the directory structure detailed in the Github repository for this 
project has been established. 

For the   DB1_DB2_multifasta_specR_KGDB_files  entry follow the steps below:-
   1. Download   DB1_DB2_multifasta_specR_KGDB_files to ...basepath/GROND
   2. Unzip and untar .  This gives a number of zipped and (mostly) tarred files.
   3. Move the blastdb1 file to the sub-directory ...basepath/GROND/blastdb1 and untar, so creating the 
      strain-only database.
   4. Analogously for the blastdb2 file to give the species-level database.
   5. Move the multifasta16S tarred file to the sub-directory ...basepath/GROND/multifasta16S and untar.
      Likewise for the .../multifasta23s and .../multifastarrn zipped files.  
   6. unzip and untar the remaining downloads and leave them in the 
      ...basepath/GROND sub-directory    
That completes setting up of the database structure and links to the King et al. (2019) Gut 
Feeling Knowledge Base.  This allows creation of the simulated reads datasets, and also carrying  
out of blastn alignments against the Walsh et al.'s (2024) GROND-derived operon databases.  

The dataset files should be placed in the ...basepath/fastq sub-directory and then unzipped 
(and untarred for the Sereika datasets).   

#### References: (full reference details in the primary document).
  * Srinivas M et al. Scientific Reports 2025
  * Sereika M et al. Nature Methods 2022
  * King C.H. et al. PLoS One 2019
  * Walsh C.J. et al. Microbial Genomics 2024

 

### Running Simulated reads dataset mockKB_rrn_C11
 1. Ensure that the following code scripts are in the **basepath** directory:-
    * badread_mKB.sh
    * RADrun.jl and juliaEnvwipe.jl
    * make_mock_operon_store.R
    * make_mock_and_denoise.R
    * blastn_call_ident_quant.R
    * juliaRADrun.sh
      
2. Determine the design of the mock microbiome.
    * The free parameters for this are the number of reads of the most abundant strain, and the region of the rrn operon (16S, 23S, rrn)
  

### Instructions for quickly running datasets 
To directly run the processing pipeline for one of the datasets proceed as follows:-

0. Set up the directory structure detailed below for **basepath**  (see section Running the Code above).
   In the example calls below  **basepath*  is /vast/projects/rrn/ASVtest
     * CaseStr            is 11 for all simulated reads examples, C03 for real rea datasets
     * meanQ              is 30 (only relevant for simulated reads datasets)
     * stdQ               is 4  (only relevant for simulated reads datasets)
     * Numops             is 291 (only relevant for simulated reads datasets)
     * Nstrain            is 59 (only relevant for simulated reads datasets)
     * ErrRateThresh      is 0.01 (or 0.015)
     * whichSubunit       is 16S or 23S for Ser datasets, rrn for Sri datasets, and 16S, 23S or rrn for mockKB datasets
     * whichMock          is only non-zero for mSerS (Sereika) datasets
     * whichPair          is only non-zero for mSri (Srinivas) datasets
     * Numops Nstrain are only non-zero for simulated read (mockKB) datasets
     * maxNop             is the maximum number of simulated reads for any operon in a mock microbiome. It is a key determinant of the library size(only relevant for simulated reads datasets).
   
   
2. Create an operon store for simulated reads mock microbiome - e.g. mockKB_23S_C11

   Rscript --vanilla make_mock_operon_store.R basepath whichMock whichSubunit whichCase
                meanQ stdQ Numops Nstrain maxNop nthreads
   
   Sample call (noting that mockKB has 59 strains and 291 operons):-
   
     Rscript --vanilla make_mock_operon_store.R /vast/projects/rrn/ASVtest 
               mockKB 23S 11 30 4 291 59 30000 12           
            
4. Create the fastq file that is the read library for a simulated reads mock microbiome.
   Denoise the microbiome created or one that is from real sequencing.
   
   Rscript --vanilla make_mock_and_denoise.R basepath whichMock whichSubunit whichCase
     whichSubMock whichPair uppThreshErrRate Numops Nstrain numcores
   
   Sample calls  
    e.g. for mockKB_23S_C11 (whichSubMock=0, whichPair=0 for any simulated reads dataset)
   
    Rscript --vanilla make_mock_and_denoise.R /vast/projects/rrn/ASVtest mockKB 23S 11 0 0
                        0.01 291 59 2
   
   
    e.g. for mSriSA2_rrn_C03 (for Srinivas data whichMock=0, whichPair>0)
   
    Rscript --vanilla make_mock_and_denoise.R /vast/projects/rrn/ASVtest mSriSA2 rrn 03 0 
                        2 0.015 0 0 2
   
   
    e.g. for mSerS3_16S_C03 (whichSubMock=2; for Sereika subsampled data whichMock is non-zero,  whichPair=0)
   
    Rscript --vanilla make_mock_and_denoise.R /vast/projects/rrn/ASVtest mSerS3 16S 03 3 0
                 0.01 0 0 2
                                    
   
5. Generate kmer spectra for the ASVs and the reads presented for denoising.

   sbatch --mem=240GB --time=00:20:00 ASVkmerSpectra.sh basepath whichMock whichSubunit
      whichCase whichSubMock whichPair uppThreshErrRate nthreads
   
   Sample call:-    
    e.g. sbatch --mem=240GB --time=00:10:00 ASVkmerSpectra.sh /vast/projects/rrn/ASVtest
            mSriSA2 rrn 03 0 2 0.015 56   


7. Run blastn on the ASVs, then analyse the alignment data to identify and quantify the
   strains, species and genera present in the mock microbiome.
   
   Rscript --vanilla blastn_call_ident_quant.R <basepath> whichMock whichSubunit whichCase
        whichSubMock whichPair uppThreshErrRate Numops Nstrain numcores
   
    e.g. for mSriSA2_rrn_03 with error rate threshold 0.0125
   
     Rscript --vanilla blastn_call_ident_quant.R /vast/projects/rrn/ASVtest 
                   mSriSA2 rrn 03 0 2 0.0125 0 0 2
                   
   Running this call for mockKB is likely to lead to exceeding the allowed resources on a
   terminal session.  Hence a slurm job is needed, via the shell script  blastnCallIDQuant.sh
   which has the same 10 parameters.
   
   sbatch --mem=100GB --time=4:00:00 blastnCallIDQuant.sh  basepath  whichMock 
            whichSubunit whichCase whichSubMock whichPair uppErrRate maxopnum nmock 56
   
   Sample call:-    
   e.g. sbatch --mem=240GB --time=48:00:00 blastnCallIDQuant.sh /vast/projects/rrn/ASVtest 
              mockKB 23S 11 0 0 0.01 291 59 16                
                           
            
