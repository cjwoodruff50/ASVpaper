#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=woodruff.c@wehi.edu.au
#SBATCH -J pipeASV4
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=56
#SBATCH --time 48:00:00
#SBATCH --mem=256G
module unload R
module load R/4.5.2
cd $1
# Create the operon store for simulated reads mock microbiome - e.g. mockKB_rrn_C11
### Rscript --vanilla /vast/projects/rrn/ASVtest/make_mock_operon_store.R $1 $2 $3 $4 ${11} ${12} $8 $9 ${13} ${10}
#
# Create the fastq file that is the read library for a simulate reads mock microbiome.  Then
# denoise that microbiome.
### Rscript --vanilla make_mock_and_denoise.R $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
#
# Generate kmer spectra for the ASVs and the reads presented for denoising.
### sbatch --mem=240GB --time=00:20:00 ASVkmerSpectra.sh $1 $2 $3 $4 $5 $6 $7 ${10}
#
# Run blastn on the ASVs, then analyse the alignment data to identify and quantify the
#  strains, species and genera present in the mock microbiome.
sbatch --time=48:00:00 blastnCallIDQuant.sh $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
# OR (if not a large computation - e.g. Sereika or Srinivas datasets)
### Rscript --vanilla blastn_call_ident_quant_NWalign_v02.R $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
#
# Pipeline Script Usage: 
#  Input parameters in order are:-
#               [1:5] basepath, whichMock, whichSubunit, whichCase, whichSubMock
#               [6:10]  whichPair, upThreshErrRate, totOperons, totalStrains, numcores,
#               [11:13]   meanQscore, sdQscore, maxReadCount
#
# IMPORTANT: Each script must be finished before starting the next. 
#            Control execution by editing the bash file, placing or removing  ###  before the appropriate calls.
#
#sbatch --mem=48GB --time=4:00:00 shellScripts/ASVpipeline.sh <basepath>  <whichMock> <whichSubunit> <whichCase> <whichSubMock> 
#                                    <whichPair> <uppThreshErrRate> <totOperons> <totalStrains> <numcores>
#                                      <meanQscore> <sdQscore> <maxReadCount>
