#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=woodruff.c@wehi.edu.au
#SBATCH -J ASVrdskmS
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=56
#SBATCH --time 48:00:00
#SBATCH --mem=120G
module unload R
module load R/4.5.1
cd /vast/projects/rrn/ASVcode
Rscript --vanilla /vast/projects/rrn/RscriptsArchive/current/ASV_reads_kmerSpectra.R $1 $2 $3 $4 $5 $6 $7
# Usage:          sbatch --mem=240GB --time=48:00:00 shellScripts/ASVkmerSpectra.sh  whichMock whichSubunit 
#                            whichCase whichSubMock whichPair 56
#  e.g.  sbatch --mem=240GB --time=48:00:00 shellScripts/ASVkmerSpectra.sh  mSriSA2 rrn 03 0 2 0.009 56