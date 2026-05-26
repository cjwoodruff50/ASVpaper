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
module load R/4.5.3
cd $1
Rscript --vanilla ASV_reads_kmerSpectra.R $1 $2 $3 $4 $5 $6 $7 $8
# Usage:          sbatch --mem=240GB --time=48:00:00 shellScripts/ASVkmerSpectra.sh  basepath whichMock whichSubunit 
#                            whichCase whichSubMock whichPair uppErrRate 56
#  e.g.  sbatch --mem=240GB --time=48:00:00 shellScripts/ASVkmerSpectra.sh  /vast/projects/rrn/ASVcode mSriSA2 rrn 03 0 2 0.009 56
