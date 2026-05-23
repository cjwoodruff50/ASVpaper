#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=woodruff.c@wehi.edu.au
#SBATCH -J mockStore
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=56
#SBATCH --time 48:00:00
#SBATCH --mem=128G
module unload R
module load R/4.5.1
cd $1
Rscript --vanilla /vast/projects/rrn/ASVtest/make_mock_operon_store.R $1 $2 $3 $4 $5
# Usage:          sbatch --mem=48GB --time=4:00:00 shellScripts/mock_store.sh <basepath>  <maxopnum> <whichSubunit> <maxTotFrags> <numcores>
