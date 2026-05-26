#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=woodruff.c@wehi.edu.au
#SBATCH -J callIDQuant
#SBATCH --nodes=1
#SBATCH --ntasks=56
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=56
#SBATCH --time 48:00:00
#SBATCH --mem=256G
module unload R
module load R/4.5.3
basepath=$1
cd $basepath
Rscript --vanilla $basepath/blastn_call_ident_quant.R $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
# Usage:          sbatch --mem=100GB --time=4:00:00 blastnCallIDQuant.sh  
#                       basepath whichMock whichSubunit whichCase whichSubMock whichPair uppErrRate maxopnum nmock numcores
#  e.g.  sbatch --mem=240GB --time=48:00:00 blastnCallIDQuant.sh /vast/projects/rrn/ASVcode mockKB 23S 11 0 0 0.01 291 59 24
