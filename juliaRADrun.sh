#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=woodruff.c@wehi.edu.au
#SBATCH -J JRADrun
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=1
#SBATCH --time 8:00:00
#SBATCH --mem=128G
cd /vast/projects/rrn/ASVcode
module load julia/1.10.4
module load python/3.13.0
export JULIA_DEPOT_PATH=/vast/projects/rrn/ASVcode/depot
julia juliaENVwipe.jl
julia RADrun.jl --whichMock=$1 --whichSubunit=$2 --whichCase=$3 $4 $5 $6
# Usage:          sbatch --mem=124GB --time=2:00:00 juliaRADrunSA.sh whichMock=mSriSA2 whichSubunit=23S whichCase=03 0 2 100
#                 Unnames Params are  whichSubMock, whichPair, 1/UppThreshErrRate
#  N.B.  The last parameter is an integer whose reciprocal is the upper threshold for read error rate.
