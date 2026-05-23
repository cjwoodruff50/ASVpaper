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
cd $1
module load julia/1.10.4
module load python/3.13.0
export JULIA_DEPOT_PATH=/vast/projects/rrn/ASVcode/depot
julia juliaENVwipe.jl
julia RADrun.jl $1 $2 $3 $4 $5 $6 $7
#julia RADrun.jl $1 --whichMock=$2 --whichSubunit=$3 --whichCase=$4 $5 $6 $7
# Usage:          sbatch --mem=124GB --time=2:00:00 juliaRADrunSA.sh /vast/projects/rrn/ASVtest mSriSA2 rrn 03 0 2 0.0125
#                 Unnamed Params are  basepath whichSubMock, whichPair, uppThreshErrRate
#  N.B.  The last parameter is an integer whose reciprocal is the upper threshold for read error rate.
