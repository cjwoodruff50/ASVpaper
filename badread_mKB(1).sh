#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=woodruff.c@wehi.edu.au
#SBATCH -J mockMake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=240GB
cd $1
total_frags=$7
whichSubunit=$2
OpLow=$3
OpHi=$4
identMode=$5
identSD=$6
echo $whichSubunit $OpLow $OpHi $identMode $identSD $total_frags
source $1/text/frag"$whichSubunit"length.sh
seqLen=("${lengths[@]}")
processes=$((OpHi-OpLow+1))
PATH=$PATH:/home/users/allstaff/woodruff.c/.local/bin
echo $PATH
for p in $(seq $OpLow $OpHi) 
do
    total_bases=$((total_frags * seqLen[$p-1]))	
    echo $total_bases
    bases=$((total_frags * seqLen[$p-1] / processes))
    log=text/badread_"$whichSubunit"_"$p".log
    reads=fastq/store/"$whichSubunit"/reads_"$whichSubunit"_Op_"$p".fastq
    badread simulate --reference fasta/mockKB_"$whichSubunit"_Op_"$p".fasta --quantity $((total_bases)) \
        --glitches 100000,1,1 --junk_reads 0 --random_reads 0 --chimeras 0 \
         --error_model nanopore2023 --qscore_model nanopore2023 \
	 --length $((seqLen[$p-1])),5 --identity $((identMode)),$((identSD)) \
           --start_adapter_seq "" --end_adapter_seq "" 2> $log 1> $reads
done

# sbatch --mem=240GB  --time=48:00:00  badread_mKB.sh <basepath>  <whichSubunit> <OpLow> <OpHi> <identMode> <identSD> <totfrags>


