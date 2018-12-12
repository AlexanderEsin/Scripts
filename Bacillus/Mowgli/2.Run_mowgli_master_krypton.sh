#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=50gb
#PBS -J 1-7298

export TMPDIR=/mnt/storage/home/aesin/temp
mkdir -p $TMPDIR

module load Anaconda3

tclsh8.6 /mnt/storage/home/aesin/Scripts/Bacillus/Mowgli/2.Run_mowgli_slave_krypton.tcl $PENALTY $PBS_ARRAY_INDEX >> /mnt/storage/home/aesin/Work/Bacillus/Mowgli/trackMow$PENALTY\.txt