#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=50gb

module load anaconda3/personal

evalue=1e-10
echo $evalue

/home/ade110/Scripts/Bacillus/RBBH_MCL/4.RBBH_collate_slave.tcl $evalue
