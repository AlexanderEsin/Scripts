#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=10gb

module load anaconda3/personal

evalue=1e-10
echo $evalue

tclsh8.6 /home/ade110/Scripts/Staph/RBBH_MCL/4.RBBH_collate_slave.tcl $evalue
