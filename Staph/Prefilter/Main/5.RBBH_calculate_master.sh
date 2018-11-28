#!/bin/sh
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-64

evalue=1e-50

/home/ade110/Scripts/Staph/Prefilter/Main/5.RBBH_calculate_slave.tcl $evalue $PBS_ARRAY_INDEX