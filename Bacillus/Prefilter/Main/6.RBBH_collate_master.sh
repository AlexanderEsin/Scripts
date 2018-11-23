#!/bin/sh
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=4:mem=50gb

evalue=1e-50

/home/ade110/Scripts/Bacillus/Prefilter/Main/6.RBBH_collate_slave.tcl $evalue