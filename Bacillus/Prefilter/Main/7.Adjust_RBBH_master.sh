#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=4:mem=20gb

evalue=1e-50

/home/ade110/Scripts/Bacillus/Prefilter/Main/7.Adjust_RBBH_slave.tcl $evalue
