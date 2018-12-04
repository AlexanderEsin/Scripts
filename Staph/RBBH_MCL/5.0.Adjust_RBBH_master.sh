#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=10gb

/home/ade110/Scripts/Staph/RBBH_MCL/5.0.Adjust_RBBH_slave.tcl $evalue
