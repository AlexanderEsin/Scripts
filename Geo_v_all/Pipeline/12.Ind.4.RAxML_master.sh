#!/bin/sh
#PBS -l walltime=500:00:00
#PBS -l select=1:ncpus=10:mem=16gb

module load raxml/2015-02-17 

org=Geo_v_all

chmod +x /home/ade110/Scripts/$org/Pipeline/12.Ind.4.RAxML_HPC_slave.tcl

tclsh8.5 /home/ade110/Scripts/$org/Pipeline/12.Ind.4.RAxML_HPC_slave.tcl 1