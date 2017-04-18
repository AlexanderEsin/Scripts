#!/bin/sh
#PBS -l walltime=40:00:00
#PBS -l select=1:ncpus=8:mem=12gb

##########################
direct=/home/ade110/Scripts
org=Geo_v_all
script=X.Indiv_align_slave.tcl
tclsh=1
bash=0
##########################

module load clustal-omega/1.2 

chmod +x $direct/$org/$script
# Make sure there is a space between each bit #
if [ "$tclsh" = 1 ]; then
	tclsh8.5 $direct/$org/$script
elif [ "$bash" = 1] ; then
	$direct/$org/$script
fi

