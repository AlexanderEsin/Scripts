#!/bin/sh
#PBS -l walltime=25:00:00
#PBS -l mem=8gb
#PBS -l ncpus=8

module load raxml/7.2.8

direct=$SCRATCH/RAxML/Geobac
cd $direct


echo 'set direct /scratch/ade110/RAxML/Geobac

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct
file mkdir Tree_Build
file mkdir Best_bootstrap_tree

cd $direct/Stitched_core_genome

set fams [glob *.phylip]

foreach fam $fams {
	regsub ".phylip" $fam {} fname
	exec raxml -f a -s $fam -n $fname.txt -m PROTCATDAYHOFF -p 1234 -x 10001 -N 1000 -T 8
}

set reduced [glob -nocomplain *.reduced]
foreach red $reduced {
	file delete $red
}

set tree_files [glob RAxML*]
foreach file $tree_files {
	file rename $file $direct/Tree_Build/$file
}

set key [glob *KEY*]
openfile $key
set key $data

cd $direct/Tree_Build

set pata {k+?[0-9]+?:+?}
set patb {\t+?.+?\n+?}
set patd {.+?\n+?}

set besttrees [glob *bipartitions.*]
foreach tree $besttrees {
	openfile $tree
	set branchids [regexp -all -inline $pata $data]
	foreach id $branchids {
		regsub {:} $id {} idtrunc
		regexp $idtrunc$patb $key hit
		regexp -nocase {\t+?[A-Z]+?.+?\n+?} $hit match
		set match [string trim $match]
		regsub $idtrunc $data $match data
	}
	set out [open Geobac_bootstrap_tree.txt w]
	puts $out $data
	close $out
}

file copy Geobac_bootstrap_tree.txt $direct/Best_bootstrap_tree' > RAxML_Geobac.tcl

chmod +x RAxML_Geobac.tcl

tclsh8.5 $direct/RAxML_Geobac.tcl

wait
