#!/usr/local/bin/tclsh

## This script takes a set of gene tree files with tips formatted as "BINOMIAL_NAME {CLASS} {PROTEIN_ID}" and converts the tip labels to "BINOMIAL_NAME". This will results in duplicate taxon names for those gene trees where there is paralogy. The tree may contain bootstrap values ##

###########################################################################
## Procs ##
source ~/Dropbox/Scripts/General_utils.tcl

###########################################################################
## Parent directory ##
set direct ~/desktop/Geo_final_trees

## Input and output directories ##
set input_direct $direct/Final_BS
set output_direct $direct/Relabelled

## Set the output file prefix (if any) ##
set prefix "Relabelled_"

## Make output directory ##
file mkdir $output_direct

## Go to input directory and get a list of trees based on the GLOB pattern ##
cd $input_direct
set trees [glob *_tree.txt*]

## Set the pattern to remove from the tip labels ##
set pata {(?:\s\{)+?.+?\}+?}

set i 1
foreach tree $trees {

	## Get the tree name (in this case the gene family number) and read in the file ##
	set tree_num [string range $tree 0 [string first \_ $tree]-1]
	set tree_data [string trim [openfile $tree]]

	## Substitute all the occurence of the pattern in the tree file ##
	regsub -all $pata $tree_data {} new_tree_data

	## Write out the files into the output directory with a prefix ##
	set out [open $output_direct/$prefix$tree_num\.nwk w]
	puts $out [string trim $new_tree_data]
	close $out

	## Show progress and update counter ##
	puts "$i / [llength $trees]"
	incr i
}

puts "Relabelling completed"