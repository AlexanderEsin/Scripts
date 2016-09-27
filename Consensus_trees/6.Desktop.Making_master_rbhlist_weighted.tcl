#!/usr/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
set direct /users/aesin/desktop/Consensus_trees

set eval -10
set rbbh_dir Rbbh_all_$eval

set taxon Geobac/Outgroups/With_excl_geos
set org $taxon/$rbbh_dir


###########################################################################
cd $direct/$org
if {[file exists Master_rbh_weight.txt] == 1} {
	file delete Master_rbh_weight.txt
}

cd $direct/$org/Bin

set rbbh_list [glob *.txt]
set i 1

foreach file $rbbh_list {
	openfile $file
	set rhead [string first \n $data]
	set rbbhs [string range $data $rhead end]
	append rbbhs "\n"
	set rbb1 [string trim $rbbhs]
	set out [open $direct/$org/Master_rbh_weight.txt a]
	puts $out $rbb1
	close $out
	puts "$i/[llength $rbbh_list]"
	incr i
}

puts "\n\n\nMASTER WEIGHTED LIST DONE\n\n\n"

# if {[file exists $direct/$taxon/Out_split] == 1} {
# 	cd $direct/$taxon/Out_split
# 	set i 1
# 	set dirs [glob -type d *]
# 	foreach dir $dirs {
# 		cd $direct/$taxon/Out_split/$dir
# 		set files [glob -nocomplain *.tsv]
# 		foreach file $files {
# 			file rename -$file $direct/$taxon/Out_all
# 			puts "$i"
# 			incr i
# 		}
# 	}

# 	cd $direct/$taxon/Out_all
# 	if {[llength [glob *tsv]] == [expr $i -1]} {
# 		file delete -force $direct/$taxon/Out_split
# 	}

# 	puts "\n\n\nALL OUT FILES TRANSFERED BACK TO OUT_ALL\n\n\n"
# }

# set i 1
# cd $direct/$taxon/Out_nw
# set files [glob -nocomplain *.tsv]
# foreach file $files {
# 	file rename $file $direct/$taxon/Out_all
# 	puts "$i"
# 	incr i
# }

