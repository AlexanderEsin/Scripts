set direct ~/desktop/Consensus_trees/Rooting
set ival 5
set eval -10
set org Bacillaceae/$eval

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

# ###########################################################################
# cd $direct/$org
# file mkdir $direct/$org/Tree_build_$ival/Indiv_aligns
# file mkdir $direct/$org/Final_trees/Indiv_aligns

# cd $direct/$org/Group_fastas_MSA_$ival

# set aligns [glob *fas]
# foreach align $aligns {
# 	regsub {.fas} $align {} fname
# 	file mkdir $direct/$org/Tree_build_$ival/Indiv_aligns/$fname
# 	file copy -force $align $direct/$org/Tree_build_$ival/Indiv_aligns/$fname

# 	cd $direct/$org/Tree_build_$ival/Indiv_aligns/$fname

# 	set reduced [glob -nocomplain *.reduced]
# 	foreach red $reduced {
# 		file delete $red
# 	}

# 	set tree_files [glob -nocomplain RAxML*]
# 	foreach file $tree_files {
# 		file delete $file
# 	}

# 	catch {exec raxml -f d -D -s $align -n $fname\.txt -m PROTCATAUTO -p 1234 -N 2 -T 6}

# 	set besttrees [glob -nocomplain *bestTree.*]
# 	file copy -force [lindex $besttrees 0] $direct/$org/Final_trees/Indiv_aligns

# 	cd $direct/$org/Group_fastas_MSA_$ival
# }
# cd $direct/$org/Final_trees/Indiv_aligns
# set indiv_trees [glob -nocomplain *txt]
# foreach tree $indiv_trees {
# 	set name [string range $tree 0 end-4]
# 	set new_extension "$name\.tree"
# 	file rename $tree $new_extension
# }

# ###########################################################################

# cd $direct/$org/Final_trees/Indiv_aligns

# if {[file exists MRC_master_tree.tre] == 1} {
# 	file delete MRC_master_tree.tre
# }

# set raxml_files [glob -nocomplain *txt]
# foreach file $raxml_files {
# 	file delete $file
# }

# set ML_trees [glob -type f *]

# set master_tree_file ""
# foreach file $ML_trees {
# 	openfile $file
# 	set tree [string trim $data]
# 	append master_tree_file "$tree\n"
# }
# set out [open MRC_master_tree.tre w]
# puts $out [string trim $master_tree_file]
# close $out

###
openfile $direct/$org/Tree_build_$ival/MRC/Species_KEY.txt
set key $data
set pata {k+[0-9]+}
set patc {k+?[0-9]+?\:+?}
set patb {\t+.+}
###
cd $direct/$org/Tree_build_$ival/MRC
set master_tree_file [glob *bootstrap*]
openfile $master_tree_file

###
regsub -all {(k+[0-9]+)} $data {\1:} data
set branchids [regexp -all -inline $patc $data]
foreach k_val $branchids {
	regsub {:} $k_val {} k_val_trunc
	regexp -line $k_val_trunc$patb $key hit
	set binomial [lindex [split [string trim $hit] \t] 1]
	#set id [lindex [split [string trim $hit] \t] 2]
	# If you use k_val_trunc below, you may substitute k17 when you want to find k1 - so searching for k1: is more accurate #
	regsub -all $k_val $data "$binomial" data
}
set out [open MRC_master_tree.tre w]
puts $out [string trim $data]
close $out

catch {exec raxml -L MR -z MRC_master_tree.tre -m GTRCAT -n MR_tree.txt -T 4}