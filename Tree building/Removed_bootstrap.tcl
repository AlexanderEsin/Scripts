source ~/Dropbox/Scripts/General_utils.tcl

set ival 2
set eval -10
set direct /users/aesin/desktop/Consensus_trees/Inside_group/$eval/Final_trees_$ival

cd $direct/Individual_trees
file mkdir Reference_nbs

set reference_set [glob *txt]

exec java -jar /users/aesin/Dropbox/Scripts/Jar_progs/PareTree.jar -nbs -t O -d $direct/Individual_trees/*txt

set nbs_reference_set [glob *nbs.txt]
foreach nbs_tree $nbs_reference_set {
	file rename -force $nbs_tree Reference_nbs/$nbs_tree

	file copy -force Reference_nbs/$nbs_tree ~/desktop/Test_mincomp/rec

	regsub {.txt} $nbs_tree {.nwk} nwk_ending
	file rename ~/desktop/Test_mincomp/rec/$nbs_tree ~/desktop/Test_mincomp/rec/$nwk_ending
}

## SAMPLE ##

# exec java -jar /users/aesin/Dropbox/Scripts/Jar_progs/PareTree.jar -nbs -t O -f $tree