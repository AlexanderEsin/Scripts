set direct ~/desktop/Geobac/Indiv_genes

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct
file mkdir Best_trees_nbs
file mkdir C2T/Good_trees

cd $direct/Best_trees_bootval
set trees [glob *.txt]

foreach tree $trees {
	catch {exec java -jar /users/aesin/desktop/Scripts/Jar_progs/PareTree.jar -f $tree -nbs}
	set temp [glob *nbs*]
	file rename $temp $direct/Best_trees_nbs/
}

###########################################################################

cd ~/desktop/Geobac/Filtered_geobac_fastas/Best_bootstrap_tree

file copy Geobac_bootstrap_tree.txt $direct
cd $direct
catch {exec java -jar /users/aesin/desktop/Scripts/Jar_progs/PareTree.jar -f Geobac_bootstrap_tree.txt -nbs}
file rename -force Geobac_bootstrap_tree_nbs.txt $direct/Best_trees_nbs/Reference_tree.tree
file delete Geobac_bootstrap_tree.txt

###########################################################################

cd $direct/Best_trees_nbs
set i 1
set compare_scores "Tree\tCompare_score\n"

set trees [glob *.txt]
foreach tree $trees {
	set x 0
	puts "$i/[llength $trees]"
	set x [exec java -Dapple.awt.UIElement=true -classpath /users/aesin/desktop/Scripts/Jar_progs/PhyloCore.jar:parallelcolt-0.9.4.jar treecomparison.Compare2Trees -f $tree Reference_tree.tree]]
	regexp {[0-9]+\.+[0-9]+} $x x
	if {$x == "Score = 1]"} {
		set x 1
	}
	if {$x > 0.95} {
		regsub "_nbs" $tree {} orig_tree
		file copy $direct/Best_trees_bootval/$orig_tree $direct/C2T/Good_trees/
	}
	incr i
	puts $x
	append compare_scores "$tree\t$x\n"
}

set compare_scores [string trim $compare_scores]
set out [open $direct/Comparison_scores_C2T.tsv w]
puts $out $compare_scores
close $out













