set direct /users/aesin/desktop/Consensus_trees
set eval -100
set rbbh_dir Rbbh_all_$eval
set org Geobac/All_geobac/$rbbh_dir

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

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
