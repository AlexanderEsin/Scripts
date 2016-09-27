set direct ~/desktop/Geobac

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct/Rbbh/Bin; ####Rbbh_all or Rbbh?

set rbbh_list [glob *.txt]
set i 1

foreach file $rbbh_list {
	openfile $file
	set rhead [string first \n $data]
	set rbbhs [string range $data $rhead end]
	append rbbhs "\n"
	set rbb1 [string trim $rbbhs]
	set out [open $direct/Rbbh/Master_rbh_weight a]; ###Same here
	puts $out $rbb1
	close $out
	puts "$i/[llength $rbbh_list]"
	incr i
}