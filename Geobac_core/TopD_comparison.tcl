set direct /users/aesin/desktop/Geobac/Indiv_genes
set org Trees

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct
openfile TopD_template.txt
set master [string trim $data]

cd $direct/$org
set dirs [glob -type d *]
set dirs [lsort -dictionary $dirs]

foreach dir $dirs {
	cd $direct/$org/$dir
	set tree [lindex [glob *bipartitions.*] 0]
	openfile $tree
	set tree_data "$dir\t[string trim $data]"
	append master "\n$tree_data"
}

set out [open $direct/TopD_comp.txt w]
puts $out $master
close $out

./topd.pl -f /users/aesin/desktop/Geobac/Indiv_genes/TopD_comp.txt -c reference -r no -p n -m split

openfile topd_results.txt
set mess $data
regsub -all {\n\n} $mess "£" mess
set results [split $mess £]
regsub -all "{}" $results {} results

set interest_list ""

foreach result $results {
	set top_line [string trim [string range $result 0 [string first \n $result]]]
	regexp {[0-9]+} $top_line file_no

	set res_line [lindex [split $result \n] end]
	set splits [string trim [string range $res_line [string last \[ $res_line] [string last / $res_line]]]
	regexp {[0-9]+} $splits split_no

	if {$split_no > 20} {
		append interest_list "$file_no\t$split_no\n"
	}
}

set out [open candidates.tsv w]
puts $out [string trim $interest_list]
close $out


#### This is the topd command to be run from bash ####

#./topd.pl -f /users/aesin/desktop/Geobac/Indiv_genes/TopD_comp.txt -c reference -r no -p n -m split


###########################################################################

##### TESTING #####



cd $direct

set fam 632.fas
set fname TEST

exec raxml -f d -s $fam -n $fname -m PROTCATAUTO -p 1234 -N 10 -T 2

set red [lindex [glob -nocomplain *reduced] 0]
file delete $red

set raxs [glob RA*]
foreach rax $raxs {
	file delete $rax
}
