#!/usr/local/bin/tclsh
###
set eval [lindex $argv 0]
set ival 2.0
set direct /users/aesin/desktop/Geo_v_all
set rbbh_dir Rbbh_all_$eval
###

source ~/Dropbox/Scripts/General_utils.tcl

package require sqlite3

###########################################################################
file mkdir $direct/$eval
### Make a filtered glist here from the mcl groups ###
cd $direct/$eval
if {[file exists Group_fastas_$ival] == 1} {
	file delete -force Group_fastas_$ival
}
file mkdir Group_fastas_$ival

### Open the sqlite database ###
cd /users/aesin/desktop/Clean_proteomes
sqlite3 db1 all_prot_geo_db

### Open the MCL-output groups file ###
cd $direct/$rbbh_dir
openfile Master_groups_weight_$ival\.txt
set glist $data

### Create Group_fastas output folder ###
cd $direct/$eval
file mkdir Group_fastas_$ival
if {[file exists Group_list_$ival\.tsv] == 1} {
	file delete Group_list_$ival\.tsv
}

set g 1
set total_genes 0

### For each group (gene family) pull out the class (Bact/Arch etc) and full gene to create a (filterable) gene list (Group_list.tsv) and fastas for each group (Group_fastas) ###
set genefams [split [string trim $glist] \n]
foreach fam $genefams {
	set genes [split [string trim $fam] \t]
	set total_genes [expr $total_genes + [llength $genes]]
	set i 1
	set desig_fam ""
	foreach gene $genes {
		db1 eval {SELECT class FROM t1 WHERE id = $gene} {
			set class "$class"
		}
		append desig_fam "$gene $class\t"
		db1 eval {SELECT gene FROM t1 WHERE id = $gene} {
			set new_gene "$gene"
		}
		set out [open $direct/$eval/Group_fastas_$ival/$g.faa a]
		puts $out $new_gene
		close $out
		puts "Group: $g/[llength $genefams] === $i/[llength $genes]"
		incr i
	}
	set out [open $direct/$eval/Group_list_$ival\.tsv a]
	puts $out "$g\t$desig_fam"
	close $out

	incr g
}

set total_genes [commas $total_genes]
if {[file exists $direct/$rbbh_dir/Rbbh_stats/Edge_number_av_eval.txt] == 1} {
	set out [open $direct/$rbbh_dir/Rbbh_stats/Edge_number_av_eval.txt a]
	puts $out "Total number of genes in groups: $total_genes"
	close $out
}
###########################################################################
### Convert the binomial names for each gene into a phylip-friendly format + fix MULTISPECIES name tags ###
### All the regsubs attempt to remove random crap from the binomial name to avoid complications downstream ###

cd $direct/$eval
file mkdir temp

cd $direct/$eval/Group_fastas_$ival

set fams [glob *.faa]
set fams [lsort -dictionary $fams]
set i 1

foreach fam $fams {
	openfile $fam
	set new_glist {}
	split_genes $data
	foreach gene $genes {
		set gene [string trim $gene]
		regsub -all {\*} $gene {} gene
		set binomial [string range [string range $gene 0 [string first \n $gene]-1] [string last \[ $gene]+1 end-1]
		### If there is a MULTISPECIES tag, then query the all_prot_geo_db with the gene_id to get the binomial of the organism ###
		if {[regexp {MULTISPECIES} $gene] == 1} {
			regexp {>+?.+?\s+?} $gene query
			set query [string trim [string range $query 1 end]]
			db1 eval {SELECT binomial FROM t1 WHERE id = $query} {
				set binomial "$binomial"
			}
			# The [binomial name] occasionally contains (bracketed terms) #
			regsub -all {\(} $binomial {} binomial
			regsub -all {\)} $binomial {} binomial
			regsub -all {\:} $binomial {} binomial
			regsub -all {\#} $binomial {} binomial
			regsub -all {\=} $binomial {} binomial
			regsub -all {\/} $binomial {} binomial
			regsub -all {'} $binomial {} binomial
			regsub -all " " $binomial {_} match
			regsub -all {\.} $match {} match
			regsub -all {\,} $match {} match
			regsub -all {\+} $match {} match
			regsub -all -- {\-} $match {_} match
			# Use the full reg syntax because if there are (curly brackets), then using a variable (e.g. $hit) will not work #
			set header_trim [string range $gene 0 [string last \[ $gene]-1]
			set new_header "$header_trim\[$match\]"
			set gene2 "$new_header\n[string trim [string range $gene [string first \n $gene] end]]"
			lappend new_glist $gene2
		} else {
			# The [binomial name] occasionally contains (bracketed terms) #
			regsub -all {\(} $binomial {} binomial
			regsub -all {\)} $binomial {} binomial
			regsub -all {\:} $binomial {} binomial
			regsub -all {\#} $binomial {} binomial
			regsub -all {\=} $binomial {} binomial
			regsub -all {\/} $binomial {} binomial
			regsub -all {'} $binomial {} binomial
			regsub -all " " $binomial {_} match
			regsub -all {\.} $match {} match
			regsub -all {\,} $match {} match
			regsub -all {\+} $match {} match
			regsub -all -- {\-} $match {_} match
			set header_trim [string range $gene 0 [string last \[ $gene]-1]
			set new_header "$header_trim\[$match\]"
			set gene2 "$new_header\n[string trim [string range $gene [string first \n $gene] end]]"
			lappend new_glist $gene2
		}
	}
	### Output the new fasta file to a temp directory (originally for debugging) ###
	set new_glist_string [string trim [join $new_glist \n]]
	set out [open $direct/$eval/temp/$i.faa w]
	puts $out $new_glist_string
	close $out
	puts "Done: $i/[llength $fams]"
	incr i
}

### Convert the temp directory into Group_fastas, and leave that as is (for checking downstream). Create a copy which will be manipulated further (Group_fastas_key_ids) ###
cd $direct/$eval
file delete -force Group_fastas_$ival
file rename -force temp Group_fastas_$ival
file copy Group_fastas_$ival Group_fastas_key_ids_$ival

###########################################################################

db1 close

puts "DONE"
puts "Total number of genes in groups: $total_genes"
