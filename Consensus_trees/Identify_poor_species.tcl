#!/usr/local/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

##
set direct /users/aesin/desktop/Consensus_trees
set org Geobac/Outgroups/With_excl_geos
set eval -10
set ival 2
set rbbh_dir Rbbh_all_$eval
###

###########################################################################

cd /users/aesin/desktop/Clean_proteomes
sqlite3 db1 all_prot_geo_ncbi_db

cd $direct/$org/$eval

openfile Group_list_$ival\.tsv
set sample_group [string trim [lindex [split [string trim $data] \n] 0]]
set scgs [lrange [split $sample_group \t] 1 end]

set bacillaceae_orgs {}
set i 1
foreach gene $scgs {
	set gene [string trim [string range $gene 0 [string first " " $gene]]]
	db1 eval {SELECT binomial FROM t1 WHERE id = $gene} {
		set binomial "$binomial"
	}
	lappend bacillaceae_orgs $binomial
	puts "$i / [llength $scgs]"
	incr i
}

###########################################################################

### Open the MCL-output groups file ###
cd $direct/$org/$rbbh_dir
openfile Master_groups_weight_$ival\.txt
set genefams [split [string trim $data] \n]

set absent_orgs {}

foreach fam $genefams {
	set genes [split $fam \t]
	## REMEMBER TO CHANGE TO NUMBER OF TOTAL ORGANISMS ##
	if {[llength $genes] < 202 && [llength $genes] > 195} {
		set test_fam_genes {}
		foreach gene $genes {
			db1 eval {SELECT binomial FROM t1 WHERE id = $gene} {
				set binomial "$binomial"
			}
			lappend test_fam_genes $binomial
		}
		foreach binomial $bacillaceae_orgs {
			if {[lsearch $test_fam_genes $binomial] == -1} {
				lappend absent_orgs $binomial
			}
		}
	} else {
		puts "Wrong size"
	}
}

set absent_orgs [join $absent_orgs "\n"]

###

db1 close

###

set orgs [split $absent_orgs \n]
set org_miss_no_list ""
set done_list {}
set i 0
foreach organ $orgs {
	if {[lsearch $done_list $organ] == -1} {
		set miss_no [regexp -all $organ $data]
		append org_miss_no_list "$organ\t$miss_no\n"
		lappend done_list $organ
		incr i
	}
}
puts "Total poor species: $i"
puts $org_miss_no_list

set out [open $direct/$org/$eval/Absent_species_$eval\_$ival\.tsv w]
puts $out [string trim $org_miss_no_list]
close $out


