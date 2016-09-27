#!/usr/local/bin/tclsh
###########################################################################
## Change Log ##
# 21 January 2016 : Replaced procs with sourced procs. Changed vars to lists. #
# 22 January 2016 : Attempted to fix this script. Unfortunately, the taxonomy db can only be interfaced using the taxa names, which are difficult to extract from the fasta files. Since Taxonomy_all_2.tsv already exists and seems correct, I deciced to not waste time on re-doing this. Added switches, so that a database could be made without other faff. #

###########################################################################
source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

set direct /users/aesin/desktop/Clean_Proteomes
set make_taxa_names_switch 0
set extract_ranks_switch 0
set make_db_switch 1

###########################################################################
if {make_taxa_names_switch == 1} {
	cd $direct/All_fastas_geo_mark

	### Make a list of the all the fasta files in the folder ###
	set files [glob *faa]
	set files [lrange [glob *faa] 0 99]

	set done_binomials {}
	set names ""
	set i 1

	## This is really a bit of a mess. It's difficult to sort out - but the file Taxonomy_all_2.tsv contains all the necessary info ##

	foreach file $files {

		regsub -all "_" $file " " new_name
		set new_name [join [lrange [split $new_name " "] 0 1]]

		set genes [split_genes_fast $data]
		set rev_genes [lreverse $genes]
		set binomial ""
		foreach gene $rev_genes {
			if {[string first "MULTISPECIES" $gene] == -1} {
				set test_binomial [string range [string range $gene 0 [string first \n $gene]-1] [string last \[ $gene]+1 end-1]
				if {$test_binomial == $new_name} {
					set binomial $test_binomial
					break
				}
			}
		}

		if {[string length $binomial] == 0} {
			foreach gene $rev_genes {
				if {[string first "MULTISPECIES" $gene] == -1} {
					set test_binomial [string range [string range $gene 0 [string first \n $gene]-1] [string last \[ $gene]+1 end-1]
					if {[string first $binomial $done_binomials] == -1} {
						set binomial $test_binomial
						break
					}
				}
			}
		}
		if {[string length $binomial] == 0} {
			puts "BINOMIAL NOT FOUND == $file"
			error
		}
		lappend done_binomials $binomial
		append names "$binomial\n"
		puts "$i / [llength $files]"
		incr i
	}

	set name_list [split [string trim $names] \n]
	puts [llength $name_list]

	set out [open $direct/Taxa_names.tsv w]
	puts $out [string trim $names]
	close $out
}

###########################################################################
if {$extract_ranks_switch == 1} {
	exec >&@stdout Rscript /users/aesin/desktop/Scripts/Proteome_preprocess/Extract_tax_ranks.R 
}

###########################################################################
## Making this database is only useful for extracting all the fasta files belonging to a certain taxanomic group from the complete list of fastas used ##

if {$make_db_switch == 1} {
	cd $direct

	sqlite3 db1 Taxonomic_ranks_db
	db1 eval {PRAGMA main.page_size = 4096}
	db1 eval {PRAGMA main.cache_size=10000}
	db1 eval {PRAGMA main.locking_mode=EXCLUSIVE}
	db1 eval {PRAGMA main.synchronous=NORMAL}
	db1 eval {PRAGMA main.journal_mode=WAL}
	db1 eval {create table t1(binomial text, ranks blob)}

	openfile Taxonomy_all_2.tsv
	set groups [split [string trim $data] \n]
	foreach group $groups {
		set ranks [split [string trim $group] \t]
		set binomial [lindex $ranks 0]
		set other_ranks [join [lrange $ranks 1 end] \t]
		regsub -all {\*} $binomial {} match
		regsub -all {\(} $match {} match
		regsub -all {\)} $match {} match
		regsub -all {\:} $match {} match
		regsub -all {\#} $match {} match
		regsub -all {\=} $match {} match
		regsub -all {\/} $match {} match
		regsub -all {'} $match {} match
		regsub -all " " $match {_} match
		regsub -all {\.} $match {} match
		regsub -all {\,} $match {} match
		regsub -all {\+} $match {} match
		regsub -all -- {\-} $match {_} binomial
		
		db1 eval {insert into t1 values($binomial,$other_ranks)}

		regsub -all {[][*():#=/',+]} $binomial {} new
		regsub -all -- {[-/]} $new " " new

	}

	db1 close
}


###########################################################################
