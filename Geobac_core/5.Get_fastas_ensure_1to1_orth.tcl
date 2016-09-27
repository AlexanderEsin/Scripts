##
set direct /users/aesin/desktop/Consensus_trees
set org Geobac/All_geobac
set eval -10
set ival 5
set rbbh_dir Rbbh_all_$eval
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

proc split_genes {fasta} {
	global genes
	regsub -all {>} [string trim $fasta] {£>} fasta
	set fasta [split $fasta £]
	regsub -all "{}" $fasta {} genes
}

package require sqlite3

###########################################################################
### Make a filtered glist here from the mcl groups ###
cd $direct/$org
if {[file exists Group_fastas_$ival] == 1} {
	file delete -force Group_fastas_$ival
}
file mkdir Group_fastas_$ival

### Open the MCL-output groups file ###
cd $direct/$org/$rbbh_dir
openfile Master_groups_weight_$ival\.txt
set glist $data

### Open the prot database & create Group_fastas output folder ###
cd /users/aesin/desktop/Clean_proteomes
sqlite3 db1 all_prot_geo_db_2
###

set g 1

### For each group (gene family) pull out the class (Bact/Arch etc) and full gene to create a (filterable) gene list (Group_list.tsv) and fastas for each group (Group_fastas) ###
set genefams [split $glist \n]
regsub -all "{}" $genefams {} genefams
foreach fam $genefams {
	set genes [split $fam \t]
	## REMEMBER TO CHANGE TO NUMBER OF TOTAL ORGANISMS ##
	if {[llength $genes] > 19} {
		puts "Too long"
	} elseif {[llength $genes] < 19} {
		puts "Too short"
	} else {
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
			set out [open $direct/$org/Group_fastas_$ival/$g.faa a]
			puts $out $new_gene
			close $out
			puts "Group: $g/[llength $genefams] === $i/[llength $genes]"
			incr i
		}
		set out [open $direct/$org/Group_list_$ival.tsv a]
		puts $out "$g\t$desig_fam"
		close $out
	}
	incr g
}

###########################################################################
### Convert the binomial names for each gene into a phylip-friendly format + fix MULTISPECIES name tags ###
### All the regsubs attempt to remove random crap from the binomial name to avoid complications downstream ###

cd $direct/$org
file mkdir temp

cd $direct/$org/Group_fastas_$ival

set fams [glob *.faa]
set fams [lsort -dictionary $fams]
set i 1

set delete_list ""

foreach fam $fams {
	set break_flag 0
	set fam_binomials {}

	openfile $fam
	set new_glist ""
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
			set binomial [string range $binomial 1 end-1]
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
			append new_glist "$gene2\n"
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
			append new_glist "$gene2\n"
		}
		## Weed out families where there is more than one copy of the same binomial ##
		if {[lsearch $fam_binomials $match] == -1} {
			lappend fam_binomials $match
		} else {
			set break_flag 1
			set binomial $match
		}
		
	}
	### Output the new fasta file to a temp directory (originally for debugging) ###
	if {$break_flag == 0} {
		set new_glist [string trim $new_glist]
		set out [open $direct/$org/temp/$i.faa w]
		puts $out $new_glist
		close $out
		puts "Done: $i/[llength $fams]"
		incr i
	} else {
		append delete_list "$fam\t$binomial\n"
		incr i
	}
}

set out [open $direct/$org/Ortholog_groups_not1to1.txt w]
puts $out [string trim $delete_list]
close $out

### Convert the temp directory into Group_fastas, and leave that as is (for checking downstream). Create a copy which will be manipulated further (Group_fastas_key_ids) ###
cd $direct/$org
file rename -force temp Group_fastas_key_ids_$ival

###########################################################################
db1 close
puts "DONE"

