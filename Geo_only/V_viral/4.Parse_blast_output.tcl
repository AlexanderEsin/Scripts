set direct /users/aesin/desktop/Geo_v_all/2.0
set org Geo_only_groups; #Anoxy_geo_only_groups or Geo_only_groups
set blast_v_what Blast_v_Gthermo
set org $org/$blast_v_what
 
proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct/$org/Out_all

set blast_output_files [glob *tsv]
set patc {\#+?.+?\n+?}
set output_table "Group\tTotal genes\tNumber of hits\tE-values"

foreach file $blast_output_files {
	openfile $file
	set data [string trim $data]

	set nohit [regexp -all "0 hits found" $data]
	set tot_genes [regexp -all "hits found" $data]

	regsub -all $patc $data {} data
	set lineparse [string trim [split [string trim $data] \n]]

	if {$nohit != $tot_genes} {
		set evals {}
		set group_no [lindex [split $file &] 0]
		set hit_no [expr $tot_genes - $nohit]

		foreach line $lineparse {
			set evalue [lindex [split $line \t] end-1]
			lappend evals $evalue
		}
		set evals_str [join $evals " "]
		append output_table "\n$group_no\t$tot_genes\t$hit_no\t$evals_str"
	}
}

set out [open $direct/$org/Significant_viral_hits.txt w]
puts $out [string trim $output_table]
close $out