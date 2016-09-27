set ival 2.0

set direct /users/aesin/desktop/Geo_v_all/$ival
set 10_eval_fasta_dir $direct/-10/Group_fastas_$ival
set missing_geo_dir $direct/Missing_geo_groups/Geo_only
# Raxml cores #
set cores 2

###########################################################################
### PROCS ###
proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

### Split a fasta file full of proteins into a list of individual proteins ###
proc split_genes {fasta} {
	global genes
	regsub -all {>} [string trim $fasta] {£>} fasta
	set genes [lrange [split $fasta £] 1 end]
    return
}

###########################################################################

## Pull out only the Geobacilli sequences ##

cd $direct/Stats
openfile Missing_geobac_families.tsv
set missing_geo_groups [lrange [split [string trim $data] \n] 1 end]

foreach group $missing_geo_groups {
	set group_no [lindex [split $group \t] 0]
	openfile $10_eval_fasta_dir/$group_no\.faa
	split_genes $data

	set geo_genes [lsearch -all -glob -inline $genes *\(Geobac\)*]
	set geo_genes [join $geo_genes {}]

	set out [open $missing_geo_dir/Fastas/$group_no\.faa w]
	puts $out [string trim $geo_genes]
	close $out
}

cd $missing_geo_dir
file mkdir Fastas_key_ids
file mkdir MSA
file mkdir Tree_build
file mkdir Final_trees

cd $missing_geo_dir/Fastas
set fasta_files [glob *faa]

## Fastas_key_ids and muscle alignment here ##

foreach file $fasta_files {
	openfile $file
	split_genes $data

	set group [string range $file 0 end-4]
	file mkdir $missing_geo_dir/Tree_build/$group
	set i 1
	set new_glist ""
	set key_list ""
	foreach gene $genes {
		set comment_line [string trim [string range $gene 0 [string first \n $gene]]]
		set id [string trim [string range $comment_line 1 [string first " " $comment_line]]]
		set class [string trim [string range $comment_line [string first \( $comment_line]+1 [string first \) $comment_line]-1]]
		set binomial [string range [string trim [string range $comment_line [string last \[ $comment_line] end]] 1 end-1]
		regsub {>+?.+?\n+?} $gene ">k$i\n" new_gene
		append new_glist $new_gene
		append key_list "k$i\t$binomial\t$class\t$id\n"
		incr i
	}

	set out [open $missing_geo_dir/Fastas_key_ids/$file w]
	puts $out [string trim $new_glist]
	close $out

	set out [open $missing_geo_dir/Fastas_key_ids/KEY_$group\.txt w]
	puts $out $key_list
	close $out
	file copy -force $missing_geo_dir/Fastas_key_ids/KEY_$group\.txt $missing_geo_dir/MSA
	file copy -force $missing_geo_dir/Fastas_key_ids/KEY_$group\.txt $missing_geo_dir/Tree_build/$group/

	set aligned_out_name "$group\.fas"
	catch {exec muscle -in $missing_geo_dir/Fastas_key_ids/$file -phyiout $missing_geo_dir/MSA/$aligned_out_name}
	file copy -force $missing_geo_dir/MSA/$aligned_out_name $missing_geo_dir/Tree_build/$group/
}

## Make the trees! ##

cd $missing_geo_dir/Tree_build
set directories [glob -type d *]


set pata {k+?[0-9]+?:+?}
set patb {\t+?.+?\n+?}


foreach family $directories {
	cd $missing_geo_dir/Tree_build/$family
	set align [lindex [glob *fas] 0]
	set out_number $family
	catch {exec raxml -f a -s $align -n $out_number\.txt -m PROTCATAUTO -p 1234 -x 10001 -N 100 -T $cores >@stdout}

	# If a file with BS values is available (small/large) - use that. Otherwise use the best ML tree. #
	set besttrees [glob -nocomplain *bipartitions.*]
	set tree [lindex $besttrees 0]

	set key [lindex [glob KEY*] 0]
	openfile $key; set key $data

	# Relabel the tips, and move the resultant folders/files to the output destinations ONLY if there is a "best tree" #
	if {[llength $tree] > 0} {
		openfile $tree
		set branchids [regexp -all -inline $pata $data]
		foreach k_val $branchids {
			regsub {:} $k_val {} k_val_trunc
			regexp -line $k_val_trunc$patb $key hit
			set binomial [lindex [split [string trim $hit] \t] 1]
			set class [lindex [split [string trim $hit] \t] 2]
			set id [lindex [split [string trim $hit] \t] 3]
			# If you use k_val_trunc below, you may substitute k17 when you want to find k1 - so searching for k1: is more accurate #
			regsub $k_val $data "$binomial \{$class\} \{$id\}:" data
		}
		set out [open $out_number\_tree.txt w]; puts $out $data; close $out
		file copy -force $out_number\_tree.txt $missing_geo_dir/Final_trees/
	}
}














