###
set eval -100
set ival 2.5
set direct /users/aesin/desktop/Geo_v_all/$eval
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

proc reverse_dict {lst} {
	global revdict
	while {[llength $lst] > 0} {
		set element [lindex $lst end]
		lappend revdict $element
		set lst [lrange $lst 0 end-1]
	}
}

proc split_genes {fasta} {
	global genes
	regsub -all {>} [string trim $fasta] {£>} fasta
	set genes [lrange [split $fasta £] 1 end]
}

###########################################################################
### Get all groups that have at least two Geobacilli and one Archaeon ###

# Split into two folders - those which will have an outgroup and those that won't #
cd $direct
file mkdir Geo_arch_analysis_$ival/Geobac_archaea_groups
file mkdir Geo_arch_analysis_$ival/Geobac_archaea_fastas/outgroup
file mkdir Geo_arch_analysis_$ival/Geobac_archaea_fastas/no_outgroup
file mkdir Geo_arch_analysis_$ival/Geobac_archaea_aligns/outgroup
file mkdir Geo_arch_analysis_$ival/Geobac_archaea_aligns/no_outgroup
file mkdir Geo_arch_analysis_$ival/Reference_for_monophyly
file mkdir Geo_arch_analysis_$ival/Results_not_mono

cd $direct
openfile Group_list_$ival.tsv

set glist [split [string trim $data] \n]
regsub -all "{}" $glist {} glist

set pata {\(+?Bact+?\)+?}
set patb {\(+?Arch+?\)+?}
set patg {\(+?Geobac+?\)+?}

set geobac_archaea ""
set i 1

foreach fam $glist {
	# The groups need to have at least two Geobacilli (patg) and two Archaea (patb) #
	if {[regexp -all $patg $fam] > 1 && [regexp -all $patb $fam] > 1 && [expr ([regexp -all $patg $fam] + [regexp -all $patb $fam])] > 3} {
		# For each such group transfer the corresponding alignment and KEY into Geobac_archaea_groups folder #
		append geobac_archaea "$fam\n"
		set group_number [string trim [string range $fam 0 [string first \t $fam]]]
		file mkdir $direct/Geo_arch_analysis_$ival/Geobac_archaea_groups/$group_number
		file copy -force $direct/Group_fastas_MSA_$ival/$group_number\.fas $direct/Geo_arch_analysis_$ival/Geobac_archaea_groups/$group_number/
		file copy -force $direct/Group_fastas_MSA_$ival/KEY_$group_number\.txt $direct/Geo_arch_analysis_$ival/Geobac_archaea_groups/$group_number/
	
		# Open the corresponding fasta file from Group_fastas and filter out the Geobacilli, Archaea, and a possible outgroup. These need to go in a separate file that can be used by is.monophyletic to assign tip labels. If there is no outgroup #
		openfile $direct/Group_fastas_$ival/$group_number\.faa
		split_genes $data

		set archaea ""
		set geobac ""

		set pull_out_MSA_list ""
		set geo_arch_only_faa ""

		set data "[string trim $data]\n>"
		set patx {.+?>+?}

		# Pull out the outgroup (subtilis) if present #
		if {[string first "Bacillus_subtilis_subsp_subtilis_str_168" $data] > -1} {
			set temp_string [string range $data 0 [string first "Bacillus_subtilis_subsp_subtilis_str_168" $data]]
			set temp_string [string range $temp_string [string last ">" $temp_string] end]
			set id [string trim [string range $temp_string 0 [string first " " $temp_string]]]
			regexp $id$patx $data outlier
			set outlier [string trim [string range $outlier 0 end-1]]
			set id [string range $id 1 end]

			set outgroup "Bacillus_subtilis_subsp_subtilis_str_168 \{$id\}"
			append pull_out_MSA_list "$id\t"
			append geo_arch_only_faa "$outlier\n"

			set out_flag 1
		} else {
			set outgroup "-"
			set out_flag 0
		}

		# Parse through all the genes in each group and pull out all Geobacilli and all Archaea #
		foreach gene $genes {
			if {[regexp $patb $gene] == 1} {
				append geo_arch_only_faa "[string trim $gene]\n"

				set binomial [string range [string range $gene 0 [string first \n $gene]-1] [string last \[ $gene]+1 end-1]
				set id [string trim [string range $gene 1 [string first " " $gene]]]

				append archaea "$binomial \{$id\}\t"

				append pull_out_MSA_list "$id\t"
			} elseif {[regexp $patg $gene] == 1} {
				append geo_arch_only_faa "[string trim $gene]\n"

				set binomial [string range [string range $gene 0 [string first \n $gene]-1] [string last \[ $gene]+1 end-1]
				set id [string trim [string range $gene 1 [string first " " $gene]]]

				append geobac "$binomial \{$id\}\t"

				append pull_out_MSA_list "$id\t"
			}
		}

		# Create the reference file for the is.monophyletic query later #
		set reference "Geobacilli:\t[string trim $geobac]\nArchaea:\t[string trim $archaea]\nOutgroup:\t[string trim $outgroup]"
		set out [open $direct/Geo_arch_analysis_$ival/Reference_for_monophyly/Reference_$group_number.txt w]
		puts $out [string trim $reference]
		close $out

		# Pull out fasta files with only the Geobacilli, Archaea, and possible outgroup sequences for corss-reference into Geobac_archaea_fastas #
		if {$out_flag == 1} {
			set out [open $direct/Geo_arch_analysis_$ival/Geobac_archaea_fastas/outgroup/$group_number\.faa w]
		} else {
			set out [open $direct/Geo_arch_analysis_$ival/Geobac_archaea_fastas/no_outgroup/$group_number\.faa w]
		}
		puts $out [string trim $geo_arch_only_faa]
		close $out
		
		# Take the relevant MSA out of interleaved phylip into aligned fasta format to make extraction of correct sequences easier #
		cd $direct/Geo_arch_analysis_$ival/Geobac_archaea_groups/$group_number/
		set fastas [glob -nocomplain *.fasta]
		if {[llength $fastas] == 0} {
			catch {exec seqret $direct/Geo_arch_analysis_$ival/Geobac_archaea_groups/$group_number/$group_number\.fas $direct/Geo_arch_analysis_$ival/Geobac_archaea_groups/$group_number/$group_number\.fasta -osformat fasta}
		}
		
		# Get all the relevant k_values from the KEY
		set patc {.+?\t+?}
		set k_val_list ""

		set pull_out [split [string trim $pull_out_MSA_list] \t]
		openfile $direct/Group_fastas_MSA_$ival/KEY_$group_number\.txt
		set key_data [string trim $data]
		foreach gene $pull_out {
			regexp -line $patc$gene $key_data key_line
			set k_value [lindex [split [string trim $key_line] \t] 0]
			append k_val_list "$k_value\t"
		}
		set k_val_list [string trim $k_val_list]

		# Extract the needed aligned sequences from the full alignment (now in fasta format) into a new alignment #

		openfile $direct/Geo_arch_analysis_$ival/Geobac_archaea_groups/$group_number/$group_number\.fasta
		set parent_align "[string trim $data]\n>"

		set daughter_align ""

		set paty {>}
		set patz {\n+?.+?>+?}
		set k_vals [split $k_val_list \t]
		set new_align_counter 0
		foreach k_val $k_vals {
			# Since we are working with the trimmed "core" alignments, certain (very diverged sequences) may have only contained gaps over this region and so would have been removed in the Finalise_core_split_for_trees.tcl script. Obviously these will not be found for pull_out in the alignment. #
			if {[regexp $paty$k_val$patz $parent_align align] == 1} {
				set align [string trim [string range $align 0 end-1]]
				append daughter_align "$align\n"
				incr new_align_counter
			}
		}

		if {$new_align_counter < 4} {
			set out_flag -1
		}

		# For each daughter alignment, output first as a fasta alignment then convert into an interleaved phylip format into the correct folder #
		if {$out_flag == 1} {
			set out [open $direct/Geo_arch_analysis_$ival/Geobac_archaea_aligns/outgroup/$group_number\.faa w]
			puts $out [string trim $daughter_align]
			close $out

			catch {exec seqret $direct/Geo_arch_analysis_$ival/Geobac_archaea_aligns/outgroup/$group_number\.faa $direct/Geo_arch_analysis_$ival/Geobac_archaea_aligns/outgroup/$group_number\.fas -osformat phylip}

			file delete -force $direct/Geo_arch_analysis_$ival/Geobac_archaea_aligns/outgroup/$group_number\.faa
			file copy -force $direct/Group_fastas_MSA_$ival/KEY_$group_number\.txt $direct/Geo_arch_analysis_$ival/Geobac_archaea_aligns/outgroup/
		} elseif {$out_flag == 0} {
			set out [open $direct/Geo_arch_analysis_$ival/Geobac_archaea_aligns/no_outgroup/$group_number\.faa w]
			puts $out [string trim $daughter_align]
			close $out

			catch {exec seqret $direct/Geo_arch_analysis_$ival/Geobac_archaea_aligns/no_outgroup/$group_number\.faa $direct/Geo_arch_analysis_$ival/Geobac_archaea_aligns/no_outgroup/$group_number\.fas -osformat phylip}

			file delete -force $direct/Geo_arch_analysis_$ival/Geobac_archaea_aligns/no_outgroup/$group_number\.faa
			file copy -force $direct/Group_fastas_MSA_$ival/KEY_$group_number\.txt $direct/Geo_arch_analysis_$ival/Geobac_archaea_aligns/no_outgroup/
		}
	}
	puts "$i / [llength $glist]"
	incr i
}

set out [open $direct/Geo_arch_analysis_$ival/Geobac_archaea_groups.tsv w]
puts $out $geobac_archaea
close $out

###########################################################################
### Run RAxML on the Geobacilli-Archaeal trees ###
cd $direct
file mkdir Geo_arch_analysis_$ival/Final_trees/outgroup
file mkdir Geo_arch_analysis_$ival/Final_trees/no_outgroup
file mkdir Geo_arch_analysis_$ival/Tree_build

cd $direct/Geo_arch_analysis_$ival/Geobac_archaea_aligns/

set directories [glob -type d *]

foreach directory $directories {
	cd $directory
	set aligns [glob *.fas]
	set aligns [lrange $aligns 96 end]

	foreach align $aligns {
		regsub {.fas} $align {} out_number
		catch {exec raxml -f d -s $align -n $out_number\.txt -m PROTCATAUTO -p 1234 -N 2 -T 4}

		openfile "KEY_$out_number\.txt"
		set key $data

		set pata {k+?[0-9]+?:+?}
		set patb {\t+?.+?\n+?}

		set besttrees [glob -nocomplain *bestTree.*]
		if {[llength $besttrees] > 0} {
			foreach tree $besttrees {
				openfile $tree
				set branchids [regexp -all -inline $pata $data]
				foreach k_val $branchids {
					regsub {:} $k_val {} k_val_trunc
					regexp -line $k_val_trunc$patb $key hit
					set binomial [lindex [split [string trim $hit] \t] 1]
					set id [lindex [split [string trim $hit] \t] 2]
					regsub $k_val $data "$binomial \{$id\}:" data
				}
				set out [open $out_number\_tree.txt w]
				puts $out $data
				close $out
				file copy $out_number\_tree.txt $direct/Geo_arch_analysis_$ival/Final_trees/$directory/
			}
		}
		file mkdir $direct/Geo_arch_analysis_$ival/Tree_build/$out_number
		set raxmls [glob RA*]
		foreach raxml $raxmls {
			file rename $raxml $direct/Geo_arch_analysis_$ival/Tree_build/$out_number/
		}
		file rename $out_number\_tree.txt $direct/Geo_arch_analysis_$ival/Tree_build/$out_number/
	}
	set reds [glob *.reduced]
	foreach red $reds {
		file delete $red
	}
	cd ..
}

###########################################################################
### Execute an R-script to test the monophyly of either the Geobacilli+outgroup or the Archaea depending on whether there was an outgroup initiall ###
# A proc to create an R vector out of a list #
proc R.c {args} {
	return c([join $args ,])
}

cd $direct/Geo_arch_analysis_$ival/Final_trees

set directories [glob -type d *]

foreach directory $directories {
	cd $direct/Geo_arch_analysis_$ival/Final_trees/$directory
	set trees [glob *.txt]
	set result_final ""
	set result_not_mono ""
	foreach tree $trees {
		set group_number [string range $tree 0 end-9]
		openfile $direct/Geo_arch_analysis_$ival/Reference_for_monophyly/Reference_$group_number\.txt
		set data "[string trim $data]\n"
		# If there is an outgroup we want to concat the Geobacilli labels with the outgroup label, then process into an R-friendly vector #
		if {$directory == "outgroup"} {
			regexp {(Geobacilli:)+?.+?\n+?} $data tips1
			set geobacs [lrange [split [string trim $tips1] \t] 1 end]
			regexp {(Outgroup:)+?.+?\n+?} $data tips2
			set outgroup [lrange [split [string trim $tips2] \t] 1 end]
			set tips_pre [concat $geobacs $outgroup]
			set tips ""
			foreach tip $tips_pre {
				regsub -all " " $tip {} tip
				set tip "\"$tip\","
				append tips $tip
			}
			set tips [string range $tips 0 end-1]
			set tips [R.c $tips]
		# Otherwise we test 
		} else {
			regexp {(Archaea:)+?.+?\n+?} $data tips1
			set tips_pre [lrange [split [string trim $tips1] \t] 1 end]
			set tips ""
			foreach tip $tips_pre {
				regsub -all " " $tip {} tip
				set tip "\"$tip\","
				append tips $tip
			}
			set tips [string range $tips 0 end-1]
			set tips [R.c $tips]
		}
		###

		set out [open $direct/Geo_arch_analysis_$ival/Test_monophyly.R w]
		puts $out "library \(ape\)\nsetwd\(\"/users/aesin/desktop/Geo_v_all/-100/Geo_arch_analysis_$ival/Final_trees/$directory/\"\)\ntree1<-read.tree\(\"$tree\"\)\nis.monophyletic\(tree1, tips = $tips, reroot = TRUE\)"
		close $out

		file attributes $direct/Geo_arch_analysis_$ival/Test_monophyly.R -permissions 0777

		set results [exec Rscript $direct/Geo_arch_analysis_$ival/Test_monophyly.R]
		set results [lindex $results 1]
		append result_final "$tree\t$results\n"
		if {[string trim $results] == "FALSE"} {
			append result_not_mono "$tree\n"
		}
	}
	set out [open $direct/Geo_arch_analysis_$ival/Results/$directory\_full.txt w]
	puts $out [string trim $result_final]
	close $out

	set out [open $direct/Geo_arch_analysis_$ival/Results/$directory\_not_mono.txt w]
	puts $out [string trim $result_not_mono]
	close $out

	cd ..
}

###########################################################################

cd $direct
file mkdir ToI_$ival

cd $direct/Trees_$ival
#openfile $direct/Geo_arch_analysis_$ival/Results/outgroup_not_mono.txt
#set trees_not_mono $data
openfile $direct/Geo_arch_analysis_$ival/Results/no_outgroup_not_mono.txt
set trees_not_mono [string trim $data]

set trees_not_mono [split [string trim $trees_not_mono] \n]
set dirs [glob -type d *]

foreach tree $trees_not_mono {
	set group_number [string range $tree 0 end-9]
	cd $direct/Trees_$ival/[lindex $dirs 0]
	if {[file exists "$tree"] == 1} {
		if {[file exists $direct/ToI_$ival/$group_number/$tree] == 0} {
			file mkdir $direct/ToI_$ival/$group_number
			file copy $direct/Tree_build_$ival/[lindex $dirs 0]/$group_number/RAxML_bestTree.$group_number\.txt  $direct/ToI_$ival/$group_number/
			continue
		}
	} else {
		cd $direct/Trees_$ival/[lindex $dirs 1]
	}
	if {[file exists "$tree"] == 1} {
		if {[file exists $direct/ToI_$ival/$group_number/$tree] == 0} {
			file mkdir $direct/ToI_$ival/$group_number
			file copy $direct/Tree_build_$ival/[lindex $dirs 1]/$group_number/RAxML_bestTree.$group_number\.txt $direct/ToI_$ival/$group_number/
			continue
		}
	} else {
		cd $direct/Trees_$ival/[lindex $dirs 2]
	}
	if {[file exists "$tree"] == 1} {
		if {[file exists $direct/ToI_$ival/$group_number/$tree] == 0} {
			file mkdir $direct/ToI_$ival/$group_number
			file copy $direct/Tree_build_$ival/[lindex $dirs 2]/$group_number/RAxML_bestTree.$group_number\.txt $direct/ToI_$ival/$group_number/
		}
	} else {
		puts "$tree not yet done"
	}
}

cd $direct/ToI_$ival

set dirs [glob -type d *]
foreach dir $dirs {
	file copy -force $direct/Group_fastas_MSA_$ival/$dir\.fas $direct/ToI_$ival/$dir/
	file copy -force $direct/Group_fastas_MSA_$ival/KEY_$dir\.txt $direct/ToI_$ival/$dir/
	cd $direct/ToI_$ival/$dir
	catch {exec raxml -s $dir\.fas -n $dir\.txt -m PROTCATAUTO -x 10001 -p 1234 -N 100 -T 6}
	catch {exec raxml -f b -m PROTCATAUTO -z RAxML_bootstrap.$dir\.txt -t RAxML_bestTree.$dir\.txt -n BS_tree_$dir\.txt -T 4}

	openfile "KEY_$dir\.txt"
	set key $data

	set pata {k+?[0-9]+?:+?}
	set patb {\t+?.+?\n+?}

	set besttrees [glob -nocomplain *bipartitions.*]
	if {[llength $besttrees] > 0} {
		foreach tree $besttrees {
			openfile $tree
			set branchids [regexp -all -inline $pata $data]
			foreach k_val $branchids {
				regsub {:} $k_val {} k_val_trunc
				regexp -line $k_val_trunc$patb $key hit
				set binomial [lindex [split [string trim $hit] \t] 1]
				set class [lindex [split [string trim $hit] \t] 2]
				set id [lindex [split [string trim $hit] \t] 3]
				regsub $k_val $data "$binomial \{$class\} \{$id\}:" data
			}
			set out [open $dir\_BS_tree.txt w]
			puts $out $data
			close $out
		}
	}
}

