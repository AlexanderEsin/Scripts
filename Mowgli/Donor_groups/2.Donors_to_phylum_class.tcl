#!/usr/local/bin/tclsh
source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

proc get_rank {species taxonomy_table {rank Phylum}} {
	global rank_out

	if {$rank eq "Phylum"} {
		set index end-3
	} elseif {$rank eq "Class"} {
		set index end-4
	} elseif {$rank eq "Order"} {
		set index end-5
	}

	set ranks [lsearch -glob -inline $taxonomy_table $species\t*]
	if {[string length $ranks] == 0} {
		# puts stderr "Entry for $species not found in taxonomy table"
		return
	} else {
		set rank_out [lindex [split $ranks \t] $index]
	}

	return $rank_out
}

proc get_most_abundant_rank {present_ranks} {
	global most_abundant_ranks

	set rank_abundance_tbl [lsort -integer -index 1 -decr [lcount $present_ranks]]
	set most_abundant_ranks [list [lindex [lindex $rank_abundance_tbl 0] 0]]

	if {[llength $rank_abundance_tbl] > 1} {
		set index 0
		## As long as the next category had an equivalent number of hits, add it to the list as well. E.g. Firmicutes Actinobacteria ##
		while 1 {
			if {[lindex [lindex $rank_abundance_tbl $index] 1] == [lindex [lindex $rank_abundance_tbl [expr $index + 1]] 1]} {
				lappend most_abundant_ranks [lindex [lindex $rank_abundance_tbl [expr $index + 1]] 0]
				incr index
			} else {
				break
			}
		}
	}

	return $most_abundant_ranks
}

###########################################################################
set direct /users/aesin/desktop

set taxa_input_dir $direct/Geo_analysis/HGT_Donors_Basic/Long_distance_HGT/Donor_Taxa
set tree_directory $direct/Mowgli/Final_BS

## Make output directories ##
set phyla_output_dir $taxa_input_dir/../Phylum
file mkdir $phyla_output_dir

set class_output_dir $taxa_input_dir/../Class
file mkdir $class_output_dir


## Pick penalty and minimum sequences in group for count to be included##
set penalty 4
set min_taxa 0

###########################################################################
## Load in the key that translates the species tree tip labels into the original taxon names. The first column is the original name, the second column is the name used in the species tree ##
set tip_fix_key_list [split [string trim [openfile $direct/Mowgli/Species_tree/Tip_fix_key.tsv]] \n]

###########################################################################
## Load in the taxonomy table so that we can identify all the extant donor species to phylum level ##
set taxonomy_table_list [split [string trim [openfile /users/aesin/desktop/Clean_Proteomes/Taxonomy_all_2.tsv]] \n]

###########################################################################

cd $taxa_input_dir
set donor_input_files [lsort -dictionary [glob *.txt]]

foreach input_file $donor_input_files {

	set penalty [string range $input_file 1 [string first "_" $input_file]-1]
	set min_taxa [string range $input_file [string last "_" $input_file]+1 end-4]

	puts stdout "\n######################################\n"
	puts stdout "For penalty $penalty with a minimum sequence cutoff of $min_taxa:"

	## Open the table with potential donor taxa. There is no header ##
	set group_donors [split [string trim [openfile $input_file]] \n]

	# set donor_set [lsearch -glob -inline $group_donors 4075*]

	set group_donor_phylum_table {}
	set group_donor_class_table {}
	set no_extant_taxa_counter 0
	set no_extant_taxa_l {}
	set id_to_class 0

	foreach donor_set $group_donors {
		set group [lindex [split $donor_set \t] 0]
		set donor_taxa [lrange [split [string trim $donor_set \t]] 1 end]

		set present_donors {}
		set present_phyla {}
		set present_classes {}

		## Open the best ML tree file from which the reconciliation was made to check whether the taxa are actually present ##
		set tree_data [string trim [openfile $tree_directory/$group\_tree.txt]]

		## Cleave off the numeric identifier appended by Mowgli, translate into the original taxon name, and if actually extant in the group, add it to a refined list (present_donors) ##
		foreach taxon $donor_taxa {
			set original_name ""
			set cleaved_name [string range $taxon 0 [string last \_ $taxon]-1]
			set original_name_index [lsearch -glob $tip_fix_key_list *\t$cleaved_name]
			set original_name [lindex $tip_fix_key_list $original_name_index 0]

			## If we can't get the original name. Error ##
			if {$original_name eq ""} {
				puts stderr "Error: Did not find original name.\nTaxon: $taxon\nGroup: $group\n"
				exit 2
			}

			## Check if the taxon is actually present in the group ##
			set search_string [escape_special $original_name]
			set present [regexp $search_string $tree_data]
			
			if {$present == 1} {
				lappend present_donors $original_name
			}
		}

		########################

		if {[llength $present_donors] == 0} {
			#puts stderr "For group : $group there seem to be NO extant donor taxa at penalty $penalty"
			incr no_extant_taxa_counter
			lappend no_extant_taxa_l $group
			continue
		}

		########################
		## For each species in the extant donor list, identify the phylum ##
		foreach donor_species $present_donors {
			set present_phylum [get_rank $donor_species $taxonomy_table_list]
			if {[string length $present_phylum] != 0} {
				lappend present_phyla $present_phylum
			}
		}
		
		## Identify the most abundant phyla. If two or more are equally represented, add both ##
		set most_abundant_phyla [get_most_abundant_rank $present_phyla]

		## Add the most abundant phylum/a for each group to the output table ##
		foreach phylum $most_abundant_phyla {
			lappend group_donor_phylum_table "$group\t$phylum"
		}

		########################
		## For each species, identify the class ##
		foreach donor_species $present_donors {
			set present_class [get_rank $donor_species $taxonomy_table_list Class]
			if {[string length $present_class] != 0} {
				lappend present_classes $present_class
			}
		}

		## Identify the most abundant phyla. If two or more are equally represented, add both ##
		set most_abundant_classes [get_most_abundant_rank $present_classes]

		## Add the most abundant class(es) for each group to the output table ##
		foreach class $most_abundant_classes {
			lappend group_donor_class_table "$group\t$class"
		}

		########################
	}

	## Write out the phyla ##
	set out [open $phyla_output_dir/T$penalty\_donor_phylum_min_$min_taxa\.tsv w]
	puts $out [join $group_donor_phylum_table \n]
	close $out

	## Write out the classes ##
	set out [open $class_output_dir/T$penalty\_donor_class_min_$min_taxa\.tsv w]
	puts $out [join $group_donor_class_table \n]
	close $out

	
	puts stdout "\tNumber of groups with OG donors: [llength $group_donors]"
	puts stdout "\tNumber of groups in which none of the edge descendants are extant: $no_extant_taxa_counter"
	puts stdout "\t\tGroups with no extant taxa: [join $no_extant_taxa_l " "]"
	puts stdout "\tNumber of groups further identified to class: $id_to_class"
	puts stdout "\tNumber of phyla included in output: [llength $group_donor_phylum_table]"
	puts stdout "\tNumber of classes included in output: [llength $group_donor_class_table]"

}



