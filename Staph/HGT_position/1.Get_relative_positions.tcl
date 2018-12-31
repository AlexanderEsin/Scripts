#!/usr/local/bin/tclsh

source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

proc GetRelativePosition {gene_pos oriPos_data hostGenomeLen} {
	set oriPos_strand	[lindex [split $oriPos_data \t] end]

	if {$oriPos_strand eq "Forward"} {
		# Position of the origin start
		set oriStart_pos	[lindex [split $oriPos_data \t] 1]

		# Get the position relative to the origin
		set relGene_pos		[expr $gene_pos - $oriStart_pos]

		# If the relative location is negative then it lies upstream of the origin
		if {$relGene_pos < 0} {
			set relGene_pos	[expr $relGene_pos + $hostGenomeLen]
		}

	} elseif {$oriPos_strand eq "Reverse"} {
		# Position of the origin start
		set oriStart_pos	[lindex [split $oriPos_data \t] 2]

		# Get the relative position #
		set relGene_pos		[expr $gene_pos - $oriStart_pos]

		# If it's negative it lies downstream of the origin. Take the absolute value and leave as is #
		# If it's positive it lies upstream. We subtract it's position from the genome size. This will result in chromStart being > chromEnd #
		if {$relGene_pos < 0} {
			set relGene_pos	[expr abs($relGene_pos)]
		} else {
			set relGene_pos [expr $hostGenomeLen - $relGene_pos]
		}				
	}

	## Get the fractional position
	set fractGene_pos		[tcl::mathfunc::roundto [expr double($relGene_pos) / double($hostGenomeLen)] 8]
	return $fractGene_pos
}

# Paths and output folders
set master_dir		/Users/aesin/Desktop/Staph
set mowgli_dir		$master_dir/Mowgli/Mowgli_output
set genome_dir		$master_dir/Core_genomes/Genome_lists
set posit_dir		$master_dir/HGT_position/Position_data

# Find and open the DB
set all_db_file		$master_dir/ForStaph_prot_db
sqlite3 allProt_db	$all_db_file

# Directory to find tipID <-> protID keys
set tipKey_dir		$master_dir/Mowgli/GeneTree_input

# Define and make output directory
set cleaned_dir		$mowgli_dir/Cleaned_events
set hgtClean_dir	$cleaned_dir/HGT_events
set verClean_dir	$cleaned_dir/Ver_events


## Positions of starting genes (dnaA and polB III)
# Input file (for each acc_ass the locus tags for dnaA and polB III)
set startLocTag_file	$genome_dir/Core_origin_gene_list.tsv
set startLocTag_data	[split [string trim [openfile $startLocTag_file]] \n]

# Read in Acc_Ass <-> Taxid translation table 
set accAssTrans_file	$genome_dir/core_AccAssTaxid_table.tsv
set accAssTrans_data	[split [string trim [openfile $accAssTrans_file]] \n]

# Output list
set oriPos_list			{}

# Remove header from the dnaA/polB locus list
set startLocTag_data	[lrange $startLocTag_data 1 end]
foreach entry $startLocTag_data {
	set acc_ass			[lindex [split $entry \t] 0]
	set dnaA_locTag		[lindex [split $entry \t] 1]

	## Translate the acc_ass to taxid
	set accAss_entry	[lsearch -all -inline $accAssTrans_data "$acc_ass\t*"]
	# Check there is only 1 hit in translation table
	if {[llength $accAss_entry] == 1} {set accAss_entry [lindex $accAss_entry 0]} else {error "Not finding exactly 1 entry for acc_ass: $acc_ass"}
	# Get the taxid
	set taxid			[lindex [split $accAss_entry \t] 1]

	set dbEntry			[allProt_db eval {SELECT taxid, gene_start, gene_end, strand FROM t1 WHERE taxid = $taxid AND locus = $dnaA_locTag}]
	set dbEntry_list	[join $dbEntry \t]
	lappend oriPos_list	$dbEntry_list
}

# Write out the start dnaA start position data
set oriPos_listHeader [join [list Taxid geneStart geneEnd Strand] \t]
set out [open $genome_dir/AG_dnaA_startEnd_positions.tsv w]
puts $out [join [concat [list $oriPos_listHeader] $oriPos_list] \n]
close $out

## Core taxids
set core_taxid_file		$master_dir/Core_genomes/Genome_lists/coreToKeep.tsv
set core_taxid_tbl		[lrange [split [string trim [openfile $core_taxid_file]] \n] 1 end]
set core_taxids		{}
foreach row $core_taxid_tbl {
	set taxid	[lindex [split $row \t] 2]
	lappend core_taxids $taxid
}

## Penalty and HGT-type lists
set penalty_list	[list 3 4]
set hgtType_list	[list "lHGT" "sHGT" "Ver" "All"]
# set hgtType_list	[list "All"]
set allDone_flag	FALSE

foreach hgtType $hgtType_list {

	foreach penalty $penalty_list {

		puts "Working on $hgtType penalty = $penalty..."

		# Prepare the output list
		set penaltyOut_head		[list "orthGroup" "protID" "locusTag" "taxid" "binomial" "COGcat" "geneStart" "geneEnd" "strand" "relGeneStart" "relGeneEnd" "relStrand"]

		if {$hgtType eq "lHGT" || $hgtType eq "sHGT"} {
			set inputEvents_file	$hgtClean_dir/T$penalty\_full_$hgtType\_events.tsv
			set inputEvents_data	[split [string trim [openfile $inputEvents_file]] \n]
			# Trim the header out < Group	Donor_nodes	Receptor_nodes	Tips >
			set inputEvents_data	[lrange $inputEvents_data 1 end]
			# Add extra column headers for the HGT data
			set penaltyOut_head		[concat $penaltyOut_head [list "donorEdge" "recepEdge" "eventIndex"]]
		} elseif {$hgtType eq "Ver"} {
			set inputEvents_file	[glob $verClean_dir/$hgtType\_const_t$penalty*.tsv]
			set inputEvents_data	[split [string trim [openfile $inputEvents_file]] \n]
		} elseif {$hgtType eq "All"} {
			if {$allDone_flag == TRUE} {
				puts "All dataset is already complete!"
				continue
			}
			
			set inputEvents_data {}
			foreach core_taxid $core_taxids {
				set input_events	[allProt_db eval {SELECT protID FROM t1 WHERE taxid = $core_taxid AND plasmid = "F"}]
				set inputEvents_data [concat $inputEvents_data $input_events]
			}
		}


		set penaltyOut_list		[list [join $penaltyOut_head \t]]
		set totalTipIDs			0
		set counter				1


		foreach event $inputEvents_data {

			puts -nonewline "\rProcessing event $counter // [llength $inputEvents_data]..."
			flush stdout

			if {$hgtType eq "lHGT" || $hgtType eq "sHGT"} {
				set eventSplit		[split $event \t]
				set eventIndex		[lindex $eventSplit 0]
				set orthGroup		[lindex $eventSplit 1]
				set donorEdge		[lindex $eventSplit 2]
				set recepEdge		[lindex $eventSplit 3]
				set tipIDs			[lindex $eventSplit 4]

				# Open the key file to translate tipID to protID
				set key_file	$tipKey_dir/$orthGroup/$orthGroup\_KEY_tips.txt
				set key_data	[lrange [split [string trim [openfile $key_file]] \n] 1 end]; # Remove header as well
			} elseif {$hgtType eq "Ver"} {
				
				set tipIDs {}
				foreach core_taxid $core_taxids {
					set tax_tipIDs	[allProt_db eval {SELECT protID FROM t1 WHERE OrthGroup = $event AND taxid = $core_taxid}]
					set tipIDs [concat $tipIDs $tax_tipIDs]
				}
			} else {
				set tipIDs		$event
			}
				
			set totalTipIDs		[expr $totalTipIDs + [llength $tipIDs]]


			foreach tipID $tipIDs {

				if {$hgtType eq "lHGT" || $hgtType eq "sHGT"} {
					set key_entry		[lsearch -all -glob -inline $key_data "*\t$tipID"]
					if {[llength $key_entry] == 1} {set key_entry [lindex $key_entry 0]} else {error "Finding too many key entries for $tipID in $orthGroup at $penalty at $hgtType"}
					set protID			[lindex [split $key_entry \t] 0]
				} else {
					set protID			$tipID
				}
				
				set prot_data		[allProt_db eval {SELECT OrthGroup, protID, locus, taxid, binomial, COGcat, gene_start, gene_end, strand, genome_l FROM t1 where protID = $protID}]

				set taxid			[lindex $prot_data 3]
				set hostGenomeLen	[lindex $prot_data end]
				set strand			[lindex $prot_data end-1]

				## Get the gene start and end positions (accounting for strand)
				if {$strand eq "Forward"} {
					set geneEnd			[lindex $prot_data end-2]
					set geneStart		[lindex $prot_data end-3]
				} elseif {$strand eq "Reverse"} {
					set geneEnd			[lindex $prot_data end-3]
					set geneStart		[lindex $prot_data end-2]
				} else {
					error "Strand must be either \'Forward\' or \'Reverse\'"
				}

				## Calculate the gene positions relative to the origin start
				set oriPos_data		[lsearch -inline $oriPos_list "$taxid\t*"]
				set relStart		[GetRelativePosition $geneStart $oriPos_data $hostGenomeLen]
				set relEnd			[GetRelativePosition $geneEnd $oriPos_data $hostGenomeLen]

				## Relative strand identity ('same' as the Ori or 'diff' to the ori)
				set oriStrand		[lindex [split $oriPos_data \t] end]
				if {$strand eq $oriStrand} {
					set relStrand	"same"
				} else {
					set relStrand	"diff"
				}

				## We don't need to output the genome length
				set prot_data		[lrange $prot_data 0 end-1]

				## Define the columns we want to add
				set colsToAdd		[list $relStart $relEnd $relStrand]
				# Add extra columns for the HGT sets
				if {$hgtType eq "lHGT" || $hgtType eq "sHGT"} {
					set colsToAdd	[concat $colsToAdd [list $donorEdge $recepEdge $eventIndex]]
				}

				## Add the columns
				set prot_data		[concat $prot_data $colsToAdd]

				## Add output to glo l list
				lappend penaltyOut_list	[join $prot_data \t]
			}
			incr counter
		}

		puts "\nProcessed a total of $totalTipIDs protein IDs\n"

		## Define and make output directory
		set output_dir			$posit_dir/$hgtType\_input
		file mkdir				$output_dir

		## Write out data
		if {$hgtType eq "All"} {
			set output_file		$output_dir/$hgtType\_positionData.tsv
			set allDone_flag	TRUE
		} else {
			set output_file		$output_dir/T$penalty\_$hgtType\_positionData.tsv
		}
		
		set output_chan		[open $output_file w]
		puts $output_chan	[join $penaltyOut_list \n]
		close $output_chan
	}
}


allProt_db close
