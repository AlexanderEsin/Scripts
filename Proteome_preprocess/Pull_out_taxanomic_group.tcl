#!/usr/local/bin/tclsh
###########################################################################
## Change Log ##
# 21 January 2016 : Replaced procs with sourced procs. Changed vars to lists. Made check for gz files #

###########################################################################
source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

set direct /users/aesin/desktop/Clean_Proteomes
set extract_groups {Bacillales}

###########################################################################

## Open databases ##
sqlite3 db1 $direct/all_prot_geo_ncbi_db
sqlite3 tax_ranks $direct/Taxonomic_ranks_db

cd $direct
set output_message {}

foreach group $extract_groups {

	## Make output folder for each taxanomic group ##
	file mkdir $direct/Extracted_groups/$group
	
	cd $direct/All_fastas_geo_mark

	## Check if any files are zipped, if they are throw an error and exit ##
	set zipped_files [glob -nocomplain *gz]
	if {[llength $zipped_files] > 0} {
		puts stderr "Fasta directory contains zipped files \([llength $zipped_files]\). Unzip all, then re-run script."
		exit 2
	}

	## Get a list of all the files ##
	set fastas [glob *faa]
	set i 1
	set x 0
	foreach file $fastas {
		db1 eval {SELECT binomial FROM t1 WHERE file_name = $file} {
		    set binomial $binomial
		}

		tax_ranks eval {SELECT ranks FROM t1 WHERE binomial = $binomial} {
			set org_entry $ranks
		}

		if {[string first "\t$group" $org_entry] > -1} {
			file copy -force $file $direct/Extracted_groups/$group
			incr x
		}

		puts "$i === $x"
		incr i
	}

	lappend output_message "$x files matching \"$group\" have been copied to $direct/Extracted_groups/$group"
}

puts [join $output_message \n]

db1 close
tax_ranks close
