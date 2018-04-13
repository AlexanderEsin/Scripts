#!/usr/local/bin/tclsh
source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

# Directories
set master_dir			"/Users/aesin/Desktop/Geo_again"
set annot_dir			$master_dir/Functional_annotation
set annotInput_dir		$annot_dir/Output_annots

# Define input annotation files
set inputAnnot_files	[glob $annotInput_dir/*emapper.annotations]
set inputAnnot_files	[lsort -dictionary $inputAnnot_files]

# Define and open all prot DB
set all_db_file		$master_dir/All_prot_db_new
sqlite3 all_prot_db $all_db_file

# Add a column to database to hold orthologous group number
set table_info		[all_prot_db eval {PRAGMA table_info(t1)}]
if {[lsearch $table_info "COGcat"] == -1} {
	puts	"Adding COGcat column to database ..."
	all_prot_db eval {ALTER TABLE t1 ADD COLUMN COGcat text}
	puts -nonewline	"\rAdding COGcat column to database ... done \n"
} else {
	puts	"COGcat column already exists in database!"
}

set total_counter	1
set noCOG_counter	0

foreach annot_file $inputAnnot_files {

	# Track progress
	puts -nonewline		"\rAdding COG functional category to database: $total_counter // [llength $inputAnnot_files] ..."
	flush stdout

	# Get the ortholog group number
	set file_name	[file tail $annot_file]
	set orthGroup	[string range $file_name 0 [string first \. $file_name]-1]

	# Read in annot file
	set annotData	[split [string trim [openfile $annot_file]] \n]

	# Ignoring first two lines (timestamp and eggNog command call), 
	# get column header (line 3) and the actual data per protID
	set annotHeader	[lindex $annotData 2]
	set annotData	[lrange $annotData 3 end]

	# Remove last three lines from data (queries scanne, run time, rate)
	set annotData	[lrange $annotData 0 end-3]

	# Unhash the header line (first character)
	set annotHeader	[string range $annotHeader 1 end]

	# For each individual entry, extract the COG assignment
	set COGcat_list	{}
	foreach nogMap $annotData {
		set dataSplit	[split $nogMap \t]
		set COGcat		[lindex $dataSplit 10]

		# Sometimes each annotation has multiple COGs associated - split by ', '
		set COGcat		[wsplit $COGcat {, }]
		set COGcat_list	[concat $COGcat_list $COGcat]
	}

	# Check how many orthGroups have multiple COGs assigned
	set uniqueCOG	[lsort -unique $COGcat_list]
	# For groups with no COG assignment, retyrb
	if {[llength $uniqueCOG] == 0} {
		set uniqueCOG	"NA"
		incr noCOG_counter
	}
	# Join any multiples with pipe
	set orthGroupCOG	[join $uniqueCOG "|"]

	# Add functional category to database
	all_prot_db eval	{UPDATE t1 SET COGcat = $orthGroupCOG WHERE OrthGroup = $orthGroup}

	incr total_counter
}

puts "\nAdding functional categories to database - done"
puts "Groups with no assigned COG: $noCOG_counter / [llength $inputAnnot_files]"

all_prot_db close




