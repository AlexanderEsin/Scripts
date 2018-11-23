#!/usr/local/bin/tclsh

source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

# Directories
set direct		/Users/aesin/Desktop/HTa_tonio
set prot_dir	$direct/Proteomes
set clean_out	$prot_dir/Proteome_clean
set genome_list	$direct/Genomes/All_genome_list.tsv

set db_dir		$direct/HTa_prot_db

# Open the database and define settings
sqlite3 db1 $db_dir
db1 eval {PRAGMA main.page_size = 4096}
db1 eval {PRAGMA main.cache_size=10000}
db1 eval {PRAGMA main.locking_mode=EXCLUSIVE}
db1 eval {PRAGMA main.synchronous=NORMAL}
db1 eval {PRAGMA main.journal_mode=WAL}

## Define the columns we want
db1 eval {
	create table t1(
		protID		text,
		locus		text,
		gene_start	integer,
		gene_end	integer,
		strand		text,
		product		text,
		sequence	text,
		acc_ass		text,
		binomial	text,
		strain		text,
		taxid		integer
	)
}


# String trim removes trailing tabs (representing the last few empty columns of the last entry)
set refseq_data_list	[split [openfile $genome_list] \n]

# Foreach clean proteome, collect the necessary information
set clean_proteomes		[glob $clean_out/*fasta]
set databased_counter 1

foreach clean_proteome $clean_proteomes {

	# Derive the accession_assembly from the file name
	set file_name	[file tail $clean_proteome]
	set acc_ass		[string range $file_name 0 end-15]

	puts "Working on $acc_ass -- $databased_counter // [llength $clean_proteomes] ..."


	# Use accession_assembly to find the taxid / binomial / strain IDs
	set refseq_entry	[lsearch -all -inline $refseq_data_list "*$acc_ass\t*"]
	# Check we only have a single hit
	if {[llength $refseq_entry] > 1 || [llength $refseq_entry] == 0} {
		error "Expecting to have one unique hit for $acc_ass"

	}

	set refseq_split	[split $refseq_entry \t]
	set taxid			[lindex $refseq_split 5]
	set binomial		[lindex $refseq_split 7]
	set strain			[lindex $refseq_split 8]
	# Remove the "strain=" from the strain format
	set strain_clean	[string map {strain= {}} $strain]


	# Now process each protein in turn from the clean proteome
	set proteome_data	[string trim [openfile $clean_proteome]]
	set protein_list	[split_genes $proteome_data]

	foreach protein $protein_list {
		# Have to take the range from the first (not 0th) character to omit the fasta ">"
		set header		[split [string trim [string range $protein 1 [string first \n $protein]]] "|"]
		set protID		[string trim [lindex $header 0]]
		set locus		[string trim [lindex $header 1]]

		set product		[string trim [lindex $header 3]]
		set location	[string trim [lindex $header 4]]

		set loc_split	[split $location " "]
		set strand		[lindex $loc_split 1]

		set boundaries	[split [lindex $loc_split 0] \:]
		set gene_start	[lindex $boundaries 0]
		set gene_end	[lindex $boundaries 1]

		# Get the AA sequence
		set sequence	[string trim [string range $protein [string first \n $protein] end]]

		# Put all our data into the database
		db1 eval {
			insert into t1 values(
				$protID,
				$locus,
				$gene_start,
				$gene_end,
				$strand,
				$product,
				$sequence,
				$acc_ass,
				$binomial,
				$strain_clean,
				$taxid
			)
		}
	}

	incr databased_counter
}


puts "\nMaking index on the protID and taxid columns ..."
## Make an index on the protID
db1 eval {create index protID_index on t1 (protID)}
db1 eval {create index taxid_index on t1 (taxid)}
## Close the database
db1 close
puts "\nAll database entry complete: see $db_dir"


	