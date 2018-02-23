#!/usr/local/bin/tclsh

## Making a database of all the Anogeo proteins and information about them
## Information to include: 
##		proteins from the clean proteomes and their names and loci
##		duplicate protein from the 'duplicate' proteomes
##		plasmid status + genome length (from the All_plasmid_tags file: 
##			see Scripts/Anogeo/Genome_processing/Plasmid_tags_&_genome_lengths.tcl)
##		species binomial and strain information from assembly_summary_refseq_13117.txt


source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

set direct				/users/aesin/Desktop/Geo_again

## Directories
set prot_clean_dir		$direct/Proteomes/Proteome_clean_fastas
set prot_dupl_dir		$direct/Proteomes/Proteome_dup_fastas
set plasmid_tag_file	$direct/Genomes/All_plasmid_tags.tsv
set genome_length_file	$direct/Genomes/Genome_lengths.tsv
set refseq_file			$direct/Genomes/assembly_summary_refseq_131117.txt
set anogeo_acc_ass_file	$direct/Genomes/Genome_lists/AG_acc_ass_names.txt

set anogeo_prot_db		$direct/All_prot_db

## Open the database and define settings
sqlite3 db1 $anogeo_prot_db
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
		plasmid		text,
		product		text,
		sequence	text,
		acc_ass		text,
		binomial	text,
		strain		text,
		taxid		integer,
		genome_l	integer,
		is_ag		integer
	)
}

## Process the global datasets:
##		plasmid tags / genome lengths / refseq annotation / Anogeo acc_ass
set plasmid_tag_list	[split [string trim [openfile $plasmid_tag_file]] \n]
set genome_length_list	[split [string trim [openfile $genome_length_file]] \n]
# String trim removes trailing tabs (representing the last few empty columns of the last entry)
set refseq_data_list	[split [openfile $refseq_file] \n]
set anogeo_acc_ass_list	[split [string trim [openfile $anogeo_acc_ass_file]] \n]

## Foreach clean proteome, collect the necessary information
set clean_proteomes		[glob $prot_clean_dir/*fasta]
set databased_counter 1

foreach clean_proteome $clean_proteomes {

	## Derive the accession_assembly from the file name
	set file_name	[file tail $clean_proteome]
	set acc_ass		[string range $file_name 0 end-15]

	## Check whether this proteome is Anogeo
	if {[lsearch $anogeo_acc_ass_list $acc_ass] != -1} {
		set is_ag	1
	} else {
		set is_ag	0
	}

	puts "Working on $acc_ass -- $databased_counter // [llength $clean_proteomes] ..."

	## Use accession_assembly to find genome length
	set genome_length	[lindex [split [lsearch -inline $genome_length_list "$acc_ass\t*"] \t] 1]

	## Use accession_assembly to find the taxid / binomial / strain IDs
	set refseq_entry	[lsearch -all -inline $refseq_data_list "*$acc_ass\t*"]
	# Check we only have a single hit
	if {[llength $refseq_entry] > 1 || [llength $refseq_entry] == 0} {
		error "Expecting to have one unique hit"

	}

	set refseq_split	[split $refseq_entry \t]
	set taxid			[lindex $refseq_split 5]
	set binomial		[lindex $refseq_split 7]
	set strain			[lindex $refseq_split 8]
	# Remove the "strain=" from the strain format
	set strain_clean	[string map {strain= {}} $strain]


	## Now process each protein in turn from the clean proteome
	## and the duplicate proteome (if it exists)
	## If there are dup proteins, add them onto the end of the clean
	## protein list
	set proteome_data	[string trim [openfile $clean_proteome]]
	set protein_list	[split_genes $proteome_data]

	if {[file exists $prot_dupl_dir/$file_name] == 1} {
		set dup_proteome_data	[string trim [openfile $prot_dupl_dir/$file_name]]
		set dup_protein_list	[split_genes $dup_proteome_data]
		set protein_list		[concat $protein_list $dup_protein_list]
	}

	foreach protein $protein_list {
		## Have to take the range from the first (not 0th) character 
		## to omit the fasta ">"
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

		## Not all the plasmid tags are found in the proteome
		## because some plasmid tags refer to genes that don't
		## have an associated protein sequence. E.g. see
		## GK_RS00045 in GCF_000009785.1_ASM978v1
		if {$is_ag == 1} {
			set plasmid		"F"
			set is_plasmid		[lsearch $plasmid_tag_list $locus]
			if {$is_plasmid != -1} {
				set plasmid "T"
			}
		} else {
			set plasmid		"NA"
		}
		

		## Get the AA sequence
		set sequence	[string trim [string range $protein [string first \n $protein] end]]

		## Put all our data into the database
		db1 eval {
			insert into t1 values(
				$protID,
				$locus,
				$gene_start,
				$gene_end,
				$strand,
				$plasmid,
				$product,
				$sequence,
				$acc_ass,
				$binomial,
				$strain_clean,
				$taxid,
				$genome_length,
				$is_ag
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
puts "\nAll database entry complete: see $anogeo_prot_db"














