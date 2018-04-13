#!/usr/local/bin/tclsh
source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

# Directories
set master_dir			"/Users/aesin/Desktop/Geo_again"
set genome_dir			$master_dir/Genomes
set nuclGenomes_dir		$genome_dir/Genome_raw_fastas

# Define input annotation files
set nuclGenomes_files	[glob $nuclGenomes_dir/*fasta]

# Define the list of acc_ass coresponding to AG genomes
set taxid_accAss_trans	[split [string trim [openfile $genome_dir/Genome_lists/Acc_ass_taxid_table.tsv]] \n]

# Define and open all_prot database
set all_db_file			$master_dir/All_prot_db_new
sqlite3 all_prot_db 	$all_db_file

# # Add a new column to the database to hold the nucleotide sequence
set table_info		[all_prot_db eval {PRAGMA table_info(t1)}]
if {[lsearch $table_info "NuclSeq"] == -1} {
	puts	"Adding NuclSeq column to database ..."
	all_prot_db eval {ALTER TABLE t1 ADD COLUMN NuclSeq text}
	puts -nonewline	"\rAdding NuclSeq column to database ... done \n"
} else {
	puts	"NuclSeq column already exists in database!"
}

# Genome counters
set genomeCounter 	1
set genomeTotal		[llength $taxid_accAss_trans]

foreach genomeEntry $taxid_accAss_trans {
	set acc_ass		[lindex $genomeEntry 0]
	set taxid		[lindex $genomeEntry 1]

	# For each accession, find the file
	set nuclGenome_file	[lsearch -inline -glob $nuclGenomes_files "*$acc_ass*"]
	
	# Check the file exists
	if {[file exists $nuclGenome_file] != 1} {
		error "\nCannot find file for accession = $acc_ass"
	}

	# Get the genes
	set nuclGenes	[split_genes [string trim [openfile $nuclGenome_file]]]

	# Per-genome counters
	set processCounter 	1
	set totalCounter	[llength $nuclGenes]

	foreach nuclGene $nuclGenes {

		# Divivde fasta entry into header + seq
		set geneHead	[string trim [string range $nuclGene 0 [string first \n $nuclGene]]]
		set geneSeq		[string trim [string range $nuclGene [string first \n $nuclGene] end]]

		# Get header items
		set headElems	[wsplit $geneHead " | "]
		set protIDtrunc	[string range [lindex $headElems 0] 1 end]
		set locusTag	[lindex $headElems 1]

		# If the protID name exists, add the nuclSeq to the database 
		if {$protIDtrunc ne "missing_protein_id_qualifer"} {
			# Remove the space in the nucleotide sequence
			set geneSeq_flat	[string map {\n {}} $geneSeq]

			# Add asterisk outside of search expression because tcl is stupid
			set protIDtruncSearch	"$protIDtrunc*"

			# See how many entries we find for this gene - should be 1?
			set foundInDB		[all_prot_db eval {SELECT COUNT(*) FROM t1 WHERE protID GLOB $protIDtruncSearch AND locus == $locusTag AND acc_ass == $acc_ass}]
			if {$foundInDB != 1} {
				error "Finding too many or too few entries in DB. Acc_ass: $acc_ass. ProtID: $protIDtrunc. Locus: $locusTag"
			}

			# Add the flat (no-gap) nucleotide sequence to the database
			all_prot_db eval	{UPDATE t1 SET NuclSeq = $geneSeq_flat WHERE protID GLOB $protIDtruncSearch AND locus == $locusTag AND acc_ass == $acc_ass}
		}

		# Write "Gene Complete" tracker output to stdout
		puts -nonewline "\rProcessing genome $acc_ass \($genomeCounter // $genomeTotal\) -- Sequence completed: $processCounter // $totalCounter..."
		flush stdout
		incr processCounter

	}
	# Write "Genome Complete" tracker output to stdout
	puts -nonewline "\rProcessing genome $acc_ass \($genomeCounter // $genomeTotal\) -- Sequence completed: $processCounter // [expr $totalCounter - 1]...done\n"
	flush stdout
	incr genomeCounter
}

all_prot_db close


# ------------------------------------------------------------------------------------- #

## Here we confirm that the length of our nucleotide sequence matches expected from DB

# In practice, the only times where the sequence doesn't match is where there is an intron or predicted ribosomal slippage
# E.g. protID: "WP_087943500.1.398511_1274894 genome: GCF_000005825.2_ASM582v2" OR "NP_230169.1.243277_552383 genome: GCF_000006745.1_ASM674v1"

# set nuclLegth	[string length $geneSeq_noGap]
# set gene_coords	[all_prot_db eval {SELECT gene_start, gene_end from t1 where protID == $protID}]
# set geneStart	[lindex $gene_coords 0]
# set geneEnd		[lindex $gene_coords 1]

# set coordLength	[expr abs($geneStart - $geneEnd) + 1]
# if {$nuclLegth != $coordLength} {
# 	## If the values are unqual this could be due to predicted ribosomal slippage - 
# 	puts "\nNucleotide sequence length does not match expected. protID: $protID genome: $acc_ass"
# }

# ------------------------------------------------------------------------------------- #
