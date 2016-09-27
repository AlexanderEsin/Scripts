###

source ~/Dropbox/Scripts/General_utils.tcl
package require sqlite3

###

set direct /users/aesin/desktop/Clean_Proteomes
set db_name all_prot_geo_ncbi_db

###########################################################################

## This was added at a later stage to make a new column called NCBI_ID - this allows every protein ID to be linked up to the correct (already adjusted for renaming of some taxa between proteome/genome downloads) genome/gbk/feature table file ##
set org_ncbi_transl_tbl [split [string trim [openfile $direct/All_fastas_ncbi_id_binomial_translation.tsv]] \n]

cd $direct/All_fastas_geo_mark

set pata {>+?.+?\]+?\n+?}
set patb {>+?.+?\)+?}
set patc {\(+?.+?\)+?}
set patd {>+?.+?\s+?}
set pate {\)+?.+\[+?}
set patf {\[+?.+?\]+?}
set patx {>+?.+?(MULTISPECIES)+?.+?\n+?}

# PRAGMA increases speed - a lot #
sqlite3 db1 $db_name
db1 eval {PRAGMA main.page_size = 4096}
db1 eval {PRAGMA main.cache_size=10000}
db1 eval {PRAGMA main.locking_mode=EXCLUSIVE}
db1 eval {PRAGMA main.synchronous=NORMAL}
db1 eval {PRAGMA main.journal_mode=WAL}
db1 eval {create table t1(id text, class text, gene blob, file_name glob, binomial text, ncbi_id text)}

set fastas [glob *.faa]
set fasta_counter 1

set done_binomials {}

foreach fasta $fastas {
	set gene_counter 1
	openfile $fasta
	##################
	### Prepare the phylip-friendly binomial names in preparation for subbing in for MULTISPCIES entries ###
	# For each proteome - list the genes in reverse (most MULTISPECIES and alient tags are at the top), and find a gene with the proper binomial #
	regsub -all "_" $fasta " " new_name
	set new_name [join [lrange [split $new_name " "] 0 1]]

	set ncbi_id ""
	set file_name_no_suff [string range $fasta 0 [string last \. $fasta]-1]
	set ncbi_id [lindex [split [lsearch -inline -glob $org_ncbi_transl_tbl *\t$file_name_no_suff] \t] 0]
	if {$ncbi_id eq ""} {
		puts $fasta
		exit
	}

	set genes [split_genes_fast $data]
	set rev_genes [lreverse $genes]
	set binomial ""
	foreach gene $rev_genes {
		if {[string first "MULTISPECIES" $gene] == -1} {
			set test_binomial [string range [string range $gene 0 [string first \n $gene]-1] [string last \[ $gene]+1 end-1]
			if {$test_binomial == $new_name} {
				set binomial $test_binomial
				break
			}
		}
	}

	if {[string length $binomial] == 0} {
		foreach gene $rev_genes {
			if {[string first "MULTISPECIES" $gene] == -1} {
				set test_binomial [string range [string range $gene 0 [string first \n $gene]-1] [string last \[ $gene]+1 end-1]
				if {[string first $binomial $done_binomials] == -1} {
					set binomial $test_binomial
					break
				}
			}
		}
	}
	if {[string length $binomial] == 0} {
		puts "BINOMIAL NOT FOUND == $file"
		error
	}
	lappend done_binomials $binomial
	# Clean up extra characters that can be present in the binomial names, this is also repeated later when assigning fasta to groups after MCL #
	regsub -all {\*} $binomial {} match
	regsub -all {\(} $match {} match
	regsub -all {\)} $match {} match
	regsub -all {\:} $match {} match
	regsub -all {\#} $match {} match
	regsub -all {\=} $match {} match
	regsub -all {\/} $match {} match
	regsub -all {'} $match {} match
	regsub -all " " $match {_} match
	regsub -all {\.} $match {} match
	regsub -all {\,} $match {} match
	regsub -all {\+} $match {} match
	regsub -all -- {\-} $match {_} binomial
	##################
	### Prepare the id, class, and gene entries ###
	# Split the fasta entry into a gene list #
	# For each gene pull out the white-space trimmed gene, id (>XXXX), class (Arch), the name of the fasta file, and the clean binomial name #
	foreach gene $genes {
		set gene [string trim $gene]
		set header_n [string first \n $gene]
		set comment_line [split [string trim [string range $gene 0 $header_n]] " "]
		set id [string range [lindex $comment_line 0] 1 end]
		set class [lindex $comment_line 1]
				##################
		# Insert the above values into a row of the database, table: t1 #
		db1 eval {insert into t1 values($id,$class,$gene,$fasta,$binomial, $ncbi_id)}
		incr gene_counter
	}
	puts "$fasta ==== $ncbi_id ==== $fasta_counter / [llength $fastas]"
	incr fasta_counter
}
# Create an index for more rapid searching and retrieval #
db1 eval {CREATE INDEX id_index ON t1 (id)}
db1 eval {CREATE INDEX file_index ON t1 (file_name)}

db1 close

# Copy the new db into the main directory #
file rename -force $db_name $direct/$db_name

puts "DONE"