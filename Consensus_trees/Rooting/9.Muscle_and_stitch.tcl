#!/usr/local/bin/tclsh
# Muscle the combined 'core' genome
set eval -10
set ival 5
set org [lindex $argv 0]
set direct /users/aesin/desktop/Consensus_trees/Rooting/$org/$eval/

puts -nonewline "Enter your name of directory on server: "
flush stdout
set server_dir_name [gets stdin]

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

#This proc remakes the 80-wide fasta entries from a flat CDS sequence

proc linebreak {s {width 80}} {
   global res
   set res ""
   while {[string length $s]>$width} {
       set res "$res[string range $s 0 79]\n"
       set s [string range $s 80 end]
   }
   set res "$res$s"
}

proc split_genes {fasta} {
	global genes
	regsub -all {>} [string trim $fasta] {£>} fasta
	set fasta [split $fasta £]
	regsub -all "{}" $fasta {} genes
}

###########################################################################
## Create namelist (easy for core as all fams have same species) ##

cd $direct
file mkdir Stitched_core_genome_$ival
file mkdir Group_fastas_MSA_$ival
cd $direct/Group_fastas_key_ids_$ival

set fastas [glob *.faa]
set fastas [lsort -dictionary $fastas]

set species_names ""
openfile [lindex $fastas 0]
split_genes $data
foreach gene $genes {
	set hit [regexp -all -inline {\[+?.+?\]+?} $gene]
	set binomial [lindex $hit end]
	regsub -all {\[} $binomial {} binomial
	regsub -all {\]} $binomial {} binomial
	append species_names "$binomial\t"
}

set namelist [split [string trim $species_names] \t]

###########################################################################
## Make namelist into a key ##

set i 1
set key ""
foreach name $namelist {
	append key "k$i\t$name\n"
	incr i
}

set key [string trim $key]
set out [open Species_KEY.txt w]
puts $out $key
close $out

set key [split $key \n]

###########################################################################
## Make reference file to backtrack the genes to original fastas ##
set patb {>.+?}

# Look at a previous MSA_reference_file.tsv for layout - but basically a header with all the fasta numbers, then each organism gets a separate line with each gene in order #
set reference_file ""
append reference_file "Key_No\tName\t"

foreach fasta $fastas {
	regsub {.faa} $fasta {} g
	append reference_file "$g\t"
}

set reference_file [string trim $reference_file]

foreach name $key {
	set name_list [split $name \t]
	set binomial [lindex $name_list 1]
	set key_no [lindex $name_list 0]
	append reference_file "\n$name\t"
	foreach fasta $fastas {
		openfile $fasta
		regexp -line $patb$binomial\]+? $data comment_line
		regsub -all {>} $comment_line {} comment_line
		append reference_file "$comment_line\t" 
		set xx [regsub -all -line $patb$binomial\]+? $data ">$key_no" data]
		if {$xx > 1 || $xx == 0} {
			puts "ERROR"
			puts $fasta
			break
		}
		set out [open $fasta w]
		puts $out $data
		close $out
	}
}
set reference_file [string trim $reference_file]
set out [open $direct/Group_fastas_key_ids_$ival/MSA_reference_file.tsv w]
puts $out $reference_file
close $out

set fastas [glob *.faa]

set muscle_counter 1

foreach file $fastas {
	regsub {.faa} $file {.fas} fas
	catch {exec muscle -in $direct/Group_fastas_key_ids_$ival/$file -out $direct/Group_fastas_MSA_$ival/$fas}
	puts "MUSCLING: $muscle_counter / [llength $fastas]"
	incr muscle_counter
}

file copy -force Species_KEY.txt $direct/Group_fastas_MSA_$ival/
file copy -force Species_KEY.txt $direct/Stitched_core_genome_$ival/
file copy -force MSA_reference_file.tsv $direct/Group_fastas_MSA_$ival/
file copy -force MSA_reference_file.tsv $direct/Stitched_core_genome_$ival/

###########################################################################
## Stich the genome via id number ##
cd $direct/Group_fastas_MSA_$ival

# The newline here is essential to pick up k1 and not k17 #
set paty {\n+?.+?>+?}

set fams [glob *.fas]
set fams [lsort -dictionary $fams]

openfile [lindex $fams 0]
set id_numbers [regexp -all -line -inline {>+.+} $data]

set check_length ""
set MSA_master ""
foreach id_number $id_numbers {
	set master_CDS ""
	append MSA_master "$id_number\n"
	foreach fam $fams {
		openfile $fam
		set data "[string trim $data]\n>\n"
		if {[regexp $id_number$paty $data entry] == 0} {
			puts "ERROR: $fam $id_number"
		}
		set entry [string trim [string range $entry 0 end-1]]
		set header [string first \n $entry]
		set CDS [string trim [string range $entry [string first \n $entry] end]]
		regsub -all "\n" $CDS {} CDS
		append master_CDS $CDS
	}
	append check_length "$id_number === [string length $master_CDS]\n"
	linebreak $master_CDS
	append MSA_master "$res\n"
}
set MSA_master [string trim $MSA_master]

set out [open $direct/Stitched_core_genome_$ival/MSA_master.faa w]
puts $out $MSA_master
close $out

cd $direct/Stitched_core_genome_$ival

catch {exec seqret MSA_master.faa MSA_master.fas -osformat phylip}

puts "DONE"
puts [string trim $check_length]

###########################################################################
exec scp -r $direct/Stitched_core_genome_$ival ade110@ax3.hpc.ic.ac.uk:/scratch/ade110/Consensus_trees/$server_dir_name >&@stdout
###########################################################################
