#!/usr/local/bin/tclsh
######
source ~/Dropbox/Scripts/General_utils.tcl
set direct /users/aesin/Desktop/Cleaning_prots
set org Bacteria_Curated
set annotation "Bact"
#####

###########################################################################

#Create output dir in marked_proteomes

file mkdir $direct/$org/Temp

#Open the input dir

cd $direct/$org/Clean_faa
set faas [glob *.faa]

#For each file, mark each gene header with either Arch, Bact, or Euk between parantheses.
#REMEMBER to change the DESIG below

foreach faa $faas {
	openfile $direct/$org/Clean_faa/$faa
	set glist [split_genes $data]
	set new_proteome {}
	foreach gene $glist {
		# What is this regsub for? Seems to remove all (bracketed terms) followed by a space: "(x) ". Is it to remove any previous annotation?
		regsub -all {\(+?.+?\)+?\s+?} [string trim $gene] {} gene

		regexp {>.+?\s+?} $gene prot_id
		regsub $prot_id $gene "$prot_id\($annotation\) " gene
		lappend new_proteome $gene
	}
	set out [open $direct/$org/Temp/$faa w]
	puts $out [join $new_proteome \n]
	close $out
}

set no_marked_fasta_files [llength [glob $direct/$org/Temp/*faa]]

if {[llength $faas] == $no_marked_fasta_files} {
	puts "All files successfully annotated"
} else {
	puts "The number of input =/= number of output"
}



###########################################################################