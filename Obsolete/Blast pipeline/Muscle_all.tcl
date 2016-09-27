###
set direct ~/desktop/True_Pos
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################
#Create  a Group_MSA folder if one does not yet exist.
cd $direct

set z [glob -nocomplain -type d *]
if {[regexp Group_MSA $z] == 1} {
} else {
	file mkdir Group_MSA
}

###########################################################################
#Substitute all the gene identifiers with sequential numbers - as they are typically longer than 10 characters and the phylip format will truncate them beyond recognition.
set pata {>+?.+?\s+?}

set fastas [glob *faas.txt]
foreach fasta $fastas {
	set i 1
	openfile $fasta
	set new_glist ""
	set key_file ""
	regsub -all ">" $data "£>" data
	set glist [split $data £]
	regsub -all "{}" $glist {} glist
	foreach gene $glist {
		regexp $pata $gene hit
		regsub -all ">" $hit {} hit
		if {$i < 10} {
			regsub $hit $gene "k$i         " new_gene
		} elseif {$i >= 10 && $i < 100} {
			regsub $hit $gene "k$i        " new_gene
		} elseif {$i >= 100 && $i < 1000} {
			regsub $hit $gene "k$i       " new_gene
		} else {
			regsub $hit $gene "k$i      " new_gene
		}
		append new_glist "$new_gene"
		append key_file "k$i\t$hit\n"
		incr i
	}
	set out [open $fasta w]
	puts $out $new_glist
	close $out
	#Produce a key which matches all the identifiers with their respective numbers.
	set out2 [open $direct/Group_MSA/KEY_$fasta a]
	puts $out2 $key_file
	close $out2
}

###########################################################################
#Run the actual muscle on the now numbered genes in each family.

set fastas [glob *faas.txt]
foreach fasta $fastas {
	regsub "_fam_faas.txt" $fasta ".fas" groupno
	catch {exec muscle -in $fasta -phyiout MSA_$groupno}
	file rename MSA_$groupno\.fas $direct/Group_MSA/MSA_$groupno\.fas
}

###########################################################################


#set gene_fams "~/Desktop/Test_Pairwise_Blast/Input/Ortho"
#set msa "~/Desktop/Test_Pairwise_Blast/Input/MSA"
#file copy $gene_fams $msa

#cd ~/Desktop/Test_Pairwise_Blast/Input/MSA
#set gene_fams [glob *.faa]
#regsub seed.faa $gene_fams {} gene_fams
#set i 0

#exec ~/desktop/software/muscle3.8.31_i86darwin64 -version

#foreach gene $gene_fams {
#	puts "$i/[llength $gene_fams]"
#	catch {exec ~/desktop/software/muscle3.8.31_i86darwin64 -in $gene -physout $gene}
#	incr i
#}
