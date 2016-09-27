###
set direct ~/desktop/True_pos
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################


#Open aligned phylip (interleaved) file

cd $direct/MSA
set filelist [glob *.fas]

foreach file $filelist {
	#Process original fasta file with organism names
	regsub ".fas" $file ".txt" fastafile
	cd $direct

	openfile $fastafile

	regsub -all > $data £> data
	set glist [split $data £]
	regsub -all "{}" $glist {} glist

	set headerlist ""

	foreach gene $glist {
		set firstn [string first \n $gene]
		set header [string trim [string range $gene 0 $firstn]]
		append headerlist "$header\n"
	}

	cd $direct/MSA
	openfile $file

	#Remove the numeric header

	set nfirst [string first \n $data]
	set phyheader [string range $data 0 $nfirst]
	set phyheader [string trim $phyheader]
	set align [string range $data $nfirst end]

	#Pull the first 'block' that contains all the gene names

	set firstnn [string first \n\n $align]
	set test [string range $align 0 $firstnn]
	set test [string trim $test]
	set test2 [split $test \n]
	foreach line $test2 {
		set name [lindex $line 0]; #Corresponds to the gene name on each line
		set pata {>+?};	#Patterns that will pull out the header from original fasta
		set patb {.+?\n+?}
		regexp $pata$name$patb $headerlist hit
		set patc {\[+?.+?\]+?}; #Pattern that identifies the organism name at the end of header
		regexp -nocase $patc $hit match
		regsub -all {\[} $match {} match
		regsub -all {\]} $match {} match
		regsub ";" $match {} match; #Trim the match
		regsub -all { } $match "_" match; #Necessary otherwise tree builder fucks up.
		regsub -all {\-} $match "_" match
		set match [string trim $match]
		regsub -all $name $align $match align
	}

	set output "$phyheader$align"

	set out [open $direct/MSA/MSA_$file a]
	puts $out $output
	close $out
}


