###
set direct ~/desktop/True_pos/

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct/Group_MSA/BT
set fams [glob *.fas]

foreach fam $fams {
	regsub "_fam_faas.fas" $fam {} fname
	exec raxml -f a -s $fam -n $fname.txt -m PROTCATDAYHOFF -p 1234 -x 10001 -N 10 -T 7
}

set pata {k+?[0-9]+?:+?}
set patb {\t+?.+?\n+?}
set patc {>+?}
set patd {.+?\n+?}
set pate {\[+?.+?\]+?}

set besttrees [glob *bestTree*]
foreach tree $besttrees {

	###
	regsub "RAxML_bestTree.MSA_" $tree {} fasta
	regsub ".fas.txt" $fasta {_fam_faas.txt} fasta
	openfile $direct/$fasta
	regsub -all > $data £> data
	set glist [split $data £]
	regsub -all "{}" $glist {} glist

	set headerlist ""

	foreach gene $glist {
		set firstn [string first \n $gene]
		set header [string trim [string range $gene 0 $firstn]]
		append headerlist "$header\n"
	}
	####
	regsub "RAxML_bestTree.MSA" $tree {KEY} keyid
	regsub ".fas.txt" $keyid {_fam_faas.txt} keyid
	openfile $keyid
	set key $data

	###
	openfile $tree
	set branchids [regexp -all -inline $pata $data]
	foreach id $branchids {
		regsub {:} $id {} idtrunc
		regexp $patc$idtrunc$patd $headerlist hit
		regexp -nocase $pate $hit match
		regsub -all {\[} $match {} match
		regsub -all {\]} $match {} match
		set match [string trim $match]
		regexp $idtrunc$patb $key hit
		set hit [lindex $hit 1]
		append match " \[$hit\]"
		regsub $idtrunc $data $match data
	}
	set out [open b$tree w]
	puts $out $data
	close $out
}





###########################################################################

###RAXML OPTIONS###
"-f a": Rapid Bootstrap analysis and search for best-scoring ML tree in one program run
"-n":   Specifies the name of the output file.
"-m PROTCATmatrixName[F|X]"       : Specified AA matrix + Optimization of substitution rates + Optimization of site-specific
                                    evolutionary rates which are categorized into numberOfCategories distinct 
                                    rate categories for greater computational efficiency.   Final tree might be evaluated
"-p":	Specify a random number seed for the parsimony inferences. This allows you to reproduce your results
        and will help me debug the program.
"-x":   Specify an integer number (random seed) and turn on rapid bootstrapping
        CAUTION: unlike in version 7.0.4 RAxML will conduct rapid BS replicates under 
        the model of rate heterogeneity you specified via "-m" and not by default under CAT
"-#|-N":	Specify the number of alternative runs on distinct starting trees


exec ~/desktop/software/standard-RAxML-8.1.1/raxmlHPC-AVX -s MSA_Pot_can_murim.fas -n can565 -m PROTCATDAYHOFF -p 10001 -b 5000 -N 50

