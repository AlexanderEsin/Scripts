proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}
array set translate {
	     atg M
	     ttt F
	     ttc F
	     tta L
	     ttg L
	     tct S
	     tcc S
	     tca S
	     tcg S
	     tat Y
	     tac Y
	     taa X
	     tag X
	     tgt C
	     tgc C
	     tga X
	     tgg W
	     ctt L
	     ctc L
	     cta L
	     ctg L
	     cct P
	     ccc P
	     cca P
	     ccg P
	     cat H
	     cac H
	     caa Q
	     cag Q
	     cgt R
	     cgc R
	     cga R
	     cgg R
	     att I
	     atc I
	     ata I
	     act T
	     acc T
	     aca T
	     acg T
	     aat N
	     aac N
	     aaa K
	     aag K
	     agt S
	     agc S
	     aga R
	     agg R
	     gtt V
	     gtc V
	     gta V
	     gtg V
	     gct A
	     gcc A
	     gca A
	     gcg A
	     gat D
	     gac D
	     gaa E
	     gag E
	     ggt G
	     ggc G
	     gga G
	     ggg G}

openfile Drosophila_notch.fa

		regsub -all {>} $data £> data
		set glist [split $data £]
		set glist [lrange $glist 1 end]

		foreach gene $glist {
			set ori [string first \n $gene]
			set cds [string range $gene $ori end]
			set header [string trim [string range $gene 0 $ori]]
			regsub -all {[^A-Z]} $cds {} cds
			set prot ""
			while {[string length $cds]>=3} {
				set codon [string tolower [string range $cds 0 2]]
				if {[regexp {[^actg]} $codon]} {
					append prot X
				} else {
				append prot $translate($codon)
				}
				set cds [string range $cds 3 end]
			}
			regsub {X} $prot {} prot
			set out [open Drosophila_notch_aa.fa a]
			puts $out "$header\n$prot"
			close $out

						}

