proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

set direct ~/desktop/Fungi-Archaea-Bact

cd $direct

file mkdir Out1
file mkdir Out2
file mkdir Out3
file mkdir Out4
file mkdir Out5
file mkdir Out6
file mkdir Out7
file mkdir Out8
file mkdir Out9
file mkdir Out10
file mkdir Out11
file mkdir Out12
file mkdir Out13
file mkdir Out14
file mkdir Out15
file mkdir Out16

cd $direct/faa_all
set inlist [glob *.faa]
set splitno [expr [llength $inlist]/16]


set pata {<.*?}
set patb {.*?>}

set i 1
cd $direct/out_all
set blasthit [glob *.tsv]


foreach input $inlist {
	set blasthit [glob *.tsv]
	regsub "_protein.faa" $input {} in
	regsub -all " " $blasthit ">\n<" bb
	set bb "<$bb>"

	set xx [regexp -all -line -inline $pata$in$patb $bb]
	regsub -all < $xx {} xx
	regsub -all > $xx {} xx

	if {$i <= $splitno} {
		foreach file $xx {
			file rename -force $file $direct/Out1/$file
		}
	} elseif {$i <= [expr 2*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out2/$file
		}
	} elseif {$i <= [expr 3*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out3/$file
		}
	} elseif {$i <= [expr 4*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out4/$file
		}
	} elseif {$i <= [expr 5*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out5/$file
		}
	} elseif {$i <= [expr 6*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out6/$file
		}
	} elseif {$i <= [expr 7*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out7/$file
		}
	} elseif {$i <= [expr 8*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out8/$file
		}
	} elseif {$i <= [expr 9*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out9/$file
		}
	} elseif {$i <= [expr 10*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out10/$file
		}
	} elseif {$i <= [expr 11*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out11/$file
		}
	} elseif {$i <= [expr 12*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out12/$file
		}
	} elseif {$i <= [expr 13*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out13/$file
		}
	} elseif {$i <= [expr 14*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out14/$file
		}
	} elseif {$i <= [expr 15*$splitno]} {
		foreach file $xx {
			file rename -force $file $direct/Out15/$file
		}
	} else {
		foreach file $xx {
			file rename -force $file $direct/Out16/$file
		}
	}
	puts "$i/[llength $inlist]"
	incr i
}

