set direct /csc/results/Warnecke/ADE/Backup
set org Geo_v_all

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################
cd $direct/$org
file mkdir Out_nw
cd $direct/$org/Out_all

set files [glob *.tsv]
set i 0
set x 1
foreach file $files {
	if {[regexp {Planococcus_halocryophilus} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Alicyclobacillus_herbarius} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Leuconostoc_pseudomesenteroides_KCTC_3652} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Thermaerobacter_subterraneus_DSM_13965} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Clostridium_tunisiense_TJ} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Oscillospiraceae_bacterium_VE202_24} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Bacillus_coahuilensis_m4_4} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Bacillus_mannanilyticus_JCM_10596} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Bacillus_vallismortis_DV1_F_3} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Lysinibacillus_boronitolerans_F1182} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Candidatus} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Clostridium_ultunense_Esp} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Clostridium_scatologenes} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Dehalobacter_sp._FTH1} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Dehalobacter_sp._E1} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Lactobacillus_apodemi_DSM_16634_JCM_16172} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Lactobacillus_fuchuensis_DSM_14340_JCM_11249} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Leuconostoc_fallax_KCTC_3537} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	} elseif {[regexp {Leuconostoc_inhae_KCTC_3774} $file] == 1} {
		file rename $file $direct/$org/Out_nw/
		incr i
	}
	puts $x
	incr x
}
	
puts "\n$i files moved out to Out_nw"
