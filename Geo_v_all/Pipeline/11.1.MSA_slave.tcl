### Procedures ### 
source ~/Dropbox/Scripts/General_utils.tcl

### Variables ###
set ival 2.0
set direct /users/aesin/desktop/Geo_analysis/Geo_v_all/$ival/Missing_geobac_reconstruction
###########################################################################

cd $direct/Individual_alignment/[lindex $argv 0]

# set dirs [glob -type d *]
# set dirs [lsort -dictionary $dirs]
# reverse_dict $dirs
# set dirs $revdict

# set count_no_core 0
# set count_core 0


# foreach dir $dirs {
# 	cd $dir
	set fastas [glob *.faa]
	if {[llength $fastas] == 1} {
		set fasta [lindex $fastas 0]
		regsub {.faa} $fasta {.fas} phy_out
		regsub {CORE_} $phy_out {} phy_out
		## Output to current directory ##
		catch {exec clustalo --auto -i $fasta -o $phy_out --outfmt=phy -v --force --log=MSA_log.txt --threads=2}
		## Version with output to Group_fastas_MSA_$ival below ##
		# catch {exec clustalo --auto -i $fasta -o $direct/Group_fastas_MSA_$ival/$phy_out --outfmt=phy -v --force --log=MSA_log.txt --threads=2}
		# set keys [glob KEY*]
		# foreach key $keys {
		# 	file copy $key $direct/Group_fastas_MSA_$ival
		# }
		# incr count_no_core
	
	} elseif {[llength $fastas] == 2} {
		set core [lindex [glob CORE*] 0]
		catch {exec clustalo --full -i $core -o core_align.phylip --outfmt=phy -v --force --log=MSA_log_core.txt --threads=2}

		#### If all the sequences in the core alignment align perfectly (no gaps {-}), then clustalo will throw up an error. So, the code below adds an X onto the end of the last sequence to force a dash in the alignment. The X in the resulting alignment is then converted into a dash, resulting in a "neutral" dash column in the core alignment. ####
		
 		openfile core_align.phylip
		if {[regexp -- {-} $data] == 0} {
			openfile $core
			set data [string trim $data]
			append data {X}
			set out [open temp_$core w]
			puts $out $data
			close $out
			catch {exec clustalo --full -i temp_$core -o core_align.phylip --outfmt=phy -v --force --log=MSA_log_core.txt --threads=2}
			file delete temp_$core
			openfile core_align.phylip
			regsub -all {X} $data {-} data
			set out [open core_align.phylip w]
			puts $out $data
			close $out
		}

		#####
		set other [lindex [glob OTHER*] 0]
		regsub {.faa} $other {.fas} phy_out
		regsub {OTHER_} $phy_out {UNTRIM_} phy_out
		openfile $other
		if {[regexp -all {>} $data] == 1} {
			catch {exec clustalo --auto --p1=core_align.phylip --p2=$other -o $phy_out -v --force --log=MSA_log_full.txt --threads=2}
			# catch {exec clustalo --auto --p1=core_align.phylip --p2=$other -o $direct/Group_fastas_MSA_$ival/$phy_out -v --force --log=MSA_log_full.txt --threads=2}
		} else {
			catch {exec clustalo --auto --p1=core_align.phylip -i $other -o $phy_out -v -v --force --log=MSA_log_full.txt --threads=2}
			# catch {exec clustalo --auto --p1=core_align.phylip -i $other -o $direct/Group_fastas_MSA_$ival/$phy_out -v --force --log=MSA_log_full.txt --threads=2}
		}
		# set keys [glob KEY*]
		# foreach key $keys {
		# 	file copy $key $direct/Group_fastas_MSA_$ival
		# }
		# set ids [glob -nocomplain win_id_*]
		# foreach id $ids {
		# 	file copy -force $id $direct/Group_fastas_MSA_$ival
		# }
		# incr count_core
	}
# 	cd ..
# }

# set total_file [expr $count_no_core + $count_core]
# if {$total_file != [llength $dirs]} {
# 	puts "ERROR: Total in =/= total out"
# }
