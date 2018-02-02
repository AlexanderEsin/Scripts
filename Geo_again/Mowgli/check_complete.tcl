set penalty		4
set direct		/scratch2/ade110/Geo_again/Mowgli

set in_direct	$direct/GeneTree_input
set out_direct	$direct/Mowgli_output/Output_$penalty

set input_dirs	[glob -type d $in_direct/*]
set recon_dirs	[glob -type d $out_direct/*]


set incomp_recon	0
set complete_recon	0

foreach recon_dir $recon_dirs {
	if {[file exists $recon_dir/Fullreconciliation.mpr] == 0} {
		error("Cannot find reconciliation file in $recon_dir")
	}

	set output_size	[file size $recon_dir/Fullreconciliation.mpr]
	if {$output_size == 0} {
		incr incomp_recon
	} else {
		incr complete_recon
	}
}

# 

puts stdout "Total number of input gene trees: [llength $input_dirs]"
puts stdout "Total number of reconciliations attempted: [llength $recon_dirs]"
puts stdout "Total number of reconciliations succeeded: $complete_recon"
puts stdout "Total number of reconciliations failed: $incomp_recon"








## /// At 50G RAM /// ##

## T3 
# Total number of input gene trees: 6712
# Total number of reconciliations attempted: 6712
# Total number of reconciliations succeeded: 6357
# Total number of reconciliations failed: 355

## T4
# Total number of input gene trees: 6712
# Total number of reconciliations attempted: 6712
# Total number of reconciliations succeeded: 6170
# Total number of reconciliations failed: 542

## T5
# Total number of input gene trees: 6712
# Total number of reconciliations attempted: 6712
# Total number of reconciliations succeeded: 6218
# Total number of reconciliations failed: 494


## /// At 100G RAM /// ##
## T3 
# Total number of input gene trees: 6712
# Total number of reconciliations attempted: 6712
# Total number of reconciliations succeeded: 6560
# Total number of reconciliations failed: 152







































