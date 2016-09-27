#!/usr/bin/perl
##########################################
#                                        #
#      TOPD-fMtS version 3.3             #
#      http://genomes.urv.es/topd        #
#      (August 2007)                     #
#                                        #
#      Author: Pere Puigbo               #
#      E-mail: ppuigbo@urv.net           #
#                                        #
##########################################

&credits();
&gethelp($ARGV[0]);
$menu = $ARGV[0];

if ($menu eq "-species") {
	my ($SP, $file_trees) = @ARGV;
	&SP($file_trees);
	exit;
}

if ($menu eq "-prune") {
	my ($prunne, $file_trees, $file_species) = @ARGV;
	&prunne_trees($file_trees, $file_species);
	exit;
}

if ($menu eq "-fmts") {

	my ($prunne, $file_trees, $type) = @ARGV;
	&fmts_trees($file_trees, $type);
	exit;
}

if ($menu eq "-quartets") {
	my ($prunne, $file_trees) = @ARGV;
	&getQuartets($file_trees);
	exit;
}

if ($menu eq "-triplets") {
	my ($prunne, $file_trees) = @ARGV;
	&getTriplets($file_trees);
	exit;
}

if ($menu eq "-qq") {
	my ($prunne, $file_trees, $nquartets) = @ARGV;
	&quickquartets($file_trees, $nquartets);
	exit;
}

if ($menu eq "-dd") {
	my ($prunne, $file_trees) = @ARGV;
	&duplidist($file_trees);
	exit;
}

if ($menu eq "") {
	&gethelp();
	exit;
}else {
	@arguments = @ARGV;
}

($LocTrees, $randomTrees, $nombreRandomTrees, $nombreSGTgenerats, $SingleMultiple, $print_prev, $fileoutput) = &getParameters(@arguments);

sleep(2);

if ($SingleMultiple eq "single") {
	&menu($LocTrees, $randomTrees, $nombreRandomTrees, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $fileoutput,$units);
}elsif ($SingleMultiple eq "multiple") {
	&menu2($LocTrees, $randomTrees, $nombreRandomTrees, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $fileoutput,$units);
}elsif ($SingleMultiple eq "reference") {
	&menu3($LocTrees, $randomTrees, $nombreRandomTrees, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $fileoutput,$units);
}

sub gethelp() {
	my $help = shift;
	
	if ( ($help eq "-help") || ($help eq "-h") || ($help eq "-HELP") || ($help eq "-Help") || ($help eq "") ){
		print "\n\n####################################################################################################################\n\n";
		print "Parameters to run TOPD/fMtS:\n\n";
		print "\t\$\/topd_v* -f [file] -m [nodal/split/i/triplets/all] -r [yes/no] -n [1-1000] -s [10-10000/all] -c [single/multiple/reference] -p [yes/no] -help\n\n";
		print "\tINPUT FILE NAME			-f	<file_name>\n";
		print "\tOUTPUT FILE NAME			-out	<file_name>\n";
		print "\tMETHOD				-m	<nodal/split/quartets/triplets/disagree/all>	default: split\n";
		print "\tNUMBER OF TRIPLETS&QUARTETS	-u	<all/random/relative>	default: all\n";
		print "\tMAXIMUM LEVEL IN THE DISAGREE METHOD -l <1/2/3/4>	default: 1\n";
		print "\tRANDOM ANALYSIS			-r	<random/guided/no>			default: random\n";
		print "\tNUMBER OF RANDOM SEQUENCES	-n	<1-1000>				default: 100\n";
		print "\tNUMBER OF fMtS SEQUENCES	-s	<10-10000/all/relative>			default: relative\n";
		print "\tCOMPARISONS			-c	<single/multiple/reference>		default: single\n";
		print "\tPRINT PREVIOUS RESULTS		-p	<yes/no>				default: yes\n";
		print "\t-------------------------------------------------------------------------------------------------------------------\n";
		print "\tGET SPECIES FROM TREES		-species	<file_name>\n";
		print "\tPRUNE TREES			-prune		<file_trees>	<file_species_to_prunne>\n";
		print "\tFMTS TREES			-fmts		<file_trees>	<random/all/relative>	default: relative\n";
		print "\tGET QUARTETS			-quartets	<file_trees>\n";
		print "\tGET TRIPLETS			-triplets	<file_trees>\n";
		print "\tDUPLIDIST			-dd		<file_trees>\n";
		print "\tQUICK QUARTETS			-qq		<file_trees>	<number_quartets>	default: 100\n";
		print "\n####################################################################################################################\n\n";
		exit;
	}
}



sub getParameters() {
	%parameters = @_;
	print "\nReading parameters ...\n\n";
	
	if ($parameters{'-f'} eq "") {
		$LocTrees = $parameters{'-f'};
		print "-f $LocTrees .......................................... error";
		die "#ERROR# Please insert the name of the input file containing the trees\n";
	}else {
		$LocTrees = $parameters{'-f'};
		print "-f $LocTrees .......................................... ok\n"
	}

	
	if ($parameters{'-out'} eq "") {
		$fileoutput = "";
	}else {
		$fileoutput = $parameters{'-out'};
		print "-out $fileoutput ...................................... ok\n";
	}

	my $nombreSGTgenerats = "";
	my $nombreRandomTrees = "";


	if ($parameters{'-m'} eq "") {
		$method = "split";
		print "-m $method .......................................... default\n"
	}else {
		$method = $parameters{'-m'};
		if ($method eq 'nodal') {
			$method = 'nodal';
			print "-m $method .......................................... ok\n"
		}
		elsif ($method eq 'split') {
			$method = 'split';
			print "-m $method .......................................... ok\n"
		}
		elsif ($method eq 'quartets') {
			$method = 'quartets';
			print "-m $method .......................................... ok\n"
		}
		elsif ($method eq 'triplets') {
			$method = 'triplets';
			print "-m $method .......................................... ok\n"
		}
		elsif ($method eq 'disagree') {
			$method = 'disagree';
			print "-m $method .......................................... ok\n"
		}
		else {
			$method = "all";
			print "-m $method .......................................... default\n"
		}
	}

	
	if ($parameters{'-l'} eq "") {
		$level = "1";
		print "-l $level .......................................... default\n"
	}else {
		$level = $parameters{'-l'};
		if ( $level < 1 ) {
			die "Fail parameter '-l'. Please, insert a number between 1 to 4\n";
		}
		elsif ( $nombreRandomTrees > 1000  ) {
			die "Fail parameter '-l'. Please, insert a number between 1 to 4\n";
		}else {
			print "-l $level .......................................... ok\n"
		}
	}

	
	if ($parameters{'-u'} eq "") {
		$units = "all";
		print "-u $units .......................................... default\n"
	}else {
		$units = $parameters{'-u'};
		if ( ($units ne "random") && ($units ne "relative") && ($units ne "all") ) {
			$units = "all";
		}
	}


	
	if ($parameters{'-r'} eq "") {
		$randomTrees = "random";
		print "-r $randomTrees .......................................... default\n"
	}else {
		$randomTrees = $parameters{'-r'};
		if ($randomTrees eq 'random') {
			$randomTrees = 'random';
			print "-r $randomTrees .......................................... ok\n"
		}
		elsif ($randomTrees eq 'guided') {
			$randomTrees = 'guided';
			print "-r $randomTrees .......................................... ok\n"
		}
		elsif ( ($randomTrees eq 'n') || ($randomTrees eq 'N') || ($randomTrees eq 'NO') || ($randomTrees eq 'no') || ($randomTrees eq 'No')) {
			$randomTrees = 'no';
			print "-r $randomTrees .......................................... ok\n"
		}
		else {
			$randomTrees = "random";
			print "-r $randomTrees .......................................... default\n"
		}
	}

	
	if ($parameters{'-n'} eq "") {
		$nombreRandomTrees = "100";
		print "-n $nombreRandomTrees .......................................... default\n"
	}else {
		$nombreRandomTrees = $parameters{'-n'};
		if ( $nombreRandomTrees < 1 ) {
			die "Fail parameter '-n'. Please, insert a number between 10 to 1000\n";
		}
		elsif ( $nombreRandomTrees > 1000  ) {
			die "Fail parameter '-n'. Please, insert a number between 10 to 1000\n";
		}else {
			print "-n $nombreRandomTrees .......................................... ok\n"
		}
	}

	
	if ($parameters{'-c'} eq "") {
		$SingleMultiple = "single";
		print "-c $SingleMultiple .......................................... default\n";
	}else {
		$SingleMultiple = $parameters{'-c'};
		if ( ($SingleMultiple eq 'single') || ($SingleMultiple eq 'Single') || ($SingleMultiple eq 's') || ($SingleMultiple eq 'S') || ($SingleMultiple eq 'SINGLE') ) {
			$SingleMultiple = 'single';
			print "-c $SingleMultiple .......................................... ok\n";
		}
		elsif ( ($SingleMultiple eq 'multiple') || ($SingleMultiple eq 'Multiple') || ($SingleMultiple eq 'm') || ($SingleMultiple eq 'M') || ($SingleMultiple eq 'MULTIPLE') ) {
			$SingleMultiple = 'multiple';	
			print "-c $SingleMultiple .......................................... ok\n";
		}
		elsif ( ($SingleMultiple eq 'reference') || ($SingleMultiple eq 'Reference') || ($SingleMultiple eq 'r') || ($SingleMultiple eq 'R') || ($SingleMultiple eq 'REFERENCE') ) {
			$SingleMultiple = 'reference';	
			print "-c $SingleMultiple .......................................... ok\n";
		}else {
			$SingleMultiple = "single";
			print "-c $SingleMultiple .......................................... default\n";
		}
	}


	
	if ($parameters{'-s'} eq "") {
		$nombreSGTgenerats = "relative";
		print "-s $nombreSGTgenerats .......................................... default\n"
	} elsif ( ($parameters{'-s'} eq "A") || ($parameters{'-s'} eq "a") ||($parameters{'-s'} eq "all") ||  ($parameters{'-s'} eq "All") ||  ($parameters{'-s'} eq "ALL") ) {
		$nombreSGTgenerats = "all";
		print "-s $nombreSGTgenerats .......................................... ok\n"
	}
	elsif ( ($parameters{'-s'} eq "R") || ($parameters{'-s'} eq "r") ||($parameters{'-s'} eq "relative") ||  ($parameters{'-s'} eq "Relative") ||  ($parameters{'-s'} eq "RELATIVE") ) {
		$nombreSGTgenerats = "relative";
		print "-s $nombreSGTgenerats .......................................... ok\n"
	}
	else {
		$nombreSGTgenerats = $parameters{'-s'};
		if ( $nombreSGTgenerats < 1 ) {
			die "Fail parameter '-s'. Please, insert a number between 10 to 10000\n";
		}
		elsif ( $nombreSGTgenerats > 100000  ) {
			die "Fail parameter '-s'. Please, insert a number between 10 to 10000\n";
		}else {
			print "-s $nombreSGTgenerats .......................................... ok\n"
		}
	}
	
	
	if ($parameters{'-p'} eq "") {
		$print_prev = "y";
		print "-p $print_prev .......................................... default\n"	
	}
	elsif ( ($parameters{'-p'} eq "y") || ($parameters{'-p'} eq "Y") || ($parameters{'-p'} eq "YES") || ($parameters{'-p'} eq "yes") || ($parameters{'-p'} eq "Yes")  ){
		$print_prev = "y";
		print "-p $print_prev .......................................... ok\n"
	}
	elsif ( ($parameters{'-p'} eq "n") || ($parameters{'-p'} eq "N") || ($parameters{'-p'} eq "NO") || ($parameters{'-p'} eq "no") || ($parameters{'-p'} eq "No")  ){
		$print_prev = "n";
		print "-p $print_prev .......................................... ok\n"
	}
	else {
		$print_prev = "s";
		print "-p $print_prev .......................................... default\n"	
	}
	return ($LocTrees, $randomTrees, $nombreRandomTrees, $nombreSGTgenerats, $SingleMultiple, $print_prev, $fileoutput,$units);
}


sub menu() { 
	($LocTrees, $randomTrees, $nombreRandomTrees, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $fileoutput,$units) = @_;
	
	my ($name1_IN, $tree1_IN, $name2_IN, $tree2_IN) = &getTrees($LocTrees);

	
	my ($Nsp1, $family1, $Nsp2, $family2) = &MoS($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);

	if ( ($family1 eq "S") && ($family2 eq "S") ) { 	
		
		($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $tree1_IN, $name2_IN, $tree2_IN, $print_prev, $method, $level, $units);
	}
	elsif ( ($family1 eq "S") && ($family2 eq "M") ) {
		
		
		($resp, $ptree1_in, $ptree2_in) = &checkPruneMGF($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);

		
		if ($resp eq "S") {
			if ($print_prev eq "n") {
			}else {
				print "\nPrunned trees are both single gene trees.\n";
			}
			($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $ptree1_in, $name2_IN, $ptree2_in, $print_prev, $method, $level, $units);
			($com, $pcCom2) = &pcCom($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
		}else {
			($topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori, $com, $pcCom2, $SplitDistSM, $SplitDistSMsd, $qSM, $qSMsd, $mSM, $mSMsd, $dif_elem, $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd, $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd, $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd) = &topdSM($name1_IN, $tree1_IN, $name2_IN, $tree2_IN, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);

		}
	}elsif ( ($family1 eq "M") && ($family2 eq "S") ) {
		
		
		($resp, $ptree2_in, $ptree1_in) = &checkPruneMGF($name2_IN, $tree2_IN, $name1_IN, $tree1_IN);
		
		if ($resp eq "S") { 
			if ($print_prev eq "n") {
			}else {
				print "\nPrunned trees are both single gene trees.\n";
			} 
			($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $ptree1_in, $name2_IN, $ptree2_in, $print_prev, $method, $level, $units);
			($com, $pcCom2) = &pcCom($name2_IN, $tree2_IN, $name1_IN, $tree1_IN);
		}else { 
			($topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori, $com, $pcCom2, $SplitDistSM, $SplitDistSMsd, $qSM, $qSMsd, $mSM, $mSMsd, $dif_elem, $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd, $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd, $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd) = &topdSM($name2_IN, $tree2_IN, $name1_IN, $tree1_IN, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
		}
	}elsif ( ($family1 eq "M") && ($family2 eq "M") ) {
		
		
		($resp, $ptree1_in, $ptree2_in) = &checkPruneMGF($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);

		
		if ($resp eq "S") {
			if ($print_prev eq "n") {
			}else {
				print "\nPrunned trees are both single gene trees.\n";
			}
			($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $ptree1_in, $name2_IN, $ptree2_in, $print_prev, $method, $level, $units);
			($com, $pcCom2) = &pcCom($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
		}else {
			($topdMM_ori, $SD_ori, $topdMMunpruned_ori, $SDunpruned_ori, $com, $pcCom2, $SplitDistMM, $SplitDistMMsd, $qMM, $qMMsd, $mMM, $mMMsd, $dif_elem, $num_sp_elimMM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1MMsd, $minMM, $minMMsd, $nquartMM, $nquartMMsd, $qsolvedMM, $qsolvedMMsd, $qdifferentMM, $qdifferentMMsd, $qr2MM, $qr2MMsd, $qr1MM, $qr1MMsd, $qUnresolvedMM, $qUnresolvedMMsd, $qDCMM, $qDCMMsd, $qEAMM, $qEAMMsd, $qSJAMM, $qSJAMMsd, $qSSJAMM, $qSSJAMMsd, $qSDMM, $qSDMMsd, $ntriplMM, $ntriplMMsd, $tsolvedMM, $tsolvedMMsd, $tdifferentMM, $tdifferentMMsd,$tr2MM, $tr2MMsd, $tr1MM, $tr1MMsd, $tUnresolvedMM, $tUnresolvedMMsd, $tDCMM, $tDCMMsd, $tEAMM, $tEAMMsd, $tSJAMM, $tSJAMMsd, $tSSJAMM, $tSSJAMMsd, $tSDMM, $tSDMMsd) = &topdMM($name1_IN, $tree1_IN, $name2_IN, $tree2_IN, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
		}
	}

	if ($fileoutput eq "") {
		open (OUTa, ">topd_result.txt");
	}else {
		open (OUTa, ">$fileoutput");
	}
	if ($com > 3) {
		if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
			
			
			my @topd_rand =();
			my @topd_unpruned_rand = ();
			my @split_dist_rand = ();
			my @Q_rand = ();
			my @M_rand = ();

			my @A_num_sp_elim_random = ();
			my @A_e_n_taxa_tree1_random = (); 
			my @A_min_random = (); 
			my %A_DisTax_random = ();
			my @A_nquart_random = ();
			my @A_qsolved_random = (); 
			my @A_qdifferent_random = (); 
			my @A_qr2_random = (); 
			my @A_qr1_random = (); 
			my @A_qUnresolved_random = (); 
			my @A_qDC_random = (); 
			my @A_qEA_random = (); 
			my @A_qSJA_random = (); 
			my @A_qSSJA_random = ();
			my @A_qSD_random = (); 
			my @A_ntripl_random = (); 
			my @A_tsolved_random = ();
			my @A_tdifferent_random = (); 
			my @A_tr1_random = (); 
			my @A_tr2_random = (); 
			my @A_tUnresolved_random = (); 
			my @A_tDC_random = (); 
			my @A_tEA_random = (); 
			my @A_tSJA_random = (); 
			my @A_tSSJA_random = (); 
			my @A_tSD_random = ();

			for ($ir=0;$ir<$nombreRandomTrees;$ir++) {
				if ($print_prev eq "n") {
				}else {
					print "\n### RANDOM $ir ###";
				}
				my ($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom) = ();
				if ($randomTrees eq 'guided') {
					($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom) = &random($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
				} elsif ($randomTrees eq 'random') {
					($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom) = &random2($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
				}
				
				my ($Nsp1random, $family1random, $Nsp2random, $family2random) = &MoS($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom);
	
				if ( ($family1random eq "S") && ($family2random eq "S") ) {
					
					my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom, $print_prev, $method, $level, $units);
					push @topd_rand, $topdrandom;
					push @topd_unpruned_rand, $topdUnprunedrandom;
					push @split_dist_rand, $SplitDistrandom;
					push @Q_rand, $qrandom;
					push @M_rand, $mrandom;

					push @A_num_sp_elim_random, $num_sp_elimrandom;
					push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
					push @A_min_random, $minrandom;
					$A_DisTax_random{$DisTaxrandom}++;
					push @A_nquart_random, $nquartrandom;
					push @A_qsolved_random, $qsolvedrandom; 
					push @A_qdifferent_random, $qdifferentrandom; 
					push @A_qr2_random, $qr2random; 
					push @A_qr1_random, $qr1random; 
					push @A_qUnresolved_random, $qUnresolvedrandom; 
					push @A_qDC_random, $qDCrandom; 
					push @A_qEA_random, $qEArandom; 
					push @A_qSJA_random, $qSJArandom; 
					push @A_qSSJA_random, $qSSJArandom;
					push @A_qSD_random, $qSDrandom; 
					push @A_ntripl_random, $ntriplrandom; 
					push @A_tsolved_random, $tsolvedrandom;
					push @A_tdifferent_random, $tdifferentrandom; 
					push @A_tr1_random, $tr1random; 
					push @A_tr2_random, $tr2random; 
					push @A_tUnresolved_random, $tUnresolvedrandom; 
					push @A_tDC_random, $tDCrandom; 
					push @A_tEA_random, $tEArandom; 
					push @A_tSJA_random, $tSJArandom; 
					push @A_tSSJA_random, $tSSJArandom; 
					push @A_tSD_random, $tSDrandom;
				}
				elsif ( ($family1random eq "S") && ($family2random eq "M") ) {
					
					
					($resp, $ptree1_in_r, $ptree2_in_r) = &checkPruneMGF($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom);
	
					
					if ($resp eq "S") {
						if ($print_prev eq "n") {
						}else {
							print "\nPrunned trees are both single gene trees.\n";
						}
						my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $ptree1_in_r, $name2_INrandom, $ptree2_in_r, $print_prev, $method, $level, $units);
						push @topd_rand, $topdrandom;
						push @topd_unpruned_rand, $topdUnprunedrandom;
						push @split_dist_rand, $SplitDistrandom;
						push @Q_rand, $qrandom;
						push @M_rand, $mrandom;

						push @A_num_sp_elim_random, $num_sp_elimrandom;
						push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
						push @A_min_random, $minrandom;
						$A_DisTax_random{$DisTaxrandom}++;
						push @A_nquart_random, $nquartrandom;
						push @A_qsolved_random, $qsolvedrandom; 
						push @A_qdifferent_random, $qdifferentrandom; 
						push @A_qr2_random, $qr2random; 
						push @A_qr1_random, $qr1random; 
						push @A_qUnresolved_random, $qUnresolvedrandom; 
						push @A_qDC_random, $qDCrandom; 
						push @A_qEA_random, $qEArandom; 
						push @A_qSJA_random, $qSJArandom; 
						push @A_qSSJA_random, $qSSJArandom;
						push @A_qSD_random, $qSDrandom; 
						push @A_ntripl_random, $ntriplrandom; 
						push @A_tsolved_random, $tsolvedrandom;
						push @A_tdifferent_random, $tdifferentrandom; 
						push @A_tr1_random, $tr1random; 
						push @A_tr2_random, $tr2random; 
						push @A_tUnresolved_random, $tUnresolvedrandom; 
						push @A_tDC_random, $tDCrandom; 
						push @A_tEA_random, $tEArandom; 
						push @A_tSJA_random, $tSJArandom; 
						push @A_tSSJA_random, $tSSJArandom; 
						push @A_tSD_random, $tSDrandom;
					}else {
						($topdSMrandom, $SDrandom, $topdSMunprunedrandom, $SDunprunedrandom, $com, $pcCom2, $SplitDistrandom, $SplitDistsdrandom, $qrandom, $qsdrandom, $mrandom, $msdrandom, $DisTaxrandom, $num_sp_elimSMrandom, $num_sp_elimSMsdrandom, $e_n_taxa_tree1SMrandom, $e_n_taxa_tree1SMsdrandom, $minSMrandom, $minSMsdrandom, $nquartSMrandom, $nquartSMsdrandom, $qsolvedSMrandom, $qsolvedSMsdrandom, $qdifferentSMrandom, $qdifferentSMsdrandom, $qr2SMrandom, $qr2SMsdrandom, $qr1SMrandom, $qr1SMsdrandom, $qUnresolvedSMrandom, $qUnresolvedSMsdrandom, $qDCSMrandom, $qDCSMsdrandom, $qEASMrandom, $qEASMsdrandom, $qSJASMrandom, $qSJASMsdrandom, $qSSJASMrandom, $qSSJASMsdrandom, $qSDSMrandom, $qSDSMsdrandom, $ntriplSMrandom, $ntriplSMsdrandom, $tsolvedSMrandom, $tsolvedSMsdrandom, $tdifferentSMrandom, $tdifferentSMsdrandom,$tr2SMrandom, $tr2SMsdrandom, $tr1SMrandom, $tr1SMsdrandom, $tUnresolvedSMrandom, $tUnresolvedSMsdrandom, $tDCSMrandom, $tDCSMsdrandom, $tEASMrandom, $tEASMsdrandom, $tSJASMrandom, $tSJASMsdrandom, $tSSJASMrandom, $tSSJASMsdrandom, $tSDSMrandom, $tSDSMsdrandom) = &topdSM($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
						push @topd_rand, $topdSMrandom;
						push @topd_unpruned_rand, $topdSMunprunedrandom;
						push @split_dist_rand, $SplitDistrandom;
						push @Q_rand, $qrandom; 
						push @M_rand, $mrandom;

						push @A_num_sp_elim_random, $num_sp_elimSMrandom;
						push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1SMrandom;
						push @A_min_random, $minSMrandom;
						$A_DisTax_random{$DisTaxrandom}++;
						push @A_nquart_random, $nquartSMrandom;
						push @A_qsolved_random, $qsolvedSMrandom; 
						push @A_qdifferent_random, $qdifferentSMrandom; 
						push @A_qr2_random, $qr2SMrandom; 
						push @A_qr1_random, $qr1SMrandom; 
						push @A_qUnresolved_random, $qUnresolvedSMrandom; 
						push @A_qDC_random, $qDCSMrandom; 
						push @A_qEA_random, $qEASMrandom; 
						push @A_qSJA_random, $qSJASMrandom; 
						push @A_qSSJA_random, $qSSJASMrandom;
						push @A_qSD_random, $qSDSMrandom; 
						push @A_ntripl_random, $ntriplSMrandom; 
						push @A_tsolved_random, $tsolvedSMrandom;
						push @A_tdifferent_random, $tdifferentSMrandom; 
						push @A_tr1_random, $tr1SMrandom; 
						push @A_tr2_random, $tr2SMrandom; 
						push @A_tUnresolved_random, $tUnresolvedSMrandom; 
						push @A_tDC_random, $tDCSMrandom; 
						push @A_tEA_random, $tEASMrandom; 
						push @A_tSJA_random, $tSJASMrandom; 
						push @A_tSSJA_random, $tSSJASMrandom; 
						push @A_tSD_random, $tSDSMrandom;
					}	
				}elsif ( ($family1random eq "M") && ($family2random eq "S") ) {
					
	
					
					($resp, $ptree1_in_r, $ptree2_in_r) = &checkPruneMGF($name2_INrandom, $tree2_INrandom, $name1_INrandom, $tree1_INrandom);
					
					
					if ($resp eq "S") {
						if ($print_prev eq "n") {
						}else {
							print "\nPrunned trees are both single gene trees.\n";
						}
						my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $ptree1_in_r, $name2_INrandom, $ptree2_in_r, $print_prev, $method, $level, $units);
						push @topd_rand, $topdrandom;
						push @topd_unpruned_rand, $topdUnprunedrandom;
						push @split_dist_rand, $SplitDistrandom;
						push @Q_rand, $qrandom; 
						push @M_rand, $mrandom;

						push @A_num_sp_elim_random, $num_sp_elimrandom;
						push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
						push @A_min_random, $minrandom;
						$A_DisTax_random{$DisTaxrandom}++;
						push @A_nquart_random, $nquartrandom;
						push @A_qsolved_random, $qsolvedrandom; 
						push @A_qdifferent_random, $qdifferentrandom; 
						push @A_qr2_random, $qr2random; 
						push @A_qr1_random, $qr1random; 
						push @A_qUnresolved_random, $qUnresolvedrandom; 
						push @A_qDC_random, $qDCrandom; 
						push @A_qEA_random, $qEArandom; 
						push @A_qSJA_random, $qSJArandom; 
						push @A_qSSJA_random, $qSSJArandom;
						push @A_qSD_random, $qSDrandom; 
						push @A_ntripl_random, $ntriplrandom; 
						push @A_tsolved_random, $tsolvedrandom;
						push @A_tdifferent_random, $tdifferentrandom; 
						push @A_tr1_random, $tr1random; 
						push @A_tr2_random, $tr2random; 
						push @A_tUnresolved_random, $tUnresolvedrandom; 
						push @A_tDC_random, $tDCrandom; 
						push @A_tEA_random, $tEArandom; 
						push @A_tSJA_random, $tSJArandom; 
						push @A_tSSJA_random, $tSSJArandom; 
						push @A_tSD_random, $tSDrandom;
					}else {
						($topdSMrandom, $SDrandom, $topdSMunprunedrandom, $SDunprunedrandom, $com, $pcCom2, $SplitDistrandom, $SplitDistsdrandom, $qrandom, $qsdrandom, $mrandom, $msdrandom, $DisTaxrandom, $num_sp_elimSMrandom, $num_sp_elimSMsdrandom, $e_n_taxa_tree1SMrandom, $e_n_taxa_tree1SMsdrandom, $minSMrandom, $minSMsdrandom, $nquartSMrandom, $nquartSMsdrandom, $qsolvedSMrandom, $qsolvedSMsdrandom, $qdifferentSMrandom, $qdifferentSMsdrandom, $qr2SMrandom, $qr2SMsdrandom, $qr1SMrandom, $qr1SMsdrandom, $qUnresolvedSMrandom, $qUnresolvedSMsdrandom, $qDCSMrandom, $qDCSMsdrandom, $qEASMrandom, $qEASMsdrandom, $qSJASMrandom, $qSJASMsdrandom, $qSSJASMrandom, $qSSJASMsdrandom, $qSDSMrandom, $qSDSMsdrandom, $ntriplSMrandom, $ntriplSMsdrandom, $tsolvedSMrandom, $tsolvedSMsdrandom, $tdifferentSMrandom, $tdifferentSMsdrandom,$tr2SMrandom, $tr2SMsdrandom, $tr1SMrandom, $tr1SMsdrandom, $tUnresolvedSMrandom, $tUnresolvedSMsdrandom, $tDCSMrandom, $tDCSMsdrandom, $tEASMrandom, $tEASMsdrandom, $tSJASMrandom, $tSJASMsdrandom, $tSSJASMrandom, $tSSJASMsdrandom, $tSDSMrandom, $tSDSMsdrandom) = &topdSM($name2_INrandom, $tree2_INrandom, $name1_INrandom, $tree1_INrandom, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);

						push @topd_rand, $topdSMrandom;
						push @topd_unpruned_rand, $topdSMunprunedrandom;
						push @split_dist_rand, $SplitDistrandom;
						push @Q_rand, $qrandom;
						push @M_rand, $mrandom;

						push @A_num_sp_elim_random, $num_sp_elimSMrandom;
						push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1SMrandom;
						push @A_min_random, $minSMrandom;
						$A_DisTax_random{$DisTaxrandom}++;
						push @A_nquart_random, $nquartSMrandom;
						push @A_qsolved_random, $qsolvedSMrandom; 
						push @A_qdifferent_random, $qdifferentSMrandom; 
						push @A_qr2_random, $qr2SMrandom; 
						push @A_qr1_random, $qr1SMrandom; 
						push @A_qUnresolved_random, $qUnresolvedSMrandom; 
						push @A_qDC_random, $qDCSMrandom; 
						push @A_qEA_random, $qEASMrandom; 
						push @A_qSJA_random, $qSJASMrandom; 
						push @A_qSSJA_random, $qSSJASMrandom;
						push @A_qSD_random, $qSDSMrandom; 
						push @A_ntripl_random, $ntriplSMrandom; 
						push @A_tsolved_random, $tsolvedSMrandom;
						push @A_tdifferent_random, $tdifferentSMrandom; 
						push @A_tr1_random, $tr1SMrandom; 
						push @A_tr2_random, $tr2SMrandom; 
						push @A_tUnresolved_random, $tUnresolvedSMrandom; 
						push @A_tDC_random, $tDCSMrandom; 
						push @A_tEA_random, $tEASMrandom; 
						push @A_tSJA_random, $tSJASMrandom; 
						push @A_tSSJA_random, $tSSJASMrandom; 
						push @A_tSD_random, $tSDSMrandom;
					}	
				}elsif ( ($family1random eq "M") && ($family2random eq "M") ) {
					
					
					($resp, $ptree1_in_r, $ptree2_in_r) = &checkPruneMGF($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom);
					
					if ($resp eq "S") {
						if ($print_prev eq "n") {
						}else {
							print "\nPrunned trees are both single gene trees.\n";
						}
						my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $ptree1_in_r, $name2_INrandom, $ptree2_in_r, $print_prev, $method, $level, $units);
						push @topd_rand, $topdrandom;
						push @topd_unpruned_rand, $topdUnprunedrandom;
						push @split_dist_rand, $SplitDistrandom;
						push @Q_rand, $qrandom;
						push @M_rand, $mrandom;

						push @A_num_sp_elim_random, $num_sp_elimrandom;
						push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
						push @A_min_random, $minrandom;
						$A_DisTax_random{$DisTaxrandom}++;
						push @A_nquart_random, $nquartrandom;
						push @A_qsolved_random, $qsolvedrandom; 
						push @A_qdifferent_random, $qdifferentrandom; 
						push @A_qr2_random, $qr2random; 
						push @A_qr1_random, $qr1random; 
						push @A_qUnresolved_random, $qUnresolvedrandom; 
						push @A_qDC_random, $qDCrandom; 
						push @A_qEA_random, $qEArandom; 
						push @A_qSJA_random, $qSJArandom; 
						push @A_qSSJA_random, $qSSJArandom;
						push @A_qSD_random, $qSDrandom; 
						push @A_ntripl_random, $ntriplrandom; 
						push @A_tsolved_random, $tsolvedrandom;
						push @A_tdifferent_random, $tdifferentrandom; 
						push @A_tr1_random, $tr1random; 
						push @A_tr2_random, $tr2random; 
						push @A_tUnresolved_random, $tUnresolvedrandom; 
						push @A_tDC_random, $tDCrandom; 
						push @A_tEA_random, $tEArandom; 
						push @A_tSJA_random, $tSJArandom; 
						push @A_tSSJA_random, $tSSJArandom; 
						push @A_tSD_random, $tSDrandom;
					}else {
						($topdMMrandom, $SDrandom, $topdMMunprunedrandom, $SDunprunedrandom, $com, $pcCom2, $SplitDistrandom, $SplitDistsdrandom, $qrandom, $qsdrandom, $mrandom, $msdrandom, $DisTaxrandom, $num_sp_elimMMrandom, $num_sp_elimMMsdrandom, $e_n_taxa_tree1MMrandom, $e_n_taxa_tree1MMsdrandom, $minMMrandom, $minMMsdrandom, $nquartMMrandom, $nquartMMsdrandom, $qsolvedMMrandom, $qsolvedMMsdrandom, $qdifferentMMrandom, $qdifferentMMsdrandom, $qr2MMrandom, $qr2MMsdrandom, $qr1MMrandom, $qr1MMsdrandom, $qUnresolvedMMrandom, $qUnresolvedMMsdrandom, $qDCMMrandom, $qDCMMsdrandom, $qEAMMrandom, $qEAMMsdrandom, $qSJAMMrandom, $qSJAMMsdrandom, $qSSJAMMrandom, $qSSJAMMsdrandom, $qSDMMrandom, $qSDMMsdrandom, $ntriplMMrandom, $ntriplMMsdrandom, $tsolvedMMrandom, $tsolvedMMsdrandom, $tdifferentMMrandom, $tdifferentMMsdrandom, $tr2MMrandom, $tr2MMsdrandom, $tr1MMrandom, $tr1MMsdrandom, $tUnresolvedMMrandom, $tUnresolvedMMsdrandom, $tDCMMrandom, $tDCMMsdrandom, $tEAMMrandom, $tEAMMsdrandom, $tSJAMMrandom, $tSJAMMsdrandom, $tSSJAMMrandom, $tSSJAMMsdrandom, $tSDMMrandom, $tSDMMsdrandom) = &topdMM($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
						push @topd_rand, $topdMMrandom;
						push @topd_unpruned_rand, $topdMMunprunedrandom;
						push @split_dist_rand, $SplitDistrandom;
						push @Q_rand, $qrandom;
						push @M_rand, $mrandom;

						push @A_num_sp_elim_random, $num_sp_elimMMrandom;
						push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1MMrandom;
						push @A_min_random, $minMMrandom;
						$A_DisTax_random{$DisTaxrandom}++;
						push @A_nquart_random, $nquartMMrandom;
						push @A_qsolved_random, $qsolvedMMrandom; 
						push @A_qdifferent_random, $qdifferentMMrandom; 
						push @A_qr2_random, $qr2MMrandom; 
						push @A_qr1_random, $qr1MMrandom; 
						push @A_qUnresolved_random, $qUnresolvedMMrandom; 
						push @A_qDC_random, $qDCMMrandom; 
						push @A_qEA_random, $qEAMMrandom; 
						push @A_qSJA_random, $qSJAMMrandom; 
						push @A_qSSJA_random, $qSSJAMMrandom;
						push @A_qSD_random, $qSDMMrandom; 
						push @A_ntripl_random, $ntriplMMrandom; 
						push @A_tsolved_random, $tsolvedMMrandom;
						push @A_tdifferent_random, $tdifferentMMrandom; 
						push @A_tr1_random, $tr1MMrandom; 
						push @A_tr2_random, $tr2MMrandom; 
						push @A_tUnresolved_random, $tUnresolvedMMrandom; 
						push @A_tDC_random, $tDCMMrandom; 
						push @A_tEA_random, $tEAMMrandom; 
						push @A_tSJA_random, $tSJAMMrandom; 
						push @A_tSSJA_random, $tSSJAMMrandom; 
						push @A_tSD_random, $tSDMMrandom;
					}
				}
			}
			
			($mean_rand,$sd_rand)=&meanSD(@topd_rand);
			($mean_unpruned_rand,$sd_unpruned_rand)=&meanSD(@topd_unpruned_rand);
			($mean_rand_split,$sd_rand_split)=&meanSD(@split_dist_rand);
			($mean_rand_q,$sd_rand_q)=&meanSD(@Q_rand);
			($mean_rand_m,$sd_rand_m)=&meanSD(@M_rand);

			($mean_num_sp_elim_random, $sd_num_sp_elim_random) = &meanSD(@A_num_sp_elim_random);
			($mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random) = &meanSD(@A_e_n_taxa_tree1_random);
			($mean_min_random, $sd_min_random) = &meanSD(@A_min_random);
			$dif_elemrandom = "";
			foreach $DisTaxSMelementrandom(sort keys(%A_DisTax_random)) {
				$aux_dist__random = $A_DisTax_random{$DisTaxSMelementrandom}; 
				$DisTaxSMelementrandom =~ s/-/ /g;
				$dif_elemrandom .= "[$DisTaxSMelementrandom("."$aux_dist__random".")] ";
			}
			($mean_nquart_random, $sd_nquart_random) = &meanSD(@A_nquart_random);
			($mean_qsolved_random, $sd_qsolved_random) = &meanSD(@A_qsolved_random);
			($mean_qdifferent_random, $sd_qdifferent_random) = &meanSD(@A_qdifferent_random);
			($mean_qr2_random, $sd_qr2_random) = &meanSD(@A_qr2_random);
			($mean_qr1_random, $sd_qr1_random) = &meanSD(@A_qr1_random);
			($mean_qUnresolved_random, $sd_qUnresolved_random) = &meanSD(@A_qUnresolved_random);
			($mean_qDC_random, $sd_qDC_random) = &meanSD(@A_qDC_random);
			($mean_qEA_random, $sd_qEA_random) = &meanSD(@A_qEA_random);
			($mean_qSJA_random, $sd_qSJA_random) = &meanSD(@A_qSJA_random);
			($mean_qSSJA_random, $sd_qSSJA_random) = &meanSD(@A_qSSJA_random);
			($mean_qSD_random, $sd_qSD_random) = &meanSD(@A_qSD_random);
			($mean_ntripl_random, $sd_ntripl_random) = &meanSD(@A_ntripl_random);
			($mean_tsolved_random, $sd_tsolved_random) = &meanSD(@A_tsolved_random);
			($mean_tdifferent_random, $sd_tdifferent_random) = &meanSD(@A_tdifferent_random);
			($mean_tr1_random, $sd_tr1_random) = &meanSD(@A_tr1_random);
			($mean_tr2_random, $sd_tr2_random) = &meanSD(@A_tr2_random);
			($mean_tUnresolved_random, $sd_tUnresolved_random) = &meanSD(@A_tUnresolved_random);
			($mean_tDC_random, $sd_tDC_random) = &meanSD(@A_tDC_random);
			($mean_tEA_random, $sd_tEA_random) = &meanSD(@A_tEA_random);
			($mean_tSJA_random, $sd_tSJA_random) = &meanSD(@A_tSJA_random);
			($mean_tSSJA_random, $sd_tSSJA_random) = &meanSD(@A_tSSJA_random);
			($mean_tSD_random, $sd_tSD_random) = &meanSD(@A_tSD_random);
		}

		
		print "\n############################################# RESULTS $name1_IN - $name2_IN ##########################################\n\n";
		if ($resp eq "S")  {
			printf "* Percentage of taxa in common: %6.1f", $pcCom2;
			print "%\n";
			printf OUTa "* Percentage of taxa in common: %6.1f", $pcCom2;
			print OUTa "%\n";
			$topdUnpruned_ori = $topd_ori * (1+(1-($pcCom2 /100)));

			if ( ($method eq 'all') || ($method eq 'nodal') ) {
				printf "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
				printf OUTa "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
			}
			if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
				print "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
				print OUTa "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
			}
			if ( ($method eq 'all') || ($method eq 'disagree') ) {
				$DisTax =~ s/-/ /g;
				print OUTa "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
				print "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
			}
			if ( ($method eq 'all') || ($method eq 'quartets') ) {
				printf "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
				printf OUTa "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
			}
			if ( ($method eq 'all') || ($method eq 'triplets') ) {
				printf "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
				printf OUTa "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
			}

			if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
				$mean_unpruned_rand = $mean_rand * (1+(1-($pcCom2 /100)));
				$sd_unpruned_rand = $sd_rand * (1+(1-($pcCom2 /100)));
				if ( ($method eq 'all') || ($method eq 'nodal') ) {
					printf "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
					printf OUTa "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
				}
				if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
					printf "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
					printf OUTa "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
				}
				if ( ($method eq 'all') || ($method eq 'disagree') ) {
					printf "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
					printf OUTa "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
				}
				if ( ($method eq 'all') || ($method eq 'quartets') ) {
					printf "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random;
					printf OUTa "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random;
				}
				if ( ($method eq 'all') || ($method eq 'triplets') ) {
					printf "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
					printf OUTa "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n",$mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
				}
			}
			print "\n";
		}
		else {
			printf "* Percentage of taxa in common: %6.1f", $pcCom2;
			print "%\n";
			printf OUTa "* Percentage of taxa in common: %6.1f", $pcCom2;
			print OUTa"%\n";
			if ( ($family1 eq "S") && ($family2 eq "S") ) {
				
				if ( ($method eq 'all') || ($method eq 'nodal') ) {
					printf "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
					printf OUTa "* Nodal Distance(Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
				}
				if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
					print "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
					print OUTa "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
				}
				if ( ($method eq 'all') || ($method eq 'disagree') ) {
					$DisTax =~ s/-/ /g;
					print OUTa "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
					print "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
				}
				if ( ($method eq 'all') || ($method eq 'quartets') ) {
					printf "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
					printf OUTa "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
				}
				if ( ($method eq 'all') || ($method eq 'triplets') ) {
					printf "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
					printf OUTa "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
				}
			}
			elsif ( ($family1 eq "S") && ($family2 eq "M") ) {
				
				if ( ($method eq 'all') || ($method eq 'nodal') ) {
					printf "\n* Nodal Distance SM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;
					printf OUTa "* Nodal Distance SM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;
	
				}
				if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
					printf "* Split Distance SM [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
					printf OUTa "* Split Distance SM [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
				}
				if ( ($method eq 'all') || ($method eq 'disagree') ) {
					printf "* Disagreement SM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
					printf OUTa "* Disagreement SM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
				}
				if ( ($method eq 'all') || ($method eq 'quartets') ) {
					printf "* Quartets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
					printf OUTa "* Quartets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
				}
				if ( ($method eq 'all') || ($method eq 'triplets') ) {
					printf "* Triplets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
					printf OUTa "* Triplets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
				}
			}elsif ( ($family1 eq "M") && ($family2 eq "S") ) {
				
				if ( ($method eq 'all') || ($method eq 'nodal') ) {
					printf "\n* Nodal Distance MS (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;
					printf OUTa "* Nodal Distance MS (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;
				}
				if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
					printf "* Split Distance MS [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
					printf OUTa "* Split Distance MS [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
				}
				if ( ($method eq 'all') || ($method eq 'disagree') ) {
					printf "* Disagreement MS [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
					printf OUTa "* Disagreement MS [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
				}
				if ( ($method eq 'all') || ($method eq 'quartets') ) {
					printf "* Quartets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
					printf OUTa "* Quartets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
				}
				if ( ($method eq 'all') || ($method eq 'triplets') ) {
					printf "* Triplets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
					printf OUTa "* Triplets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
				}
			}elsif ( ($family1 eq "M") && ($family2 eq "M") ) {
				
				if ( ($method eq 'all') || ($method eq 'nodal') ) {
					printf "\n* Nodal Distance MM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdMM_ori, $SD_ori, $topdMMunpruned_ori, $SDunpruned_ori;
					printf OUTa "* Nodal Distance MM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdMM_ori, $SD_ori, $topdMMunpruned_ori, $SDunpruned_ori;
				}
				if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
					printf "* Split Distance MM [differents/possibles]: $SplitDistMM +/- %5.3f [ $qMM +/- %5.3f \/ $mMM +/- %5.3f ]\n", $SplitDistMMsd, $qMMsd, $mMMsd ;
					printf OUTa "* Split Distance MM [differents/possibles]: $SplitDistMM +/- %5.3f [ $qMM +/- %5.3f \/ $mMM +/- %5.3f ]\n", $SplitDistMMsd, $qMMsd, $mMMsd ;
				}
				if ( ($method eq 'all') || ($method eq 'disagree') ) {
					printf "* Disagreement MM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimMM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1MMsd, $minMM, $minMMsd ;
					printf OUTa "* Disagreement MM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimMM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1MMsd, $minMM, $minMMsd ;
				}
				if ( ($method eq 'all') || ($method eq 'quartets') ) {
					printf "* Quartets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartMM, $nquartMMsd, $qsolvedMM, $qsolvedMMsd, $qdifferentMM, $qdifferentMMsd, $qr2MM, $qr2MMsd, $qr1MM, $qr1MMsd, $qUnresolvedMM, $qUnresolvedMMsd, $qDCMM, $qDCMMsd, $qEAMM, $qEAMMsd, $qSJAMM, $qSJAMMsd, $qSSJAMM, $qSSJAMMsd, $qSDMM, $qSDMMsd;
					printf OUTa "* Quartets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartMM, $nquartMMsd, $qsolvedMM, $qsolvedMMsd, $qdifferentMM, $qdifferentMMsd, $qr2MM, $qr2MMsd, $qr1MM, $qr1MMsd, $qUnresolvedMM, $qUnresolvedMMsd, $qDCMM, $qDCMMsd, $qEAMM, $qEAMMsd, $qSJAMM, $qSJAMMsd, $qSSJAMM, $qSSJAMMsd, $qSDMM, $qSDMMsd;
				}
				if ( ($method eq 'all') || ($method eq 'triplets') ) {
					printf "* Triplets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplMM, $ntriplMMsd, $tsolvedMM, $tsolvedMMsd, $tdifferentMM, $tdifferentMMsd,$tr2MM, $tr2MMsd, $tr1MM, $tr1MMsd, $tUnresolvedMM, $tUnresolvedMMsd, $tDCMM, $tDCMMsd, $tEAMM, $tEAMMsd, $tSJAMM, $tSJAMMsd, $tSSJAMM, $tSSJAMMsd, $tSDMM, $tSDMMsd;
					printf OUTa "* Triplets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplMM, $ntriplMMsd, $tsolvedMM, $tsolvedMMsd, $tdifferentMM, $tdifferentMMsd,$tr2MM, $tr2MMsd, $tr1MM, $tr1MMsd, $tUnresolvedMM, $tUnresolvedMMsd, $tDCMM, $tDCMMsd, $tEAMM, $tEAMMsd, $tSJAMM, $tSJAMMsd, $tSSJAMM, $tSSJAMMsd, $tSDMM, $tSDMMsd;
				}
			}
			if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
				if ( ($method eq 'all') || ($method eq 'nodal') ) {
					printf "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
					printf OUTa "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
				}
				if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
					printf "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
					printf OUTa "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
				}
				if ( ($method eq 'all') || ($method eq 'disagree') ) {
					printf "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
					printf OUTa "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
				}
				if ( ($method eq 'all') || ($method eq 'quartets') ) {
					printf "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random;
					printf OUTa "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random;
				}
				if ( ($method eq 'all') || ($method eq 'triplets') ) {
					printf "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
					printf OUTa "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n",$mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
				}
			}
			print "\n";
		}
	}else {
		print "The overlap in the data between $name1_IN and $name2_IN is too separate. ($LocTrees) \n";
		print OUTa "The overlap in the data between $name1_IN and $name2_IN is too separate. ($LocTrees) \n";
		if ( ($method eq 'all') || ($method eq 'nodal') ) {
			printf "* Nodal Distance (Pruned/Unpruned): 0.000000 / 0.000000\n";
			printf OUTa "* Nodal Distance (Pruned/Unpruned): 0.000000 / 0.000000\n";
		}
		if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
			print "* Split Distance [differents/possibles]: 0 [ 0 \/ 0 ]\n";
			print OUTa "* Split Distance [differents/possibles]: 0 [ 0 \/ 0 ]\n";
		}
		if ( ($method eq 'all') || ($method eq 'disagree') ) {
			print "* Disagreement [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
			print OUTa "* Disagreement [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
		}
		if ( ($method eq 'all') || ($method eq 'quartets') ) {
			print "* Quartets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
			print OUTa "* Quartets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
		}
		if ( ($method eq 'all') || ($method eq 'triplets') ) {
			print "* Triplets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
			print OUTa "* Triplets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
		}

		if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
			if ( ($method eq 'all') || ($method eq 'nodal') ) {
				print "* Nodal Distance random (Pruned/Unpruned): ( 0.000000 +/- 0.000000 ) / ( 0.000000 +/- 0.000000 )\n";
				print OUTa "* Nodal Distance random (Pruned/Unpruned): ( 0.000000 +/- 0.000000 ) / ( 0.000 +/- 0.000000 )\n";
			}
			if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
				print "* Split Distance random [differents/possibles]: 0 [ 0 \/ 0 ]\n";
				print OUTa "* Split Distance random [differents/possibles]: 0 [ 0 \/ 0 ]\n";
			}
			if ( ($method eq 'all') || ($method eq 'disagree') ) {
				print "* Disagreement random [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
				print OUTa "* Disagreement random [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
			}
			if ( ($method eq 'all') || ($method eq 'quartets') ) {
				print "* Quartets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
				print OUTa "* Quartets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
			}
			if ( ($method eq 'all') || ($method eq 'triplets') ) {
				print "* Triplets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
				print OUTa "* Triplets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
			}
		}
	}
	close OUTa;
	print "\nThanks for using TOPD-fMtS .................................. Bye.\n\n";
}

sub menu2() {
	($LocTrees, $randomTrees, $nombreRandomTrees, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $fileoutput,$units) = @_;
	
	if ($fileoutput eq "") {	
		open (OUTb, ">topd_result_multiple.txt");
	}else {
		open (OUTb, ">$fileoutput");	
	}

	
	my %name_tree = &getTrees($LocTrees);
	my %fets = "";
	foreach $name1_IN(sort keys(%name_tree)) {
		$tree1_IN = $name_tree{$name1_IN};
		foreach $name2_IN(sort keys(%name_tree)) {
			$tree2_IN = $name_tree{$name2_IN};
			if ( ($name1_IN ne $name2_IN) && ($fets{$name1_IN}{$name2_IN} eq "") ) {
				$fets{$name1_IN}{$name2_IN} = "s";
				$fets{$name2_IN}{$name1_IN} = "s";
				print OUTb "##################################### topd $name1_IN - $name2_IN #######################################\n";

				
				my ($Nsp1, $family1, $Nsp2, $family2) = &MoS($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
				
				if ( ($family1 eq "S") && ($family2 eq "S") ) {
					
					($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $tree1_IN, $name2_IN, $tree2_IN, $print_prev, $method, $level, $units);
			
				}
				elsif ( ($family1 eq "S") && ($family2 eq "M") ) {
					
					
					($resp, $ptree1_in, $ptree2_in) = &checkPruneMGF($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
			
					
					if ($resp eq "S") {
						if ($print_prev eq "n") {
						}else {
							print "\nPrunned trees are both single gene trees.\n";
						}
						($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $ptree1_in, $name2_IN, $ptree2_in, $print_prev, $method, $level, $units);
						($com, $pcCom2) = &pcCom($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
					}else {
						($topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori, $com, $pcCom2, $SplitDistSM, $SplitDistSMsd, $qSM, $qSMsd, $mSM, $mSMsd, $dif_elem, $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd, $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd, $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd) = &topdSM($name1_IN, $tree1_IN, $name2_IN, $tree2_IN, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
					}
				}elsif ( ($family1 eq "M") && ($family2 eq "S") ) {
					
					
					($resp, $ptree2_in, $ptree1_in) = &checkPruneMGF($name2_IN, $tree2_IN, $name1_IN, $tree1_IN);
					
					if ($resp eq "S") {
						if ($print_prev eq "n") {
						}else {
							print "\nPrunned trees are both single gene trees.\n";
						}
						($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $ptree1_in, $name2_IN, $ptree2_in, $print_prev, $method, $level, $units);
						($com, $pcCom2) = &pcCom($name2_IN, $tree2_IN, $name1_IN, $tree1_IN);
					}else {
						($topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori, $com, $pcCom2, $SplitDistSM, $SplitDistSMsd, $qSM, $qSMsd, $mSM, $mSMsd, $dif_elem, $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd, $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd, $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd) = &topdSM($name2_IN, $tree2_IN, $name1_IN, $tree1_IN, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
					}
				}elsif ( ($family1 eq "M") && ($family2 eq "M") ) {
					
					
					($resp, $ptree1_in, $ptree2_in) = &checkPruneMGF($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
			
					
					if ($resp eq "S") {
						if ($print_prev eq "n") {
						}else {
							print "\nPrunned trees are both single gene trees.\n";
						}
						($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $ptree1_in, $name2_IN, $ptree2_in, $print_prev, $method, $level, $units);
						($com, $pcCom2) = &pcCom($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
					}else {
						($topdMM_ori, $SD_ori, $topdMMunpruned_ori, $SDunpruned_ori, $com, $pcCom2, $SplitDistMM, $SplitDistMMsd, $qMM, $qMMsd, $mMM, $mMMsd, $dif_elem, $num_sp_elimMM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1MMsd, $minMM, $minMMsd, $nquartMM, $nquartMMsd, $qsolvedMM, $qsolvedMMsd, $qdifferentMM, $qdifferentMMsd, $qr2MM, $qr2MMsd, $qr1MM, $qr1MMsd, $qUnresolvedMM, $qUnresolvedMMsd, $qDCMM, $qDCMMsd, $qEAMM, $qEAMMsd, $qSJAMM, $qSJAMMsd, $qSSJAMM, $qSSJAMMsd, $qSDMM, $qSDMMsd, $ntriplMM, $ntriplMMsd, $tsolvedMM, $tsolvedMMsd, $tdifferentMM, $tdifferentMMsd,$tr2MM, $tr2MMsd, $tr1MM, $tr1MMsd, $tUnresolvedMM, $tUnresolvedMMsd, $tDCMM, $tDCMMsd, $tEAMM, $tEAMMsd, $tSJAMM, $tSJAMMsd, $tSSJAMM, $tSSJAMMsd, $tSDMM, $tSDMMsd) = &topdMM($name1_IN, $tree1_IN, $name2_IN, $tree2_IN, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
					}
				}
				if ($com > 3) {
					if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
						
					
						my @topd_rand =();
						my @topd_unpruned_rand = ();
						my @split_dist_rand = ();
						my @Q_rand = ();
						my @M_rand = ();

						my @A_num_sp_elim_random = ();
						my @A_e_n_taxa_tree1_random = (); 
						my @A_min_random = (); 
						my %A_DisTax_random = ();
						my @A_nquart_random = ();
						my @A_qsolved_random = (); 
						my @A_qdifferent_random = (); 
						my @A_qr2_random = (); 
						my @A_qr1_random = (); 
						my @A_qUnresolved_random = (); 
						my @A_qDC_random = (); 
						my @A_qEA_random = (); 
						my @A_qSJA_random = (); 
						my @A_qSSJA_random = ();
						my @A_qSD_random = (); 
						my @A_ntripl_random = (); 
						my @A_tsolved_random = ();
						my @A_tdifferent_random = (); 
						my @A_tr1_random = (); 
						my @A_tr2_random = (); 
						my @A_tUnresolved_random = (); 
						my @A_tDC_random = (); 
						my @A_tEA_random = (); 
						my @A_tSJA_random = (); 
						my @A_tSSJA_random = (); 
						my @A_tSD_random = ();

						for ($ir=0;$ir<$nombreRandomTrees;$ir++) {
							if ($print_prev eq "n") {
							}else {
								print "\n### RANDOM $ir ###";
							}
							my ($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom) = ();
							if ($randomTrees eq 'guided') {
								($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom) = &random($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
							} elsif ($randomTrees eq 'random') {
								($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom) = &random2($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
							}

							
							my ($Nsp1random, $family1random, $Nsp2random, $family2random) = &MoS($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom);
							if ( ($family1random eq "S") && ($family2random eq "S") ) {
								
								my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom, $print_prev, $method, $level, $units);
								push @topd_rand, $topdrandom;
								push @topd_unpruned_rand, $topdUnprunedrandom;
								push @split_dist_rand, $SplitDistrandom;
								push @Q_rand, $qrandom;
								push @M_rand, $mrandom;

								push @A_num_sp_elim_random, $num_sp_elimrandom;
								push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
								push @A_min_random, $minrandom;
								$A_DisTax_random{$DisTaxrandom}++;
								push @A_nquart_random, $nquartrandom;
								push @A_qsolved_random, $qsolvedrandom; 
								push @A_qdifferent_random, $qdifferentrandom; 
								push @A_qr2_random, $qr2random; 
								push @A_qr1_random, $qr1random; 
								push @A_qUnresolved_random, $qUnresolvedrandom; 
								push @A_qDC_random, $qDCrandom; 
								push @A_qEA_random, $qEArandom; 
								push @A_qSJA_random, $qSJArandom; 
								push @A_qSSJA_random, $qSSJArandom;
								push @A_qSD_random, $qSDrandom; 
								push @A_ntripl_random, $ntriplrandom; 
								push @A_tsolved_random, $tsolvedrandom;
								push @A_tdifferent_random, $tdifferentrandom; 
								push @A_tr1_random, $tr1random; 
								push @A_tr2_random, $tr2random; 
								push @A_tUnresolved_random, $tUnresolvedrandom; 
								push @A_tDC_random, $tDCrandom; 
								push @A_tEA_random, $tEArandom; 
								push @A_tSJA_random, $tSJArandom; 
								push @A_tSSJA_random, $tSSJArandom; 
								push @A_tSD_random, $tSDrandom;
							}
							elsif ( ($family1random eq "S") && ($family2random eq "M") ) {
								
							
								($resp, $ptree1_in_r, $ptree2_in_r) = &checkPruneMGF($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom);
				
								
								if ($resp eq "S") {
									if ($print_prev eq "n") {
									}else {
										print "\nPrunned trees are both single gene trees.\n";
									}
									my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $ptree1_in_r, $name2_INrandom, $ptree2_in_r, $print_prev, $method, $level, $units);
									push @topd_rand, $topdrandom;
									push @topd_unpruned_rand, $topdUnprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;

									push @A_num_sp_elim_random, $num_sp_elimrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
									push @A_min_random, $minrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartrandom;
									push @A_qsolved_random, $qsolvedrandom; 
									push @A_qdifferent_random, $qdifferentrandom; 
									push @A_qr2_random, $qr2random; 
									push @A_qr1_random, $qr1random; 
									push @A_qUnresolved_random, $qUnresolvedrandom; 
									push @A_qDC_random, $qDCrandom; 
									push @A_qEA_random, $qEArandom; 
									push @A_qSJA_random, $qSJArandom; 
									push @A_qSSJA_random, $qSSJArandom;
									push @A_qSD_random, $qSDrandom; 
									push @A_ntripl_random, $ntriplrandom; 
									push @A_tsolved_random, $tsolvedrandom;
									push @A_tdifferent_random, $tdifferentrandom; 
									push @A_tr1_random, $tr1random; 
									push @A_tr2_random, $tr2random; 
									push @A_tUnresolved_random, $tUnresolvedrandom; 
									push @A_tDC_random, $tDCrandom; 
									push @A_tEA_random, $tEArandom; 
									push @A_tSJA_random, $tSJArandom; 
									push @A_tSSJA_random, $tSSJArandom; 
									push @A_tSD_random, $tSDrandom;
								}else {
									($topdSMrandom, $SDrandom, $topdSMunprunedrandom, $SDunprunedrandom, $com, $pcCom2, $SplitDistrandom, $SplitDistsdrandom, $qrandom, $qsdrandom, $mrandom, $msdrandom, $dif_elemrandom, $num_sp_elimSMrandom, $num_sp_elimSMsdrandom, $e_n_taxa_tree1SMrandom, $e_n_taxa_tree1SMsdrandom, $minSMrandom, $minSMsdrandom, $nquartSMrandom, $nquartSMsdrandom, $qsolvedSMrandom, $qsolvedSMsdrandom, $qdifferentSMrandom, $qdifferentSMsdrandom, $qr2SMrandom, $qr2SMsdrandom, $qr1SMrandom, $qr1SMsdrandom, $qUnresolvedSMrandom, $qUnresolvedSMsdrandom, $qDCSMrandom, $qDCSMsdrandom, $qEASMrandom, $qEASMsdrandom, $qSJASMrandom, $qSJASMsdrandom, $qSSJASMrandom, $qSSJASMsdrandom, $qSDSMrandom, $qSDSMsdrandom, $ntriplSMrandom, $ntriplSMsdrandom, $tsolvedSMrandom, $tsolvedSMsdrandom, $tdifferentSMrandom, $tdifferentSMsdrandom, $tr2SMrandom, $tr2SMsdrandom, $tr1SMrandom, $tr1SMsdrandom, $tUnresolvedSMrandom, $tUnresolvedSMsdrandom, $tDCSMrandom, $tDCSMsdrandom, $tEASMrandom, $tEASMsdrandom, $tSJASMrandom, $tSJASMsdrandom, $tSSJASMrandom, $tSSJASMsdrandom, $tSDSMrandom, $tSDSMsdrandom) = &topdSM($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
									push @topd_rand, $topdSMrandom;
									push @topd_unpruned_rand, $topdSMunprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;

									push @A_num_sp_elim_random, $num_sp_SMelimrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1SMrandom;
									push @A_min_random, $minSMrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartSMrandom;
									push @A_qsolved_random, $qsolvedSMrandom; 
									push @A_qdifferent_random, $qdifferentSMrandom; 
									push @A_qr2_random, $qr2SMrandom; 
									push @A_qr1_random, $qr1SMrandom; 
									push @A_qUnresolved_random, $qUnresolvedSMrandom; 
									push @A_qDC_random, $qDCSMrandom; 
									push @A_qEA_random, $qEASMrandom; 
									push @A_qSJA_random, $qSJASMrandom; 
									push @A_qSSJA_random, $qSSJASMrandom;
									push @A_qSD_random, $qSDSMrandom; 
									push @A_ntripl_random, $ntriplSMrandom; 
									push @A_tsolved_random, $tsolvedSMrandom;
									push @A_tdifferent_random, $tdifferentSMrandom; 
									push @A_tr1_random, $tr1SMrandom; 
									push @A_tr2_random, $tr2SMrandom; 
									push @A_tUnresolved_random, $tUnresolvedSMrandom; 
									push @A_tDC_random, $tDCSMrandom; 
									push @A_tEA_random, $tEASMrandom; 
									push @A_tSJA_random, $tSJASMrandom; 
									push @A_tSSJA_random, $tSSJASMrandom; 
									push @A_tSD_random, $tSDSMrandom;

								}	
							}elsif ( ($family1random eq "M") && ($family2random eq "S") ) {
								
				
								
								($resp, $ptree1_in_r, $ptree2_in_r) = &checkPruneMGF($name2_INrandom, $tree2_INrandom, $name1_INrandom, $tree1_INrandom);
								
								
								if ($resp eq "S") {
									if ($print_prev eq "n") {
									}else {
										print "\nPrunned trees are both single gene trees.\n";
									}
									my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $ptree1_in_r, $name2_INrandom, $ptree2_in_r, $print_prev, $method, $level, $units);
									push @topd_rand, $topdrandom;
									push @topd_unpruned_rand, $topdUnprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;

									push @A_num_sp_elim_random, $num_sp_elimrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
									push @A_min_random, $minrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartrandom;
									push @A_qsolved_random, $qsolvedrandom; 
									push @A_qdifferent_random, $qdifferentrandom; 
									push @A_qr2_random, $qr2random; 
									push @A_qr1_random, $qr1random; 
									push @A_qUnresolved_random, $qUnresolvedrandom; 
									push @A_qDC_random, $qDCrandom; 
									push @A_qEA_random, $qEArandom; 
									push @A_qSJA_random, $qSJArandom; 
									push @A_qSSJA_random, $qSSJArandom;
									push @A_qSD_random, $qSDrandom; 
									push @A_ntripl_random, $ntriplrandom; 
									push @A_tsolved_random, $tsolvedrandom;
									push @A_tdifferent_random, $tdifferentrandom; 
									push @A_tr1_random, $tr1random; 
									push @A_tr2_random, $tr2random; 
									push @A_tUnresolved_random, $tUnresolvedrandom; 
									push @A_tDC_random, $tDCrandom; 
									push @A_tEA_random, $tEArandom; 
									push @A_tSJA_random, $tSJArandom; 
									push @A_tSSJA_random, $tSSJArandom; 
									push @A_tSD_random, $tSDrandom;
								}else {
									($topdSMrandom, $SDrandom, $topdSMunprunedrandom, $SDunprunedrandom, $com, $pcCom2, $SplitDistrandom, $SplitDistsdrandom, $qrandom, $qsdrandom, $mrandom, $msdrandom, $DisTaxrandom, $num_sp_elimSMrandom, $num_sp_elimSMsdrandom, $e_n_taxa_tree1SMrandom, $e_n_taxa_tree1SMsdrandom, $minSMrandom, $minSMsdrandom, $nquartSMrandom, $nquartSMsdrandom, $qsolvedSMrandom, $qsolvedSMsdrandom, $qdifferentSMrandom, $qdifferentSMsdrandom, $qr2SMrandom, $qr2SMsdrandom, $qr1SMrandom, $qr1SMsdrandom, $qUnresolvedSMrandom, $qUnresolvedSMsdrandom, $qDCSMrandom, $qDCSMsdrandom, $qEASMrandom, $qEASMsdrandom, $qSJASMrandom, $qSJASMsdrandom, $qSSJASMrandom, $qSSJASMsdrandom, $qSDSMrandom, $qSDSMsdrandom, $ntriplSMrandom, $ntriplSMsdrandom, $tsolvedSMrandom, $tsolvedSMsdrandom, $tdifferentSMrandom, $tdifferentSMsdrandom, $tr2SMrandom, $tr2SMsdrandom, $tr1SMrandom, $tr1SMsdrandom, $tUnresolvedSMrandom, $tUnresolvedSMsdrandom, $tDCSMrandom, $tDCSMsdrandom, $tEASMrandom, $tEASMsdrandom, $tSJASMrandom, $tSJASMsdrandom, $tSSJASMrandom, $tSSJASMsdrandom, $tSDSMrandom, $tSDSMsdrandom) = &topdSM($name2_INrandom, $tree2_INrandom, $name1_INrandom, $tree1_INrandom, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
									push @topd_rand, $topdSMrandom;
									push @topd_unpruned_rand, $topdSMunprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;

									push @A_num_sp_elim_random, $num_sp_elimSMrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1SMrandom;
									push @A_min_random, $minSMrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartSMrandom;
									push @A_qsolved_random, $qsolvedSMrandom; 
									push @A_qdifferent_random, $qdifferentSMrandom; 
									push @A_qr2_random, $qr2SMrandom; 
									push @A_qr1_random, $qr1SMrandom; 
									push @A_qUnresolved_random, $qUnresolvedSMrandom; 
									push @A_qDC_random, $qDCSMrandom; 
									push @A_qEA_random, $qEASMrandom; 
									push @A_qSJA_random, $qSJASMrandom; 
									push @A_qSSJA_random, $qSSJASMrandom;
									push @A_qSD_random, $qSDSMrandom; 
									push @A_ntripl_random, $ntriplSMrandom; 
									push @A_tsolved_random, $tsolvedSMrandom;
									push @A_tdifferent_random, $tdifferentSMrandom; 
									push @A_tr1_random, $tr1SMrandom; 
									push @A_tr2_random, $tr2SMrandom; 
									push @A_tUnresolved_random, $tUnresolvedSMrandom; 
									push @A_tDC_random, $tDCSMrandom; 
									push @A_tEA_random, $tEASMrandom; 
									push @A_tSJA_random, $tSJASMrandom; 
									push @A_tSSJA_random, $tSSJASMrandom; 
									push @A_tSD_random, $tSDSMrandom;
								}	
							}elsif ( ($family1random eq "M") && ($family2random eq "M") ) {
								
								
								($resp, $ptree1_in_r, $ptree2_in_r) = &checkPruneMGF($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom);
								
								if ($resp eq "S") {
									if ($print_prev eq "n") {
									}else {
										print "\nPrunned trees are both single gene trees.\n";
									}
									my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $ptree1_in_r, $name2_INrandom, $ptree2_in_r, $print_prev, $method, $level, $units);
									push @topd_rand, $topdrandom;
									push @topd_unpruned_rand, $topdUnprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;

									push @A_num_sp_elim_random, $num_sp_elimrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
									push @A_min_random, $minrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartrandom;
									push @A_qsolved_random, $qsolvedrandom; 
									push @A_qdifferent_random, $qdifferentrandom; 
									push @A_qr2_random, $qr2random; 
									push @A_qr1_random, $qr1random; 
									push @A_qUnresolved_random, $qUnresolvedrandom; 
									push @A_qDC_random, $qDCrandom; 
									push @A_qEA_random, $qEArandom; 
									push @A_qSJA_random, $qSJArandom; 
									push @A_qSSJA_random, $qSSJArandom;
									push @A_qSD_random, $qSDrandom; 
									push @A_ntripl_random, $ntriplrandom; 
									push @A_tsolved_random, $tsolvedrandom;
									push @A_tdifferent_random, $tdifferentrandom; 
									push @A_tr1_random, $tr1random; 
									push @A_tr2_random, $tr2random; 
									push @A_tUnresolved_random, $tUnresolvedrandom; 
									push @A_tDC_random, $tDCrandom; 
									push @A_tEA_random, $tEArandom; 
									push @A_tSJA_random, $tSJArandom; 
									push @A_tSSJA_random, $tSSJArandom; 
									push @A_tSD_random, $tSDrandom;
								}else {
									($topdMMrandom, $SDrandom, $topdMMunprunedrandom, $SDunprunedrandom, $com, $pcCom2, $SplitDistrandom, $SplitDistsdrandom, $qrandom, $qsdrandom, $mrandom, $msdrandom, $DisTaxrandom, $num_sp_elimMMrandom, $num_sp_elimMMsdrandom, $e_n_taxa_tree1MMrandom, $e_n_taxa_tree1MMsdrandom, $minMMrandom, $minMMsdrandom, $nquartMMrandom, $nquartMMsdrandom, $qsolvedMMrandom, $qsolvedMMsdrandom, $qdifferentMMrandom, $qdifferentMMsdrandom, $qr2MMrandom, $qr2MMsdrandom, $qr1MMrandom, $qr1MMsdrandom, $qUnresolvedMMrandom, $qUnresolvedMMsdrandom, $qDCMMrandom, $qDCMMsdrandom, $qEAMMrandom, $qEAMMsdrandom, $qSJAMMrandom, $qSJAMMsdrandom, $qSSJAMMrandom, $qSSJAMMsdrandom, $qSDMMrandom, $qSDMMsdrandom, $ntriplMMrandom, $ntriplMMsdrandom, $tsolvedMMrandom, $tsolvedMMsdrandom, $tdifferentMMrandom, $tdifferentMMsdrandom, $tr2MMrandom, $tr2MMsdrandom, $tr1MMrandom, $tr1MMsdrandom, $tUnresolvedMMrandom, $tUnresolvedMMsdrandom, $tDCMMrandom, $tDCMMsdrandom, $tEAMMrandom, $tEAMMsdrandom, $tSJAMMrandom, $tSJAMMsdrandom, $tSSJAMMrandom, $tSSJAMMsdrandom, $tSDMMrandom, $tSDMMsdrandom) = &topdMM($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
									push @topd_rand, $topdMMrandom;
									push @topd_unpruned_rand, $topdMMunprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;

									push @A_num_sp_elim_random, $num_sp_elimMMrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1MMrandom;
									push @A_min_random, $minMMrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartMMrandom;
									push @A_qsolved_random, $qsolvedMMrandom; 
									push @A_qdifferent_random, $qdifferentMMrandom; 
									push @A_qr2_random, $qr2MMrandom; 
									push @A_qr1_random, $qr1MMrandom; 
									push @A_qUnresolved_random, $qUnresolvedMMrandom; 
									push @A_qDC_random, $qDCMMrandom; 
									push @A_qEA_random, $qEAMMrandom; 
									push @A_qSJA_random, $qSJAMMrandom; 
									push @A_qSSJA_random, $qSSJAMMrandom;
									push @A_qSD_random, $qSDMMrandom; 
									push @A_ntripl_random, $ntriplMMrandom; 
									push @A_tsolved_random, $tsolvedMMrandom;
									push @A_tdifferent_random, $tdifferentMMrandom; 
									push @A_tr1_random, $tr1MMrandom; 
									push @A_tr2_random, $tr2MMrandom; 
									push @A_tUnresolved_random, $tUnresolvedMMrandom; 
									push @A_tDC_random, $tDCMMrandom; 
									push @A_tEA_random, $tEAMMrandom; 
									push @A_tSJA_random, $tSJAMMrandom; 
									push @A_tSSJA_random, $tSSJAMMrandom; 
									push @A_tSD_random, $tSDMMrandom;
								}
							}
						}
						
						($mean_rand,$sd_rand)=&meanSD(@topd_rand);
						($mean_unpruned_rand,$sd_unpruned_rand)=&meanSD(@topd_unpruned_rand);
						($mean_rand_split,$sd_rand_split)=&meanSD(@split_dist_rand);
						($mean_rand_q,$sd_rand_q)=&meanSD(@Q_rand);
						($mean_rand_m,$sd_rand_m)=&meanSD(@M_rand);

						($mean_num_sp_elim_random, $sd_num_sp_elim_random) = &meanSD(@A_num_sp_elim_random);
						($mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random) = &meanSD(@A_e_n_taxa_tree1_random);
						($mean_min_random, $sd_min_random) = &meanSD(@A_min_random);
						$dif_elemrandom = "";
						foreach $DisTaxSMelementrandom(sort keys(%A_DisTax_random)) {
							$aux_dist__random = $A_DisTax_random{$DisTaxSMelementrandom}; 
							$DisTaxSMelementrandom =~ s/-/ /g;
							$dif_elemrandom .= "[$DisTaxSMelementrandom("."$aux_dist__random".")] ";
						}
						($mean_nquart_random, $sd_nquart_random) = &meanSD(@A_nquart_random);
						($mean_qsolved_random, $sd_qsolved_random) = &meanSD(@A_qsolved_random);
						($mean_qdifferent_random, $sd_qdifferent_random) = &meanSD(@A_qdifferent_random);
						($mean_qr2_random, $sd_qr2_random) = &meanSD(@A_qr2_random);
						($mean_qr1_random, $sd_qr1_random) = &meanSD(@A_qr1_random);
						($mean_qUnresolved_random, $sd_qUnresolved_random) = &meanSD(@A_qUnresolved_random);
						($mean_qDC_random, $sd_qDC_random) = &meanSD(@A_qDC_random);
						($mean_qEA_random, $sd_qEA_random) = &meanSD(@A_qEA_random);
						($mean_qSJA_random, $sd_qSJA_random) = &meanSD(@A_qSJA_random);
						($mean_qSSJA_random, $sd_qSSJA_random) = &meanSD(@A_qSSJA_random);
						($mean_qSD_random, $sd_qSD_random) = &meanSD(@A_qSD_random);
						($mean_ntripl_random, $sd_ntripl_random) = &meanSD(@A_ntripl_random);
						($mean_tsolved_random, $sd_tsolved_random) = &meanSD(@A_tsolved_random);
						($mean_tdifferent_random, $sd_tdifferent_random) = &meanSD(@A_tdifferent_random);
						($mean_tr1_random, $sd_tr1_random) = &meanSD(@A_tr1_random);
						($mean_tr2_random, $sd_tr2_random) = &meanSD(@A_tr2_random);
						($mean_tUnresolved_random, $sd_tUnresolved_random) = &meanSD(@A_tUnresolved_random);
						($mean_tDC_random, $sd_tDC_random) = &meanSD(@A_tDC_random);
						($mean_tEA_random, $sd_tEA_random) = &meanSD(@A_tEA_random);
						($mean_tSJA_random, $sd_tSJA_random) = &meanSD(@A_tSJA_random);
						($mean_tSSJA_random, $sd_tSSJA_random) = &meanSD(@A_tSSJA_random);
						($mean_tSD_random, $sd_tSD_random) = &meanSD(@A_tSD_random);
					}
			
					
					print "\n############################################ RESULTS $name1_IN - $name2_IN ###################################\n\n";
					if ($resp eq "S")  {
						printf "* Percentage of taxa in common: %6.1f", $pcCom2;
						print "%\n";
						printf OUTb "* Percentage of taxa in common: %6.1f", $pcCom2;
						print OUTb "%\n";
						$topdUnpruned_ori = $topd_ori * (1+(1-($pcCom2 /100)));

						if ( ($method eq 'all') || ($method eq 'nodal') ) {
							printf "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
							printf OUTb "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
						}
						if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
							print "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
							print OUTb "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
						}
						if ( ($method eq 'all') || ($method eq 'disagree') ) {
							$DisTax =~ s/-/ /g;
							print OUTb "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
							print "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
						}
						if ( ($method eq 'all') || ($method eq 'quartets') ) {
							printf "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
							printf OUTb "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
						}
						if ( ($method eq 'all') || ($method eq 'triplets') ) {
							printf "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
							printf OUTb "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
						}

						if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
							$mean_unpruned_rand = $mean_rand * (1+(1-($pcCom2 /100)));
							$sd_unpruned_rand = $sd_rand * (1+(1-($pcCom2 /100)));
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
								printf OUTb "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								printf "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
								printf OUTb "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								printf "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
								printf OUTb "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
							}
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random;
								printf OUTb "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random; 
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
								printf OUTb "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n",$mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
							}
						}
						print "\n";
					}
					else {
						printf "* Percentage of taxa in common: %6.1f", $pcCom2;
						print "%\n";
						printf OUTb "* Percentage of taxa in common: %6.1f", $pcCom2;
						print OUTb "%\n";
						if ( ($family1 eq "S") && ($family2 eq "S") ) {
							
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
								printf OUTb "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								print "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
								print OUTb "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								$DisTax =~ s/-/ /g;
								print OUTb "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
								print "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
							} 
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
								printf OUTb "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
								printf OUTb "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
							}
						}
						elsif ( ($family1 eq "S") && ($family2 eq "M") ) {
							
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "\n* Nodal Distance SM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;
								printf OUTb "* Nodal Distance SM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;
							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								printf "* Split Distance SM [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
								printf OUTb "* Split Distance SM [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								printf "* Disagreement SM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
								printf OUTb "* Disagreement SM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ; 
							}
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
								printf OUTb "* Quartets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
								printf OUTb "* Triplets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
							}
						}elsif ( ($family1 eq "M") && ($family2 eq "S") ) {
							
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "\n* Nodal Distance MS (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;
								printf OUTb "* Nodal Distance MS (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;
							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								printf "* Split Distance MS [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
								printf OUTb "* Split Distance MS [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								printf "* Disagreement MS [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
								printf OUTb "* Disagreement MS [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
							}
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
								printf OUTb "* Quartets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd; 
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
								printf OUTb "* Triplets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
							}
						}elsif ( ($family1 eq "M") && ($family2 eq "M") ) {
							
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "* Nodal Distance MM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdMM_ori, $SD_ori, $topdMMunpruned_ori, $SDunpruned_ori;
								printf OUTb "* Nodal Distance MM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdMM_ori, $SD_ori, $topdMMunpruned_ori, $SDunpruned_ori;
							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								printf "* Split Distance MM [differents/possibles]: $SplitDistMM +/- %5.3f [ $qMM +/- %5.3f \/ $mMM +/- %5.3f ]\n", $SplitDistMMsd, $qMMsd, $mMMsd ;
								printf OUTb "* Split Distance MM [differents/possibles]: $SplitDistMM +/- %5.3f [ $qMM +/- %5.3f \/ $mMM +/- %5.3f ]\n", $SplitDistMMsd, $qMMsd, $mMMsd ;
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								printf "* Disagreement MM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimMM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1MMsd, $minMM, $minMMsd ;
								printf OUTb "* Disagreement MM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimMM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1MMsd, $minMM, $minMMsd ;
							}
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartMM, $nquartMMsd, $qsolvedMM, $qsolvedMMsd, $qdifferentMM, $qdifferentMMsd, $qr2MM, $qr2MMsd, $qr1MM, $qr1MMsd, $qUnresolvedMM, $qUnresolvedMMsd, $qDCMM, $qDCMMsd, $qEAMM, $qEAMMsd, $qSJAMM, $qSJAMMsd, $qSSJAMM, $qSSJAMMsd, $qSDMM, $qSDMMsd;
								printf OUTb "* Quartets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartMM, $nquartMMsd, $qsolvedMM, $qsolvedMMsd, $qdifferentMM, $qdifferentMMsd, $qr2MM, $qr2MMsd, $qr1MM, $qr1MMsd, $qUnresolvedMM, $qUnresolvedMMsd, $qDCMM, $qDCMMsd, $qEAMM, $qEAMMsd, $qSJAMM, $qSJAMMsd, $qSSJAMM, $qSSJAMMsd, $qSDMM, $qSDMMsd;
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplMM, $ntriplMMsd, $tsolvedMM, $tsolvedMMsd, $tdifferentMM, $tdifferentMMsd,$tr2MM, $tr2MMsd, $tr1MM, $tr1MMsd, $tUnresolvedMM, $tUnresolvedMMsd, $tDCMM, $tDCMMsd, $tEAMM, $tEAMMsd, $tSJAMM, $tSJAMMsd, $tSSJAMM, $tSSJAMMsd, $tSDMM, $tSDMMsd;
								printf OUTb "* Triplets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplMM, $ntriplMMsd, $tsolvedMM, $tsolvedMMsd, $tdifferentMM, $tdifferentMMsd,$tr2MM, $tr2MMsd, $tr1MM, $tr1MMsd, $tUnresolvedMM, $tUnresolvedMMsd, $tDCMM, $tDCMMsd, $tEAMM, $tEAMMsd, $tSJAMM, $tSJAMMsd, $tSSJAMM, $tSSJAMMsd, $tSDMM, $tSDMMsd;
							}
						}
						if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
								printf OUTb "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								printf "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
								printf OUTb "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								printf "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
								printf OUTb "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
							}
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random;
								printf OUTb "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random;
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
								printf OUTb "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n",$mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
							}
						}
						print "\n";
					}
				}else {
					print "The overlap in the data between $name1_IN and $name2_IN is too separate. ($LocTrees) \n";
					print OUTb "The overlap in the data between $name1_IN and $name2_IN is too separate. ($LocTrees) \n";
					if ( ($method eq 'all') || ($method eq 'nodal') ) {
						printf "* Nodal Distance (Pruned/Unpruned): 0.000000 / 0.000000\n";
						printf OUTb "* Nodal Distance (Pruned/Unpruned): 0.000000 / 0.000000\n";
					}
					if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
						print "* Split Distance [differents/possibles]: 0 [ 0 \/ 0 ]\n";
						print OUTb "* Split Distance [differents/possibles]: 0 [ 0 \/ 0 ]\n";
					}
					if ( ($method eq 'all') || ($method eq 'disagree') ) {
						print "* Disagreement [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
						print OUTb "* Disagreement [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
					}
					if ( ($method eq 'all') || ($method eq 'quartets') ) {
						print "* Quartets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
						print OUTb "* Quartets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
					}
					if ( ($method eq 'all') || ($method eq 'triplets') ) {
						print "* Triplets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
						print OUTb "* Triplets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
					}

					if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
						if ( ($method eq 'all') || ($method eq 'nodal') ) {
							print "* Nodal Distance random (Pruned/Unpruned): ( 0.000000 +/- 0.000000 ) / ( 0.000000 +/- 0.000000 )\n";
							print OUTb "* Nodal Distance random (Pruned/Unpruned): ( 0.000000 +/- 0.000000 ) / ( 0.000000 +/- 0.000000 )\n";
						}
						if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
							print "* Split Distance random [differents/possibles]: 0 [ 0 \/ 0 ]\n";
							print OUTb "* Split Distance random [differents/possibles]: 0 [ 0 \/ 0 ]\n";
						}
						if ( ($method eq 'all') || ($method eq 'disagree') ) {
							print "* Disagreement random [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
							print OUTb "* Disagreement random [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
						}
						if ( ($method eq 'all') || ($method eq 'quartets') ) {
							print "* Quartets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
							print OUTb "* Quartets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
						}
						if ( ($method eq 'all') || ($method eq 'triplets') ) {
							print "* Triplets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
							print OUTb "* Triplets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
						}
					}
				}
				print OUTb "\n";
			}
		}
	}
	close OUTb;
	print "\nThanks for using TOPD-fMtS .................................. Bye.\n\n";
}


sub menu3() {
	($LocTrees, $randomTrees, $nombreRandomTrees, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $fileoutput,$units) = @_;

	if ($fileoutput eq "") {
		open (OUTb, ">topd_result_reference.txt");
	}else {
		open (OUTb, ">$fileoutput");
	}

	
	my %name_ref_tree = &getRefTree($LocTrees);

	
	my %name_tree = &getTrees($LocTrees); 

	my %fets = "";
	foreach $name1_IN(sort keys(%name_ref_tree)) {
		$tree1_IN = $name_tree{$name1_IN};
		foreach $name2_IN(sort keys(%name_tree)) {
			$tree2_IN = $name_tree{$name2_IN};
			if ( ($name1_IN ne $name2_IN) && ($fets{$name1_IN}{$name2_IN} eq "") ) {
				$fets{$name1_IN}{$name2_IN} = "s";
				$fets{$name2_IN}{$name1_IN} = "s";
				print OUTb "##################################### topd $name1_IN - $name2_IN #######################################\n";

				
				my ($Nsp1, $family1, $Nsp2, $family2) = &MoS($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
				
				if ( ($family1 eq "S") && ($family2 eq "S") ) {
					
					($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $tree1_IN, $name2_IN, $tree2_IN, $print_prev, $method, $level, $units);
			
				}
				elsif ( ($family1 eq "S") && ($family2 eq "M") ) {
					
					
					($resp, $ptree1_in, $ptree2_in) = &checkPruneMGF($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
			
					
					if ($resp eq "S") {
						if ($print_prev eq "n") {
						}else {
							print "\nPrunned trees are both single gene trees.\n";
						}
						($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $ptree1_in, $name2_IN, $ptree2_in, $print_prev, $method, $level, $units);
						($com, $pcCom2) = &pcCom($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
					}else {
						($topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori, $com, $pcCom2, $SplitDistSM, $SplitDistSMsd, $qSM, $qSMsd, $mSM, $mSMsd, $dif_elem, $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd, $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd, $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd) = &topdSM($name1_IN, $tree1_IN, $name2_IN, $tree2_IN, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
					}
				}elsif ( ($family1 eq "M") && ($family2 eq "S") ) {
					
					
					($resp, $ptree2_in, $ptree1_in) = &checkPruneMGF($name2_IN, $tree2_IN, $name1_IN, $tree1_IN);
					
					if ($resp eq "S") {
						if ($print_prev eq "n") {
						}else {
							print "\nPrunned trees are both single gene trees.\n";
						}
						($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $ptree1_in, $name2_IN, $ptree2_in, $print_prev, $method, $level, $units);
						($com, $pcCom2) = &pcCom($name2_IN, $tree2_IN, $name1_IN, $tree1_IN);
					}else {
						($topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori, $com, $pcCom2, $SplitDistSM, $SplitDistSMsd, $qSM, $qSMsd, $mSM, $mSMsd, $dif_elem, $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd, $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd, $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd) = &topdSM($name2_IN, $tree2_IN, $name1_IN, $tree1_IN, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
					}
				}elsif ( ($family1 eq "M") && ($family2 eq "M") ) {
					
					
					($resp, $ptree1_in, $ptree2_in) = &checkPruneMGF($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
			
					
					if ($resp eq "S") {
						if ($print_prev eq "n") {
						}else {
							print "\nPrunned trees are both single gene trees.\n";
						}
						($topd_ori, $topdUnpruned_ori, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_IN, $ptree1_in, $name2_IN, $ptree2_in, $print_prev, $method, $level, $units);
						($com, $pcCom2) = &pcCom($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
					}else {
						($topdMM_ori, $SD_ori, $topdMMunpruned_ori, $SDunpruned_ori, $com, $pcCom2, $SplitDistMM, $SplitDistMMsd, $qMM, $qMMsd, $mMM, $mMMsd, $dif_elem, $num_sp_elimMM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1MMsd, $minMM, $minMMsd, $nquartMM, $nquartMMsd, $qsolvedMM, $qsolvedMMsd, $qdifferentMM, $qdifferentMMsd, $qr2MM, $qr2MMsd, $qr1MM, $qr1MMsd, $qUnresolvedMM, $qUnresolvedMMsd, $qDCMM, $qDCMMsd, $qEAMM, $qEAMMsd, $qSJAMM, $qSJAMMsd, $qSSJAMM, $qSSJAMMsd, $qSDMM, $qSDMMsd, $ntriplMM, $ntriplMMsd, $tsolvedMM, $tsolvedMMsd, $tdifferentMM, $tdifferentMMsd,$tr2MM, $tr2MMsd, $tr1MM, $tr1MMsd, $tUnresolvedMM, $tUnresolvedMMsd, $tDCMM, $tDCMMsd, $tEAMM, $tEAMMsd, $tSJAMM, $tSJAMMsd, $tSSJAMM, $tSSJAMMsd, $tSDMM, $tSDMMsd) = &topdMM($name1_IN, $tree1_IN, $name2_IN, $tree2_IN, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
					}
				}
			
				if ($com > 3) {
					if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
						
						
						my @topd_rand =();
						my @topd_unpruned_rand = ();
						my @split_dist_rand = ();
						my @Q_rand = ();
						my @M_rand = ();

						my @A_num_sp_elim_random = ();
						my @A_e_n_taxa_tree1_random = (); 
						my @A_min_random = (); 
						my %A_DisTax_random = (); 
						my @A_nquart_random = ();
						my @A_qsolved_random = (); 
						my @A_qdifferent_random = (); 
						my @A_qr2_random = (); 
						my @A_qr1_random = (); 
						my @A_qUnresolved_random = (); 
						my @A_qDC_random = (); 
						my @A_qEA_random = (); 
						my @A_qSJA_random = (); 
						my @A_qSSJA_random = ();
						my @A_qSD_random = (); 
						my @A_ntripl_random = (); 
						my @A_tsolved_random = ();
						my @A_tdifferent_random = (); 
						my @A_tr1_random = (); 
						my @A_tr2_random = (); 
						my @A_tUnresolved_random = (); 
						my @A_tDC_random = (); 
						my @A_tEA_random = (); 
						my @A_tSJA_random = (); 
						my @A_tSSJA_random = (); 
						my @A_tSD_random = ();

						for ($ir=0;$ir<$nombreRandomTrees;$ir++) {
							if ($print_prev eq "n") {
							}else {
								print "\n### RANDOM $ir ###";
							}
							my ($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom) = ();
							if ($randomTrees eq 'guided') {
								($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom) = &random($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
							} elsif ($randomTrees eq 'random') {
								($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom) = &random2($name1_IN, $tree1_IN, $name2_IN, $tree2_IN);
							}

							
							my ($Nsp1random, $family1random, $Nsp2random, $family2random) = &MoS($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom);
							if ( ($family1random eq "S") && ($family2random eq "S") ) {
								
								my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom, $print_prev, $method, $level, $units);
								push @topd_rand, $topdrandom;
								push @topd_unpruned_rand, $topdUnprunedrandom;
								push @split_dist_rand, $SplitDistrandom;
								push @Q_rand, $qrandom;
								push @M_rand, $mrandom;

								push @A_num_sp_elim_random, $num_sp_elimrandom;
								push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
								push @A_min_random, $minrandom;
								$A_DisTax_random{$DisTaxrandom}++;
								push @A_nquart_random, $nquartrandom;
								push @A_qsolved_random, $qsolvedrandom; 
								push @A_qdifferent_random, $qdifferentrandom; 
								push @A_qr2_random, $qr2random; 
								push @A_qr1_random, $qr1random; 
								push @A_qUnresolved_random, $qUnresolvedrandom; 
								push @A_qDC_random, $qDCrandom; 
								push @A_qEA_random, $qEArandom; 
								push @A_qSJA_random, $qSJArandom; 
								push @A_qSSJA_random, $qSSJArandom;
								push @A_qSD_random, $qSDrandom; 
								push @A_ntripl_random, $ntriplrandom; 
								push @A_tsolved_random, $tsolvedrandom;
								push @A_tdifferent_random, $tdifferentrandom; 
								push @A_tr1_random, $tr1random; 
								push @A_tr2_random, $tr2random; 
								push @A_tUnresolved_random, $tUnresolvedrandom; 
								push @A_tDC_random, $tDCrandom; 
								push @A_tEA_random, $tEArandom; 
								push @A_tSJA_random, $tSJArandom; 
								push @A_tSSJA_random, $tSSJArandom; 
								push @A_tSD_random, $tSDrandom;
							}
							elsif ( ($family1random eq "S") && ($family2random eq "M") ) {
								
								
								($resp, $ptree1_in_r, $ptree2_in_r) = &checkPruneMGF($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom);
				
								
								if ($resp eq "S") {
									if ($print_prev eq "n") {
									}else {
										print "\nPrunned trees are both single gene trees.\n";
									}
									my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $ptree1_in_r, $name2_INrandom, $ptree2_in_r, $print_prev, $method, $level, $units);
									push @topd_rand, $topdrandom;
									push @topd_unpruned_rand, $topdUnprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;

									push @A_num_sp_elim_random, $num_sp_elimrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
									push @A_min_random, $minrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartrandom;
									push @A_qsolved_random, $qsolvedrandom; 
									push @A_qdifferent_random, $qdifferentrandom; 
									push @A_qr2_random, $qr2random; 
									push @A_qr1_random, $qr1random; 
									push @A_qUnresolved_random, $qUnresolvedrandom; 
									push @A_qDC_random, $qDCrandom; 
									push @A_qEA_random, $qEArandom; 
									push @A_qSJA_random, $qSJArandom; 
									push @A_qSSJA_random, $qSSJArandom;
									push @A_qSD_random, $qSDrandom; 
									push @A_ntripl_random, $ntriplrandom; 
									push @A_tsolved_random, $tsolvedrandom;
									push @A_tdifferent_random, $tdifferentrandom; 
									push @A_tr1_random, $tr1random; 
									push @A_tr2_random, $tr2random; 
									push @A_tUnresolved_random, $tUnresolvedrandom; 
									push @A_tDC_random, $tDCrandom; 
									push @A_tEA_random, $tEArandom; 
									push @A_tSJA_random, $tSJArandom; 
									push @A_tSSJA_random, $tSSJArandom; 
									push @A_tSD_random, $tSDrandom;
								}else {
									($topdSMrandom, $SDrandom, $topdSMunprunedrandom, $SDunprunedrandom, $com, $pcCom2, $SplitDistrandom, $SplitDistsdrandom, $qrandom, $qsdrandom, $mrandom, $msdrandom, $DisTaxrandom, $num_sp_elimSMrandom, $num_sp_elimSMsdrandom, $e_n_taxa_tree1SMrandom, $e_n_taxa_tree1SMsdrandom, $minSMrandom, $minSMsdrandom, $nquartSMrandom, $nquartSMsdrandom, $qsolvedSMrandom, $qsolvedSMsdrandom, $qdifferentSMrandom, $qdifferentSMsdrandom, $qr2SMrandom, $qr2SMsdrandom, $qr1SMrandom, $qr1SMsdrandom, $qUnresolvedSMrandom, $qUnresolvedSMsdrandom, $qDCSMrandom, $qDCSMsdrandom, $qEASMrandom, $qEASMsdrandom, $qSJASMrandom, $qSJASMsdrandom, $qSSJASMrandom, $qSSJASMsdrandom, $qSDSMrandom, $qSDSMsdrandom, $ntriplSMrandom, $ntriplSMsdrandom, $tsolvedSMrandom, $tsolvedSMsdrandom, $tdifferentSMrandom, $tdifferentSMsdrandom, $tr2SMrandom, $tr2SMsdrandom, $tr1SMrandom, $tr1SMsdrandom, $tUnresolvedSMrandom, $tUnresolvedSMsdrandom, $tDCSMrandom, $tDCSMsdrandom, $tEASMrandom, $tEASMsdrandom, $tSJASMrandom, $tSJASMsdrandom, $tSSJASMrandom, $tSSJASMsdrandom, $tSDSMrandom, $tSDSMsdrandom) = &topdSM($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
									push @topd_rand, $topdSMrandom;
									push @topd_unpruned_rand, $topdSMunprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;

									push @A_num_sp_elim_random, $num_sp_elimSMrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1SMrandom;
									push @A_min_random, $minSMrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartSMrandom;
									push @A_qsolved_random, $qsolvedSMrandom; 
									push @A_qdifferent_random, $qdifferentSMrandom; 
									push @A_qr2_random, $qr2SMrandom; 
									push @A_qr1_random, $qr1SMrandom; 
									push @A_qUnresolved_random, $qUnresolvedSMrandom; 
									push @A_qDC_random, $qDCSMrandom; 
									push @A_qEA_random, $qEASMrandom; 
									push @A_qSJA_random, $qSJASMrandom; 
									push @A_qSSJA_random, $qSSJASMrandom;
									push @A_qSD_random, $qSDSMrandom; 
									push @A_ntripl_random, $ntriplSMrandom; 
									push @A_tsolved_random, $tsolvedSMrandom;
									push @A_tdifferent_random, $tdifferentSMrandom; 
									push @A_tr1_random, $tr1SMrandom; 
									push @A_tr2_random, $tr2SMrandom; 
									push @A_tUnresolved_random, $tUnresolvedSMrandom; 
									push @A_tDC_random, $tDCSMrandom; 
									push @A_tEA_random, $tEASMrandom; 
									push @A_tSJA_random, $tSJASMrandom; 
									push @A_tSSJA_random, $tSSJASMrandom; 
									push @A_tSD_random, $tSDSMrandom;
								}	
							}elsif ( ($family1random eq "M") && ($family2random eq "S") ) {
								
								
								
								($resp, $ptree1_in_r, $ptree2_in_r) = &checkPruneMGF($name2_INrandom, $tree2_INrandom, $name1_INrandom, $tree1_INrandom);
								
								
								if ($resp eq "S") {
									if ($print_prev eq "n") {
									}else {
										print "\nPrunned trees are both single gene trees.\n";
									}
									my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $ptree1_in_r, $name2_INrandom, $ptree2_in_r, $print_prev, $method, $level, $units);
									push @topd_rand, $topdrandom;
									push @topd_unpruned_rand, $topdUnprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;

									push @A_num_sp_elim_random, $num_sp_elimrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
									push @A_min_random, $minrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartrandom;
									push @A_qsolved_random, $qsolvedrandom; 
									push @A_qdifferent_random, $qdifferentrandom; 
									push @A_qr2_random, $qr2random; 
									push @A_qr1_random, $qr1random; 
									push @A_qUnresolved_random, $qUnresolvedrandom; 
									push @A_qDC_random, $qDCrandom; 
									push @A_qEA_random, $qEArandom; 
									push @A_qSJA_random, $qSJArandom; 
									push @A_qSSJA_random, $qSSJArandom;
									push @A_qSD_random, $qSDrandom; 
									push @A_ntripl_random, $ntriplrandom; 
									push @A_tsolved_random, $tsolvedrandom;
									push @A_tdifferent_random, $tdifferentrandom; 
									push @A_tr1_random, $tr1random; 
									push @A_tr2_random, $tr2random; 
									push @A_tUnresolved_random, $tUnresolvedrandom; 
									push @A_tDC_random, $tDCrandom; 
									push @A_tEA_random, $tEArandom; 
									push @A_tSJA_random, $tSJArandom; 
									push @A_tSSJA_random, $tSSJArandom; 
									push @A_tSD_random, $tSDrandom;
								}else {
									($topdSMrandom, $SDrandom, $topdSMunprunedrandom, $SDunprunedrandom, $com, $pcCom2, $SplitDistrandom, $SplitDistsdrandom, $qSMrandom, $qsdrandom, $mrandom, $msdrandom, $DisTaxrandom, $num_sp_elimSMrandom, $num_sp_elimSMsdrandom, $e_n_taxa_tree1SMrandom, $e_n_taxa_tree1SMsdrandom, $minSMrandom, $minSMsdrandom, $nquartSMrandom, $nquartSMsdrandom, $qsolvedSMrandom, $qsolvedSMsdrandom, $qdifferentSMrandom, $qdifferentSMsdrandom, $qr2SMrandom, $qr2SMsdrandom, $qr1SMrandom, $qr1SMsdrandom, $qUnresolvedSMrandom, $qUnresolvedSMsdrandom, $qDCSMrandom, $qDCSMsdrandom, $qEASMrandom, $qEASMsdrandom, $qSJASMrandom, $qSJASMsdrandom, $qSSJASMrandom, $qSSJASMsdrandom, $qSDSMrandom, $qSDSMsdrandom, $ntriplSMrandom, $ntriplSMsdrandom, $tsolvedSMrandom, $tsolvedSMsdrandom, $tdifferentSMrandom, $tdifferentSMsdrandom, $tr2SMrandom, $tr2SMsdrandom, $tr1SMrandom, $tr1SMsdrandom, $tUnresolvedSMrandom, $tUnresolvedSMsdrandom, $tDCSMrandom, $tDCSMsdrandom, $tEASMrandom, $tEASMsdrandom, $tSJASMrandom, $tSJASMsdrandom, $tSSJASMrandom, $tSSJASMsdrandom, $tSDSMrandom, $tSDSMsdrandom) = &topdSM($name2_INrandom, $tree2_INrandom, $name1_INrandom, $tree1_INrandom, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
									push @topd_rand, $topdSMrandom;
									push @topd_unpruned_rand, $topdSMunprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;

									push @A_num_sp_elim_random, $num_sp_elimSMrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1SMrandom;
									push @A_min_random, $minSMrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartSMrandom;
									push @A_qsolved_random, $qsolvedSMrandom; 
									push @A_qdifferent_random, $qdifferentSMrandom; 
									push @A_qr2_random, $qr2SMrandom; 
									push @A_qr1_random, $qr1SMrandom; 
									push @A_qUnresolved_random, $qUnresolvedSMrandom; 
									push @A_qDC_random, $qDCSMrandom; 
									push @A_qEA_random, $qEASMrandom; 
									push @A_qSJA_random, $qSJASMrandom; 
									push @A_qSSJA_random, $qSSJASMrandom;
									push @A_qSD_random, $qSDSMrandom; 
									push @A_ntripl_random, $ntriplSMrandom; 
									push @A_tsolved_random, $tsolvedSMrandom;
									push @A_tdifferent_random, $tdifferentSMrandom; 
									push @A_tr1_random, $tr1SMrandom; 
									push @A_tr2_random, $tr2SMrandom; 
									push @A_tUnresolved_random, $tUnresolvedSMrandom; 
									push @A_tDC_random, $tDCSMrandom; 
									push @A_tEA_random, $tEASMrandom; 
									push @A_tSJA_random, $tSJASMrandom; 
									push @A_tSSJA_random, $tSSJASMrandom; 
									push @A_tSD_random, $tSDSMrandom;
								}	 
							}elsif ( ($family1random eq "M") && ($family2random eq "M") ) {
								
								
								($resp, $ptree1_in_r, $ptree2_in_r) = &checkPruneMGF($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom);
								
								if ($resp eq "S") {
									if ($print_prev eq "n") {
									}else {
										print "\nPrunned trees are both single gene trees.\n";
									}
									my ($topdrandom, $topdUnprunedrandom, $com, $pcCom2, $SplitDistrandom, $qrandom, $mrandom, $num_sp_elimrandom, $e_n_taxa_tree1random, $minrandom, $DisTaxrandom, $nquartrandom, $qsolvedrandom, $qdifferentrandom, $qr2random, $qr1random, $qUnresolvedrandom, $qDCrandom, $qEArandom, $qSJArandom, $qSSJArandom, $qSDrandom, $ntriplrandom, $tsolvedrandom, $tdifferentrandom, $tr1random, $tr2random, $tUnresolvedrandom, $tDCrandom, $tEArandom, $tSJArandom, $tSSJArandom, $tSDrandom) = &topd($name1_INrandom, $ptree1_in_r, $name2_INrandom, $ptree2_in_r, $print_prev, $method, $level, $units);
									push @topd_rand, $topdrandom;
									push @topd_unpruned_rand, $topdUnprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;

									push @A_num_sp_elim_random, $num_sp_elimrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1random;
									push @A_min_random, $minrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartrandom;
									push @A_qsolved_random, $qsolvedrandom; 
									push @A_qdifferent_random, $qdifferentrandom; 
									push @A_qr2_random, $qr2random; 
									push @A_qr1_random, $qr1random; 
									push @A_qUnresolved_random, $qUnresolvedrandom; 
									push @A_qDC_random, $qDCrandom; 
									push @A_qEA_random, $qEArandom; 
									push @A_qSJA_random, $qSJArandom; 
									push @A_qSSJA_random, $qSSJArandom;
									push @A_qSD_random, $qSDrandom; 
									push @A_ntripl_random, $ntriplrandom; 
									push @A_tsolved_random, $tsolvedrandom;
									push @A_tdifferent_random, $tdifferentrandom; 
									push @A_tr1_random, $tr1random; 
									push @A_tr2_random, $tr2random; 
									push @A_tUnresolved_random, $tUnresolvedrandom; 
									push @A_tDC_random, $tDCrandom; 
									push @A_tEA_random, $tEArandom; 
									push @A_tSJA_random, $tSJArandom; 
									push @A_tSSJA_random, $tSSJArandom; 
									push @A_tSD_random, $tSDrandom;
								}else {
									($topdMMrandom, $SDrandom, $topdMMunprunedrandom, $SDunprunedrandom, $com, $pcCom2, $SplitDistrandom, $SplitDistsdrandom, $qrandom, $qsdrandom, $mrandom, $msdrandom, $DisTaxrandom, $num_sp_elimMMrandom, $num_sp_elimMMsdrandom, $e_n_taxa_tree1MMrandom, $e_n_taxa_tree1MMsdrandom, $minMMrandom, $minMMsdrandom, $nquartMMrandom, $nquartMMsdrandom, $qsolvedMMrandom, $qsolvedMMsdrandom, $qdifferentMMrandom, $qdifferentMMsdrandom, $qr2MMrandom, $qr2MMsdrandom, $qr1MMrandom, $qr1MMsdrandom, $qUnresolvedMMrandom, $qUnresolvedMMsdrandom, $qDCMMrandom, $qDCMMsdrandom, $qEAMMrandom, $qEAMMsdrandom, $qSJAMMrandom, $qSJAMMsdrandom, $qSSJAMMrandom, $qSSJAMMsdrandom, $qSDMMrandom, $qSDMMsdrandom, $ntriplMMrandom, $ntriplMMsdrandom, $tsolvedMMrandom, $tsolvedMMsdrandom, $tdifferentMMrandom, $tdifferentMMsdrandom, $tr2MMrandom, $tr2MMsdrandom, $tr1MMrandom, $tr1MMsdrandom, $tUnresolvedMMrandom, $tUnresolvedMMsdrandom, $tDCMMrandom, $tDCMMsdrandom, $tEAMMrandom, $tEAMMsdrandom, $tSJAMMrandom, $tSJAMMsdrandom, $tSSJAMMrandom, $tSSJAMMsdrandom, $tSDMMrandom, $tSDMMsdrandom) = &topdMM($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level, $units);
									push @topd_rand, $topdMMrandom;
 									push @topd_unpruned_rand, $topdMMunprunedrandom;
									push @split_dist_rand, $SplitDistrandom;
									push @Q_rand, $qrandom;
									push @M_rand, $mrandom;
				
									push @A_num_sp_elim_random, $num_sp_elimMMrandom;
									push @A_e_n_taxa_tree1_random, $e_n_taxa_tree1MMrandom;
									push @A_min_random, $minMMrandom;
									$A_DisTax_random{$DisTaxrandom}++;
									push @A_nquart_random, $nquartMMrandom;
									push @A_qsolved_random, $qsolvedMMrandom; 
									push @A_qdifferent_random, $qdifferentMMrandom; 
									push @A_qr2_random, $qr2MMrandom; 
									push @A_qr1_random, $qr1MMrandom; 
									push @A_qUnresolved_random, $qUnresolvedMMrandom; 
									push @A_qDC_random, $qDCMMrandom; 
									push @A_qEA_random, $qEAMMrandom; 
									push @A_qSJA_random, $qSJAMMrandom; 
									push @A_qSSJA_random, $qSSJAMMrandom;
									push @A_qSD_random, $qSDMMrandom; 
									push @A_ntripl_random, $ntriplMMrandom; 
									push @A_tsolved_random, $tsolvedMMrandom;
									push @A_tdifferent_random, $tdifferentMMrandom; 
									push @A_tr1_random, $tr1MMrandom; 
									push @A_tr2_random, $tr2MMrandom; 
									push @A_tUnresolved_random, $tUnresolvedMMrandom; 
									push @A_tDC_random, $tDCMMrandom; 
									push @A_tEA_random, $tEAMMrandom; 
									push @A_tSJA_random, $tSJAMMrandom; 
									push @A_tSSJA_random, $tSSJAMMrandom; 
									push @A_tSD_random, $tSDMMrandom;
								}
							}
						}
						
						($mean_rand,$sd_rand)=&meanSD(@topd_rand);
						($mean_unpruned_rand,$sd_unpruned_rand)=&meanSD(@topd_unpruned_rand);
						($mean_rand_split,$sd_rand_split)=&meanSD(@split_dist_rand);
						($mean_rand_q,$sd_rand_q)=&meanSD(@Q_rand);
						($mean_rand_m,$sd_rand_m)=&meanSD(@M_rand);

						($mean_num_sp_elim_random, $sd_num_sp_elim_random) = &meanSD(@A_num_sp_elim_random);
						($mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random) = &meanSD(@A_e_n_taxa_tree1_random);
						($mean_min_random, $sd_min_random) = &meanSD(@A_min_random);
						$dif_elemrandom = "";
						foreach $DisTaxSMelementrandom(sort keys(%A_DisTax_random)) {
							$aux_dist__random = $A_DisTax_random{$DisTaxSMelementrandom}; 
							$DisTaxSMelementrandom =~ s/-/ /g;
							$dif_elemrandom .= "[$DisTaxSMelementrandom("."$aux_dist__random".")] ";
						}
						($mean_nquart_random, $sd_nquart_random) = &meanSD(@A_nquart_random);
						($mean_qsolved_random, $sd_qsolved_random) = &meanSD(@A_qsolved_random);
						($mean_qdifferent_random, $sd_qdifferent_random) = &meanSD(@A_qdifferent_random);
						($mean_qr2_random, $sd_qr2_random) = &meanSD(@A_qr2_random);
						($mean_qr1_random, $sd_qr1_random) = &meanSD(@A_qr1_random);
						($mean_qUnresolved_random, $sd_qUnresolved_random) = &meanSD(@A_qUnresolved_random);
						($mean_qDC_random, $sd_qDC_random) = &meanSD(@A_qDC_random);
						($mean_qEA_random, $sd_qEA_random) = &meanSD(@A_qEA_random);
						($mean_qSJA_random, $sd_qSJA_random) = &meanSD(@A_qSJA_random);
						($mean_qSSJA_random, $sd_qSSJA_random) = &meanSD(@A_qSSJA_random);
						($mean_qSD_random, $sd_qSD_random) = &meanSD(@A_qSD_random);
						($mean_ntripl_random, $sd_ntripl_random) = &meanSD(@A_ntripl_random);
						($mean_tsolved_random, $sd_tsolved_random) = &meanSD(@A_tsolved_random);
						($mean_tdifferent_random, $sd_tdifferent_random) = &meanSD(@A_tdifferent_random);
						($mean_tr1_random, $sd_tr1_random) = &meanSD(@A_tr1_random);
						($mean_tr2_random, $sd_tr2_random) = &meanSD(@A_tr2_random);
						($mean_tUnresolved_random, $sd_tUnresolved_random) = &meanSD(@A_tUnresolved_random);
						($mean_tDC_random, $sd_tDC_random) = &meanSD(@A_tDC_random);
						($mean_tEA_random, $sd_tEA_random) = &meanSD(@A_tEA_random);
						($mean_tSJA_random, $sd_tSJA_random) = &meanSD(@A_tSJA_random);
						($mean_tSSJA_random, $sd_tSSJA_random) = &meanSD(@A_tSSJA_random);
						($mean_tSD_random, $sd_tSD_random) = &meanSD(@A_tSD_random);
					}
			
					
					print "\n############################################ RESULTS $name1_IN - $name2_IN ###################################\n\n";
					if ($resp eq "S")  {
						printf "* Percentage of taxa in common: %6.1f", $pcCom2;
						print "%\n";
						printf OUTb "* Percentage of taxa in common: %6.1f", $pcCom2;
						print OUTb "%\n";
						$topdUnpruned_ori = $topd_ori * (1+(1-($pcCom2 /100)));

						if ( ($method eq 'all') || ($method eq 'nodal') ) {
							printf "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
							printf OUTb "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
						}
						if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
							print "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
							print OUTb "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
						}
						if ( ($method eq 'all') || ($method eq 'disagree') ) {
							$DisTax =~ s/-/ /g;
							print OUTb "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
							print "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
						}
						if ( ($method eq 'all') || ($method eq 'quartets') ) {
							printf "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
							printf OUTb "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
						}
						if ( ($method eq 'all') || ($method eq 'triplets') ) {
							printf "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
							printf OUTb "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
						}

						if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
							$mean_unpruned_rand = $mean_rand * (1+(1-($pcCom2 /100)));
							$sd_unpruned_rand = $sd_rand * (1+(1-($pcCom2 /100)));
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
								printf OUTb "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								printf "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
								printf OUTb "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								printf "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
								printf OUTb "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
							}
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random;
								printf OUTb "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random;
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
								printf OUTb "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n",$mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
							}
						}
						print "\n";
					}
					else {
						printf "* Percentage of taxa in common: %6.1f", $pcCom2;
						print "%\n";
						printf OUTb "* Percentage of taxa in common: %6.1f", $pcCom2;
						print OUTb "%\n";
						if ( ($family1 eq "S") && ($family2 eq "S") ) {
							
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
								printf OUTb "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f\n", $topd_ori, $topdUnpruned_ori;
							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								print "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
								print OUTb "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								$DisTax =~ s/-/ /g;
								print OUTb "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
								print "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( $DisTax)\n";	
							}
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
								printf OUTb "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
								printf OUTb "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
							}
						}
						elsif ( ($family1 eq "S") && ($family2 eq "M") ) {
							
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "\n* Nodal Distance SM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;
								printf OUTb "* Nodal Distance SM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;
							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								printf "* Split Distance SM [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
								printf OUTb "* Split Distance SM [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								printf "* Disagreement SM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
								printf OUTb "* Disagreement SM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
							}
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
								printf OUTb "* Quartets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
								printf OUTb "* Triplets SM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
							}
						}elsif ( ($family1 eq "M") && ($family2 eq "S") ) {
							
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "\n* Nodal Distance MS (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;
								printf OUTb "* Nodal Distance MS (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdSM_ori, $SD_ori, $topdSMunpruned_ori, $SDunpruned_ori;							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								printf "* Split Distance MS [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
								printf OUTb "* Split Distance MS [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								printf "* Disagreement MS [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
								printf OUTb "* Disagreement MS [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
							}
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
								printf OUTb "* Quartets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
								printf OUTb "* Triplets MS: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
							}
						}elsif ( ($family1 eq "M") && ($family2 eq "M") ) {
							
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "* Nodal Distance MM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdMM_ori, $SD_ori, $topdMMunpruned_ori, $SDunpruned_ori;
								printf OUTb "* Nodal Distance MM (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n",$topdMM_ori, $SD_ori, $topdMMunpruned_ori, $SDunpruned_ori;
							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								printf "* Split Distance MM [differents/possibles]: $SplitDistMM +/- %5.3f [ $qMM +/- %5.3f \/ $mMM +/- %5.3f ]\n", $SplitDistMMsd, $qMMsd, $mMMsd ;
								printf OUTb "* Split Distance MM [differents/possibles]: $SplitDistMM +/- %5.3f [ $qMM +/- %5.3f \/ $mMM +/- %5.3f ]\n", $SplitDistMMsd, $qMMsd, $mMMsd ;
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								printf "* Disagreement MM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimMM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1MMsd, $minMM, $minMMsd ;
								printf OUTb "* Disagreement MM [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimMM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1MMsd, $minMM, $minMMsd ;
							}
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartMM, $nquartMMsd, $qsolvedMM, $qsolvedMMsd, $qdifferentMM, $qdifferentMMsd, $qr2MM, $qr2MMsd, $qr1MM, $qr1MMsd, $qUnresolvedMM, $qUnresolvedMMsd, $qDCMM, $qDCMMsd, $qEAMM, $qEAMMsd, $qSJAMM, $qSJAMMsd, $qSSJAMM, $qSSJAMMsd, $qSDMM, $qSDMMsd;
								printf OUTb "* Quartets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartMM, $nquartMMsd, $qsolvedMM, $qsolvedMMsd, $qdifferentMM, $qdifferentMMsd, $qr2MM, $qr2MMsd, $qr1MM, $qr1MMsd, $qUnresolvedMM, $qUnresolvedMMsd, $qDCMM, $qDCMMsd, $qEAMM, $qEAMMsd, $qSJAMM, $qSJAMMsd, $qSSJAMM, $qSSJAMMsd, $qSDMM, $qSDMMsd;
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplMM, $ntriplMMsd, $tsolvedMM, $tsolvedMMsd, $tdifferentMM, $tdifferentMMsd,$tr2MM, $tr2MMsd, $tr1MM, $tr1MMsd, $tUnresolvedMM, $tUnresolvedMMsd, $tDCMM, $tDCMMsd, $tEAMM, $tEAMMsd, $tSJAMM, $tSJAMMsd, $tSSJAMM, $tSSJAMMsd, $tSDMM, $tSDMMsd;
								printf OUTb "* Triplets MM: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplMM, $ntriplMMsd, $tsolvedMM, $tsolvedMMsd, $tdifferentMM, $tdifferentMMsd,$tr2MM, $tr2MMsd, $tr1MM, $tr1MMsd, $tUnresolvedMM, $tUnresolvedMMsd, $tDCMM, $tDCMMsd, $tEAMM, $tEAMMsd, $tSJAMM, $tSJAMMsd, $tSSJAMM, $tSSJAMMsd, $tSDMM, $tSDMMsd;
							}
						}
						if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
							if ( ($method eq 'all') || ($method eq 'nodal') ) {
								printf "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
								printf OUTb "* Nodal Distance random (Pruned/Unpruned): ( %6.6f +/- %6.6f ) / ( %6.6f +/- %6.6f )\n", $mean_rand, $sd_rand,$mean_unpruned_rand,$sd_unpruned_rand;
							}
							if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
								printf "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
								printf OUTb "* Split Distance random [differents/possibles]: $mean_rand_split +/- %5.3f [ $mean_rand_q +/- %5.3f \/ $mean_rand_m +/- %5.3f ]\n", $sd_rand_split, $sd_rand_q, $sd_rand_m ;
							}
							if ( ($method eq 'all') || ($method eq 'disagree') ) {
								printf "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
								printf OUTb "* Disagreement random [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elemrandom\n", $mean_num_sp_elim_random, $sd_num_sp_elim_random, $mean_e_n_taxa_tree1_random, $sd_e_n_taxa_tree1_random, $mean_min_random, $sd_min_random ;
							}
							if ( ($method eq 'all') || ($method eq 'quartets') ) {
								printf "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random;
								printf OUTb "* Quartets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_nquart_random, $sd_nquart_random, $mean_qsolved_random, $sd_qsolved_random, $mean_qdifferent_random, $sd_qdifferent_random, $mean_qr2_random, $sd_qr2_random, $mean_qr1_random, $sd_qr1_random, $mean_qUnresolved_random, $sd_qUnresolved_random, $mean_qDC_random, $sd_qDC_random, $mean_qEA_random, $sd_qEA_random, $mean_qSJA_random, $sd_qSJA_random, $mean_qSSJA_random, $sd_qSSJA_random, $mean_qSD_random, $sd_qSD_random;
							}
							if ( ($method eq 'all') || ($method eq 'triplets') ) {
								printf "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
								printf OUTb "* Triplets random: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n",$mean_ntripl_random, $sd_ntripl_random, $mean_tsolved_random, $sd_tsolved_random, $mean_tdifferent_random, $sd_tdifferent_random, $mean_tr1_random, $sd_tr1_random, $mean_tr2_random, $sd_tr2_random, $mean_tUnresolved_random, $sd_tUnresolved_random, $mean_tDC_random, $sd_tDC_random, $mean_tEA_random, $sd_tEA_random, $mean_tSJA_random, $sd_tSJA_random, $mean_tSSJA_random, $sd_tSSJA_random, $mean_tSD_random, $sd_tSD_random;
							}
						}
						print "\n";
					}
				}else {
					print "The overlap in the data between $name1_IN and $name2_IN is too separate. ($LocTrees) \n";
					print OUTb "The overlap in the data between $name1_IN and $name2_IN is too separate. ($LocTrees) \n";
					if ( ($method eq 'all') || ($method eq 'nodal') ) {
						printf "* Nodal Distance (Pruned/Unpruned): 0.000000 / 0.000000\n";
						printf OUTb "* Nodal Distance (Pruned/Unpruned): 0.000000 / 0.000000\n";
					}
					if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
						print "* Split Distance [differents/possibles]: $SplitDist [ $q \/ $m ]\n";
						print OUTb "* Split Distance [differents/possibles]: 0 [ 0 \/ 0 ]\n";
					}
					if ( ($method eq 'all') || ($method eq 'disagree') ) {
						print "* Disagreement [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
						print OUTb "* Disagreement [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
					}
					if ( ($method eq 'all') || ($method eq 'quartets') ) {
						print "* Quartets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
						print OUTb "* Quartets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
					}
					if ( ($method eq 'all') || ($method eq 'triplets') ) {
						print "* Triplets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
						print OUTb "* Triplets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
					}

					if ( ($randomTrees eq 'guided') || ($randomTrees eq 'random') ) {
						if ( ($method eq 'all') || ($method eq 'nodal') ) {
							print "* Nodal Distance random (Pruned/Unpruned): ( 0.000000 +/- 0.000000 ) / ( 0.000000 +/- 0.000000 )\n";
							print OUTb "* Nodal Distance random (Pruned/Unpruned): ( 0.000000 +/- 0.000000 ) / ( 0.000000 +/- 0.000000 )\n";
						}
						if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
							print "* Split Distance random [differents/possibles]: 0 [ 0 \/ 0 ]\n";
							print OUTb "* Split Distance random [differents/possibles]: 0 [ 0 \/ 0 ]\n";
						}
						if ( ($method eq 'all') || ($method eq 'disagree') ) {
							print "* Disagreement random [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
							print OUTb "* Disagreement random [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: (none)\n";
						}
						if ( ($method eq 'all') || ($method eq 'quartets') ) {
							print "* Quartets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
							print OUTb "* Quartets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
						}
						if ( ($method eq 'all') || ($method eq 'triplets') ) {
							print "* Triplets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
							print OUTb "* Triplets random: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
						}
					}
				}
				print OUTb "\n";
			}
		}
	}
	close OUTb;
	print "\nThanks for using TOPD-fMtS .................................. Bye.\n\n";
}


sub getTrees() {
	my $LocTrees = shift;
	my $topd = 1;
	my @TREES = ();

	
	open (IN, "<$LocTrees") || die "\n#ERROR# Cannot open $LocTrees\n";
	print "\n################################################ INPUT ##############################################\n";
	$numtrees = 0;
	while (<IN>) {
		$numtrees++;
		chomp;
		($name, $tree) = split '\t', $_;
		if ( ($name eq '') || ( $tree eq '') ) {
			$tree = $name;
			$name = "Tree"."$numtrees";
		}
		push @TREES, $name;
		$tree =~ s/e-\d+//g;
		$tree =~ s/:-\d+\.\d+//g;
		$tree =~ s/:\d+\.\d+//g;
		$tree =~ s/-\d+\.\d+//g;
		$tree =~ s/\d+\.\d+//g;
		$tree =~ s/:\d+//g;
		$tree =~ s/://g;
		$tree =~ s/-//g;
		$tree =~ s/_//g;
		$tree =~ s/\s+//;
		$tree =~ s/\t+//;
		$tree =~ s/;.+/;/;
		$tree =~ s/\)\d+/\)/g;
		push @TREES, $tree;
		print "\n* $name: $tree\n";
	}
	close IN;
	return (@TREES);
}


sub getRefTree() {
	my $LocTrees = shift;
	my $topd = 1;
	my @TREES = ();

	
	open (IN, "<$LocTrees") || die "\n#ERROR# Cannot open $LocTrees\n";
	print "\n################################################ REFERENCE TREE ######################################\n";
	$numtrees = 0;
	while (<IN>) {
		$numtrees++;
		chomp;
		($name, $tree) = split '\t', $_;
		if ( ($name eq '') || ( $tree eq '') ) {
			$tree = $name;
			$name = "Tree"."$numtrees";
		}
		push @TREES, $name;
		$tree =~ s/e-\d+//g;
		$tree =~ s/:-\d+\.\d+//g;
		$tree =~ s/:\d+\.\d+//g;
		$tree =~ s/-\d+\.\d+//g;
		$tree =~ s/\d+\.\d+//g;
		$tree =~ s/:\d+//g;
		$tree =~ s/://g;
		$tree =~ s/-//g;
		$tree =~ s/_//g;
		$tree =~ s/\s+//;
		$tree =~ s/\t+//;
		$tree =~ s/;.+/;/;
		$tree =~ s/\)\d+/\)/g;
		push @TREES, $tree;
		print "\n* $name: $tree\n";
		last;
	}
	close IN;
	return (@TREES);
}


sub checkPruneMGF {
	my ($name1_IN, $tree1_IN, $name2_IN, $tree2_IN) = @_;
	my $resp = "M";
	my %trees = ();	
	$trees{$name1_IN} = $tree1_IN;
	$trees{$name2_IN} = $tree2_IN;
	
	
	my $sp_tree1 = $tree1_IN;
	$sp_tree1 =~ s/\W+/\t/g;
	my @sptree1 = split '\t', $sp_tree1;

	my $sp_tree2 = $tree2_IN;
	$sp_tree2 =~ s/\W+/\t/g;
	my @sptree2 = split '\t', $sp_tree2;
	my @spNoCom=();
	my $comCheck=0;
	
	foreach $tr1(@sptree1) {
		$trobat = "false";
		foreach $tr2(@sptree2) {
			if ($tr1 eq $tr2) {
				$trobat = "true";
			}
		}
		if ($trobat eq "false") {
			push @spNoCom, $tr1;
		}
	}

	foreach $tr2(@sptree2) {
		$trobat = "false";
		foreach $tr1(@sptree1) {
			if ($tr1 eq $tr2) {
				$trobat = "true";
				$comCheck++
			}
		}
		if ($trobat eq "false") {
			push @spNoCom, $tr2;
		}
	}

	$comCheck /= 2;
	my ($Nsp1, $family1, $Nsp2, $family2) = ();
	if ($comCheck >=4) {
		
		foreach $tree(sort keys(%trees)) {
			$pTree = &pruneTree($trees{$tree}, @spNoCom);
			
			$trees{$tree} = $pTree;
		}
		
		
		foreach $tree(sort keys(%trees)) {
			$pTree = &unRootTrees($trees{$tree});
			$trees{$tree} = $pTree;
			
		}
	
		
		($Nsp1, $family1, $Nsp2, $family2) = &MoS($name1_IN, $trees{$name1_IN}, $name2_IN, $trees{$name2_IN});
		if ( ($family1 eq "S") && ($family2 eq "S") ) {
			$resp = "S";
		}
	}else {
		$resp = "S";
	}

	return ($resp, $trees{$name1_IN}, $trees{$name2_IN});
}


sub pcCom() {
	my ($name1_IN, $tree1_IN, $name2_IN, $tree2_IN) = @_;
	my %trees = ();	
	$trees{$name1_IN} = $tree1_IN;
	$trees{$name2_IN} = $tree2_IN;

	
	my %sp1 = ();
	my %sp2=();
	my %spTotal=();
	my $i=1;
	my $j=0;

	foreach $tree(sort keys(%trees)) {
		$treeSP = $trees{$tree};
		
		while ($treeSP =~ /(\w+)/) {
			$sp = $1;
			$treeSP =~ s/$sp//;
			if ($i == 1) {
				$sp1{$sp} = "$i";	
				$j++;
				
			}elsif ($i == 2) {
				$sp2{$sp} = "$i";
				$j++;
				
			}
			$spTotal{$sp} = $i;
		}
		
		$j=0;
		$i=2;
	}

	
	my %spCom=();
	my %spNoCom=();
	my @spNoCom=();
	my $com = "";
	my $nocom = "";
	foreach $sp_(keys%spTotal) {
		if ( ($sp1{$sp_} == 1) && ($sp2{$sp_} == 2) ){
			$spCom{$sp_} = "12";
			$com++;
		}elsif ( ($sp1{$sp_} == 1) ){
			push @spNoCom, $sp_;
			$nocom++;
		}elsif ( ($sp2{$sp_} == 2) ){
			push @spNoCom, $sp_;
			$nocom++;
		}
		else {
		}
	}
	
	
	

	
	

	
	my $pcCom2 = ( ($com) / (($com) + $nocom) ) * 100;

	return ($com, $pcCom2);
}


sub topd() {
	my ($name1_IN, $tree1_IN, $name2_IN, $tree2_IN, $print_prev, $method, $level,$units) = @_;
	my %trees = ();	
	$trees{$name1_IN} = $tree1_IN;
	$trees{$name2_IN} = $tree2_IN;

	my $topd = 1;
	my ($SplitDist, $p, $q, $m) = ();

	
	my %sp1 = ();
	my %sp2=();
	my %spTotal=();
	my $i=1;
	my $j=0;
	if ($print_prev eq 'n') {
	}else{
		print "\n###################################### topd $name1_IN - $name2_IN #########################################\n";
	}
	foreach $tree(sort keys(%trees)) {
		$treeSP = $trees{$tree};
		
		while ($treeSP =~ /(\w+)/) {
			$sp = $1;
			$treeSP =~ s/$sp//;
			if ($i == 1) {
				$sp1{$sp} = "$i";	
				$j++;
				
			}elsif ($i == 2) {
				$sp2{$sp} = "$i";
				$j++;
				
			}
			$spTotal{$sp} = $i;
		}
		
		$j=0;
		$i=2;
	}

	
	my %spCom=();
	my %spNoCom=();
	my @spNoCom=();
	my $com = "";
	my $nocom = "";
	if ($print_prev eq 'n') {
	}else{
		print "\n* TAXA IN COMMON: ";
	}
	foreach $sp_(keys%spTotal) {
		if ( ($sp1{$sp_} == 1) && ($sp2{$sp_} == 2) ){
			if ($print_prev eq "n") {
			}else {
				print "$sp_ ";
			}
			$spCom{$sp_} = "12";
			$com++;
		}elsif ( ($sp1{$sp_} == 1) ){
			push @spNoCom, $sp_;
			$nocom++;
		}elsif ( ($sp2{$sp_} == 2) ){
			push @spNoCom, $sp_;
			$nocom++;
		}
		else {
		}
	}

	
	
	

	
	

	
	my $pcCom2 = ( ($com) / (($com) + $nocom) ) * 100;

	if ($print_prev eq 'n') {
	}else{
		printf "\n* Percentage of taxa in common: %6.1f \% \n\n", $pcCom2;
	}

	
	
	
	my ($min, $num_sp_elim, $e_n_taxa_tree1) = (0,0,0);
	my $DisTax = "";
	my @elim_e = ();
	my ($SplitDist, $p, $m) = (0,0,0);
	my ($nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD) = (0,0,0,0,0,0,0,0,0,0,0,);
	my ($ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD)  = (0,0,0,0,0,0,0,0,0,0,0,);
	my @splittrees = ();

	if ($com <= 3) {
		($topd, $topdUnpruned, $q, $m, $SplitDist, $num_sp_elim, $e_n_taxa_tree1, $min, $nquart, $qsolved, $qdifferent, $qr1, $qr2, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = (0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
		@elim_e = ();

		if ($print_prev eq 'n') {
		}else{
			print "\n";
			printf "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f \n\n", $topd, $topdUnpruned;
			printf "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f \n\n", $topd, $topdUnpruned;
			printf "* Split Distance [differents/possibles]: %6.6f [ $q \/ $m ]\n\n", $SplitDist;
			print "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( @elim_e )\n";
			printf "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
			printf "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
		}
	} else {
		
		foreach $tree(sort keys(%trees)) {
			$pTree = &pruneTree($trees{$tree}, @spNoCom);
			
			$trees{$tree} = $pTree;
		}
		
		
		@splittrees = ();
		foreach $tree(sort keys(%trees)) {
			$pTree = &unRootTrees($trees{$tree});
			$trees{$tree} = $pTree;
			push @splittrees, $pTree;
		if ($print_prev eq 'n') {
			}else{
				print "* PRUNED $tree: $trees{$tree}\n";
				}
			}
	
			
			%sp_fetes_dist = ();
			$contCOMP=0;
			$DIST12 = 0;

			if ($print_prev eq 'n') {
			}else{
				print "\n";
			}
			
			
			
			if ( ($method eq 'all') || ($method eq 'nodal') ) {	

				
				foreach $s1(sort keys(%spCom)) {
					foreach $s2(sort keys(%spCom)) {
						if ($sp_fetes_dist{$s2}{$s1} ne 'y') {
		
							$sp_fetes_dist{$s1}{$s2} = 'y';
		
							
							$tree_m1 = $trees{$name1_IN};
							if ($s1 eq $s2) {
								
								$dist1{$s1}{$s2} = 0;
								
							}
							elsif ($tree_m1 =~ /($s1),($s2)/) {
								
								$dist1{$s1}{$s2} = 1;
								
							}
							elsif ($tree_m1 =~ /($s2),($s1)/) {
								
								$dist1{$s1}{$s2} = 1;
								
							}
							else {
								
								
								if ($tree_m1 =~ /($s1)(\W.*\W)($s2)/) {
									$ds = $2;
									
									while ( $ds =~ /(\([\w*,]+\))/)  {
										$ds =~ s/\($1\)//;
										
									}
									$pos = $ds =~ s/\)/\)/g;
									$pos += $ds =~ s/\(/\(/g;
			
									$valor = $pos+1;
				
									$dist1{$s1}{$s2} = $valor;
									
								}
								elsif ($tree_m1 =~ /($s2)(\W.*\W)($s1)/) {
									$ds = $2;
									
									while ( $ds =~ /(\([\w*,]+\))/)  {
										$ds =~ s/\($1\)//;
										
									}
									$pos = $ds =~ s/\)/\)/g;
									$pos += $ds =~ s/\(/\(/g;
			
									$valor = $pos+1;
				
									$dist1{$s1}{$s2} = $valor;
									
								}
							}
			
							
		
							$tree_m2 = $trees{$name2_IN};
		
							if ($s1 eq $s2) {
								
								$dist2{$s1}{$s2} = 0;
								
							}
							elsif ($tree_m2 =~ /($s1),($s2)/) {
								
								$dist2{$s1}{$s2} = 1;
								
							}
							elsif ($tree_m2 =~ /($s2),($s1)/) {
								
								$dist2{$s1}{$s2} = 1;
								
							}
							else {
								
								
								if ($tree_m2 =~ /($s1)(\W.*\W)($s2)/) {
									$ds = $2;
									
									while ( $ds =~ /(\([\w*,]+\))/)  {
										$ds =~ s/\($1\)//;
										
									}
									$pos = $ds =~ s/\)/\)/g;
									$pos += $ds =~ s/\(/\(/g;
			
									$valor = $pos+1;
				
									$dist2{$s1}{$s2} = $valor;
									
								}
								elsif ($tree_m2 =~ /($s2)(\W.*\W)($s1)/) {
									$ds = $2;
									
									while ( $ds =~ /(\([\w*,]+\))/)  {
										$ds =~ s/\($1\)//;
										
									}
									$pos = $ds =~ s/\)/\)/g;
									$pos += $ds =~ s/\(/\(/g;
			
									$valor = $pos+1;
				
									$dist2{$s1}{$s2} = $valor;
									
								}
							}
						}else {
							
							$dist12 = ($dist1{$s2}{$s1}-$dist2{$s2}{$s1})**2;
							$DIST12 += $dist12;
							$contCOMP++;
							
						}
					}
				}
				
				$topd = sqrt($DIST12 / $contCOMP);
				$topdUnpruned = $topd *	(1+(1-($pcCom2 /100)));
				
				
				if ($print_prev eq 'n') {
				}else{
					printf "* Nodal Distance (Pruned/Unpruned): %6.6f / %6.6f \n\n", $topd, $topdUnpruned;
				}
			}

		
		
		
 		if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
			($SplitDist, $p, $m) = &splitDistance(@splittrees);
			$q = $m-$p;

			
			if ($print_prev eq 'n') {
			}else{
				printf "* Split Distance [differents/possibles]: %6.6f [ $q \/ $m ]\n\n", $SplitDist;
			}
		}

		
		
		
		if ( ($method eq 'all') || ($method eq 'disagree') ) {
			($min, $num_sp_elim, $e_n_taxa_tree1, @elim_e) = &elimTaxa($SplitDist, $level, @splittrees, $print_prev);

			
			$DisTax = "";
			foreach $e_lim_(@elim_e) {
				$DisTax .= "$e_lim_-";
			}
			
			if ($print_prev eq 'n') {
			}else{
				print "* Disagreement [ taxa disagree / all taxa ]: [ $num_sp_elim / $e_n_taxa_tree1 ], New Split Distance: $min, Taxa disagree: ( @elim_e )\n\n";
			}
		}

		
		
		
		if ( ($method eq 'all') || ($method eq 'quartets') ) {
			($nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD) = &quartetsDistance(@splittrees, $units);

			
			if ($print_prev eq 'n') {
			}else{
				printf "* Quartets: $nquart, s: $qsolved, d: $qdifferent, r1: $qr1, r2: $qr2, u: $qUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n\n", $qDC, $qEA, $qSJA, $qSSJA, $qSD;
			}
		}

		
		
		
		if ( ($method eq 'all') || ($method eq 'triplets') ) {
			($ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &tripletsDistance(@splittrees, $units);
			
			
			if ($print_prev eq 'n') {
			}else{
				printf "* Triplets: $ntripl, s: $tsolved, d: $tdifferent, r1: $tr1, r2: $tr2, u: $tUnresolved, DC: %6.6f, EA: %6.6f, SJA: %6.6f, SSJA: %6.6f, SD: %6.6f\n\n", $tDC, $tEA, $tSJA, $tSSJA, $tSD;
			}
		}
	}
	
	
	
	return ($topd, $topdUnpruned, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD);
}

sub elimTaxa() {
	my ($SplitDiste, $level, $treeElim1, $treeElim2, $print_prev) = @_;

	my $num_sp_elim  =0;
	
	
	my @e_taxa_tree1 = sort split '\W+', $treeElim1;
	my $e_n_taxa_tree1 = (@e_taxa_tree1 - 1);
	

	
	my @e_taxa_tree2 = sort split '\W+', $treeElim2;
	my $e_n_taxa_tree2 = (@e_taxa_tree2 - 1);
	


	
	
	my $min = $SplitDiste;
	my $ant = "";
	my @elim_e = ();
	my @elim_e1 = ();
	my %sp_feta = ();
	while ($min > 0) {
		
		$num_sp_elim = @elim_e;
		$num_sp_elim++;
		
		if ($e_n_taxa_tree1 - $num_sp_elim < 4) {
			if (@elim_e[0] == 0) {
				$min = 0;
				@elim_e = @e_taxa_tree1;
			}
			last;
		}
		foreach $te1(@e_taxa_tree1) {
			if ( ($te1 ne "")  && ($sp_feta{$te1} eq "") ){
				
				@auxelimsp = ();
				push @auxelimsp, $te1;
				push @auxelimsp, @elim_e;
				my $epTree1 = &pruneTree($treeElim1, @auxelimsp);
				my $epTree2 = &pruneTree($treeElim2, @auxelimsp);
				$epTree1 = &unRootTrees($epTree1);
				$epTree2 = &unRootTrees($epTree2);
				my @splittreesE = ($epTree1, $epTree2);

				
				($AUXSplitDist, $AUXp, $AUXm) = &splitDistance(@splittreesE);
				

				
				if ($AUXSplitDist <= $min) {
					if ($AUXSplitDist < $min) {
						$min = $AUXSplitDist;
						foreach $aaa(@elim_e1) {
							$sp_feta{$aaa} = "";
						}
						@elim_e1 = ();
						if ($sp_feta{$te1} eq "") {
							$sp_feta{$te1} = "Y";
							push @elim_e1, $te1;
						}
					}
					elsif ($AUXSplitDist == $min) {
						if ($sp_feta{$te1} eq "") {
							$sp_feta{$te1} = "Y";
							push @elim_e1, $te1;
						}
					}
					
				}
			}
		}

		
		push @elim_e, @elim_e1;
		

		if ($min == $ant) {
			
			
			last;
		}

		
		$num_sp_elim = @elim_e;
		
		if ($e_n_taxa_tree1 - $num_sp_elim < 4) {
			if (@elim_e[0] == 0) {
				$min = 0;
				@elim_e = @e_taxa_tree1;
			}
			last;
		}

		
		
		@auxelimsp = ();
		push @auxelimsp, @elim_e;
		my $epTree1 = &pruneTree($treeElim1, @auxelimsp);
		my $epTree2 = &pruneTree($treeElim2, @auxelimsp);
		$epTree1 = &unRootTrees($epTree1);
		$epTree2 = &unRootTrees($epTree2);
		my @splittreesE = ($epTree1, $epTree2);

		
		($AUXSplitDist, $AUXp, $AUXm) = &splitDistance(@splittreesE);
		$min = $AUXSplitDist;
		
		if ($AUXSplitDist == 0) {
			last;
		}
		$ant = $min;
	}
	
	if ($print_prev eq 'n') {
	}else{
		print "New Split Distance at level 1: $min\n";
	}

	
	if ( ($min > 0) && ($level > 1)) { 
		$min = $SplitDiste;
		$ant = "";
		@elim_e = ();
		@elim_e1 = ();
		%sp_feta = ();
		while ($min > 0) {
			
			$num_sp_elim = @elim_e;
			$num_sp_elim++;
			$num_sp_elim++;
			
			if ($e_n_taxa_tree1 - $num_sp_elim < 4) {
				if (@elim_e[0] == 0) {
					$min = 0;
					@elim_e = @e_taxa_tree1;
				}
				last;
			}
			foreach $te1(@e_taxa_tree1) {
				foreach  $te2(@e_taxa_tree1) {
					if ( ($te1 ne "")  && ($sp_feta{$te1} eq "") ){
						if ( ($te2 ne "")  && ($sp_feta{$te2} eq "") ){
							if ($te1 ne $te2) {
								
								@auxelimsp = ();
								push @auxelimsp, $te1;
								push @auxelimsp, $te2;
								push @auxelimsp, @elim_e;
								my $epTree1 = &pruneTree($treeElim1, @auxelimsp);
								my $epTree2 = &pruneTree($treeElim2, @auxelimsp);
								$epTree1 = &unRootTrees($epTree1);
								$epTree2 = &unRootTrees($epTree2);
								my @splittreesE = ($epTree1, $epTree2);
				
								
								($AUXSplitDist, $AUXp, $AUXm) = &splitDistance(@splittreesE);
								
				
								
								if ($AUXSplitDist <= $min) {
									if ($AUXSplitDist < $min) {
										$min = $AUXSplitDist;
										foreach $aaa(@elim_e1) {
											$sp_feta{$aaa} = "";
										}
										@elim_e1 = ();
										if ($sp_feta{$te1} eq "") {
											$sp_feta{$te1} = "Y";
											push @elim_e1, $te1;
										}
										if ($sp_feta{$te2} eq "") {
											$sp_feta{$te2} = "Y";
											push @elim_e1, $te2;
										}
									}
									elsif ($AUXSplitDist == $min) {
										if ($sp_feta{$te1} eq "") {
											$sp_feta{$te1} = "Y";
											push @elim_e1, $te1;
										}
										if ($sp_feta{$te2} eq "") {
											$sp_feta{$te2} = "Y";
											push @elim_e1, $te2;
										}
									}
									
								}
							}
						}
					}
				}
			}
	
			
			push @elim_e, @elim_e1;
			
	
			if ($min == $ant) {
				
				
				last;
			}
	
			
			$num_sp_elim = @elim_e;
			
			if ($e_n_taxa_tree1 - $num_sp_elim < 4) {
				if (@elim_e[0] == 0) {
					$min = 0;
					@elim_e = @e_taxa_tree1;
				}
				last;
			}
	
			
			
			@auxelimsp = ();
			push @auxelimsp, @elim_e;
			my $epTree1 = &pruneTree($treeElim1, @auxelimsp);
			my $epTree2 = &pruneTree($treeElim2, @auxelimsp);
			$epTree1 = &unRootTrees($epTree1);
			$epTree2 = &unRootTrees($epTree2);
			my @splittreesE = ($epTree1, $epTree2);
	
			
			($AUXSplitDist, $AUXp, $AUXm) = &splitDistance(@splittreesE);
			$min = $AUXSplitDist;
			
			if ($AUXSplitDist == 0) {
				last;
			}
	
			$ant = $min;
		}

		if ($print_prev eq 'n') {
		}else{
			print "New Split Distance at Level 2: $min\n";
		}
	}
	

	
	if (($min > 0) && ($level > 2)) { 
		$min = $SplitDiste;
		$ant = "";
		@elim_e = ();
		@elim_e1 = ();
		%sp_feta = ();
		while ($min > 0) {
			
			$num_sp_elim = @elim_e;
			$num_sp_elim++;
			$num_sp_elim++;
			$num_sp_elim++;
			
			if ($e_n_taxa_tree1 - $num_sp_elim < 4) {
				if (@elim_e[0] == 0) {
					$min = 0;
					@elim_e = @e_taxa_tree1;
				}
				last;
			}
			foreach $te1(@e_taxa_tree1) {
				foreach  $te2(@e_taxa_tree1) {
					foreach  $te3(@e_taxa_tree1) {
						if ( ($te1 ne "")  && ($sp_feta{$te1} eq "") ){
							if ( ($te2 ne "")  && ($sp_feta{$te2} eq "") ){
								if ( ($te3 ne "")  && ($sp_feta{$te3} eq "") ){
									if ( ($te1 ne $te2) && ($te1 ne $te3) && ($te2 ne $te3) ){
										
										@auxelimsp = ();
										push @auxelimsp, $te1;
										push @auxelimsp, $te2;
										push @auxelimsp, $te3;
										push @auxelimsp, @elim_e;
										my $epTree1 = &pruneTree($treeElim1, @auxelimsp);
										my $epTree2 = &pruneTree($treeElim2, @auxelimsp);
										$epTree1 = &unRootTrees($epTree1);
										$epTree2 = &unRootTrees($epTree2);
										my @splittreesE = ($epTree1, $epTree2);
										
										
										($AUXSplitDist, $AUXp, $AUXm) = &splitDistance(@splittreesE);
										
						
										
										if ($AUXSplitDist <= $min) {
											if ($AUXSplitDist < $min) {
												$min = $AUXSplitDist;
												foreach $aaa(@elim_e1) {
													$sp_feta{$aaa} = "";
												}
												@elim_e1 = ();
												if ($sp_feta{$te1} eq "") {
													$sp_feta{$te1} = "Y";
													push @elim_e1, $te1;
												}
												if ($sp_feta{$te2} eq "") {
													$sp_feta{$te2} = "Y";
													push @elim_e1, $te2;
												}
												if ($sp_feta{$te3} eq "") {
													$sp_feta{$te3} = "Y";
													push @elim_e1, $te3;
												}
											}
											elsif ($AUXSplitDist == $min) {
												if ($sp_feta{$te1} eq "") {
													$sp_feta{$te1} = "Y";
													push @elim_e1, $te1;
												}
												if ($sp_feta{$te2} eq "") {
													$sp_feta{$te2} = "Y";
													push @elim_e1, $te2;
												}
												if ($sp_feta{$te3} eq "") {
													$sp_feta{$te3} = "Y";
													push @elim_e1, $te3;
												}
											}
											
										}
									}
								}
							}
						}
					}
				}
			}
	
			
			push @elim_e, @elim_e1;
		
	
			if ($min == $ant) {
		
	
				last;
			}
	

			$num_sp_elim = @elim_e;


			if ($e_n_taxa_tree1 - $num_sp_elim < 4) {
				if (@elim_e[0] == 0) {
					$min = 0;
					@elim_e = @e_taxa_tree1;
				}
				last;
			}
	
			
			
			@auxelimsp = ();
			push @auxelimsp, @elim_e;

			

			my $epTree1 = &pruneTree($treeElim1, @auxelimsp);
			my $epTree2 = &pruneTree($treeElim2, @auxelimsp);

			

			$epTree1 = &unRootTrees($epTree1);
			$epTree2 = &unRootTrees($epTree2);
			my @splittreesE = ($epTree1, $epTree2);

			

			
			($AUXSplitDist, $AUXp, $AUXm) = &splitDistance(@splittreesE);
			$min = $AUXSplitDist;
			

			if ($AUXSplitDist == 0) {
				last;
			}
	
			$ant = $min;
		}

		if ($print_prev eq 'n') {
		}else{
			print "New Split Distance at Level 3: $min\n";
		}
	}
	

	
	if (($min > 0) && ($level > 3)) { 
		$min = $SplitDiste;
		$ant = "";
		@elim_e = ();
		@elim_e1 = ();
		%sp_feta = ();
		while ($min > 0) {
			
			$num_sp_elim = @elim_e;
			$num_sp_elim++;
			$num_sp_elim++;
			$num_sp_elim++;
			$num_sp_elim++;
			
			if ($e_n_taxa_tree1 - $num_sp_elim < 4) {
				if (@elim_e[0] == 0) {
					$min = 0;
					@elim_e = @e_taxa_tree1;
				}
				last;
			}
			foreach $te1(@e_taxa_tree1) {
				foreach  $te2(@e_taxa_tree1) {
					foreach  $te3(@e_taxa_tree1) {
						foreach  $te4(@e_taxa_tree1) {
							if ( ($te1 ne "")  && ($sp_feta{$te1} eq "") ){
								if ( ($te2 ne "")  && ($sp_feta{$te2} eq "") ){
									if ( ($te3 ne "")  && ($sp_feta{$te3} eq "") ){
										if ( ($te4 ne "")  && ($sp_feta{$te4} eq "") ){
											if ( ($te1 ne $te2) && ($te1 ne $te3) && ($te2 ne $te3) && ($te1 ne $te4) && ($te2 ne $te4) && ($te3 ne $te4)){
												
												@auxelimsp = ();
												push @auxelimsp, $te1;
												push @auxelimsp, $te2;
												push @auxelimsp, $te3;
												push @auxelimsp, $te4;
												push @auxelimsp, @elim_e;
												my $epTree1 = &pruneTree($treeElim1, @auxelimsp);
												my $epTree2 = &pruneTree($treeElim2, @auxelimsp);
												$epTree1 = &unRootTrees($epTree1);
												$epTree2 = &unRootTrees($epTree2);
												my @splittreesE = ($epTree1, $epTree2);
												
												
												($AUXSplitDist, $AUXp, $AUXm) = &splitDistance(@splittreesE);
												
								
												
												if ($AUXSplitDist <= $min) {
													if ($AUXSplitDist < $min) {
														$min = $AUXSplitDist;
														foreach $aaa(@elim_e1) {
															$sp_feta{$aaa} = "";
														}
														@elim_e1 = ();
														if ($sp_feta{$te1} eq "") {
															$sp_feta{$te1} = "Y";
															push @elim_e1, $te1;
														}
														if ($sp_feta{$te2} eq "") {
															$sp_feta{$te2} = "Y";
															push @elim_e1, $te2;
														}
														if ($sp_feta{$te3} eq "") {
															$sp_feta{$te3} = "Y";
															push @elim_e1, $te3;
														}
														if ($sp_feta{$te4} eq "") {
															$sp_feta{$te4} = "Y";
															push @elim_e1, $te4;
														}
													}
													elsif ($AUXSplitDist == $min) {
														if ($sp_feta{$te1} eq "") {
															$sp_feta{$te1} = "Y";
															push @elim_e1, $te1;
														}
														if ($sp_feta{$te2} eq "") {
															$sp_feta{$te2} = "Y";
															push @elim_e1, $te2;
														}
														if ($sp_feta{$te3} eq "") {
															$sp_feta{$te3} = "Y";
															push @elim_e1, $te3;
														}
														if ($sp_feta{$te4} eq "") {
															$sp_feta{$te4} = "Y";
															push @elim_e1, $te4;
														}
													}
													
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
	
			
			push @elim_e, @elim_e1;
			
	
			if ($min == $ant) {
				
				
				last;
			}
	
			
			$num_sp_elim = @elim_e;
			

			if ($e_n_taxa_tree1 - $num_sp_elim < 4) {
				if (@elim_e[0] == 0) {
					$min = 0;
					@elim_e = @e_taxa_tree1;
				}
				last;
			}
	
			
			
			@auxelimsp = ();
			push @auxelimsp, @elim_e;

			

			my $epTree1 = &pruneTree($treeElim1, @auxelimsp);
			my $epTree2 = &pruneTree($treeElim2, @auxelimsp);

			

			$epTree1 = &unRootTrees($epTree1);
			$epTree2 = &unRootTrees($epTree2);
			my @splittreesE = ($epTree1, $epTree2);

			

			
			($AUXSplitDist, $AUXp, $AUXm) = &splitDistance(@splittreesE);
			$min = $AUXSplitDist;
			

			if ($AUXSplitDist == 0) {
				last;
			}
	
			$ant = $min;
		}

		if ($print_prev eq 'n') {
		}else{
			print "New Split distance Level 4: $min\n";
		}
	}
	

	
	my %elim_e_hash = ();
	my @aux_elim_e = ();
	$num_sp_elim = 0;
	foreach my $hola(@elim_e) {
		if ($hola ne '') {
			
			if ($elim_e_hash{$hola} eq 's') {
			}else {
				$elim_e_hash{$hola} = 's';
				push @aux_elim_e, $hola;
				$num_sp_elim++;
			}
		}
	}
	@elim_e = @aux_elim_e;

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
	if (@elim_e == "") {
		@elim_e = "none";
	}

	
	print "\n";
	return ($min, $num_sp_elim, $e_n_taxa_tree1, @elim_e);
}

sub tripletsDistance() {
	my ($treeTriplet1, $treeTriplet2, $units) = @_;
	
	
	
	my @t_taxa_tree1 = sort split '\W+', $treeTriplet1;
	my $t_n_taxa_tree1 = (@t_taxa_tree1 - 1);
	

	
	my @t_taxa_tree2 = sort split '\W+', $treeTriplet2;
	my $t_n_taxa_tree2 = (@t_taxa_tree2 - 1);
	

	my $elim_taxa = $t_n_taxa_tree1 - 4;
	my $number_triplets = $t_n_taxa_tree1*($t_n_taxa_tree1-1)*($t_n_taxa_tree1-2) / 6;
	
	

	
	my $numberoftriplets = 100;
	if ($units eq "all") {
		$ferunits = "all";
	}elsif ($units eq "random") {
		$ferunits = "random";
	}elsif ($units eq "relative") {
		if ($number_triplets <= 100) {
			$ferunits = "all";
		}else {
			$ferunits = "random";
		}
	}

	
	my $tsolved = 0;
	my $tdifferent = 0;
	my $tUnresolved = 0;
	my $tr1 = 0;
	my $tr2 = 0;
	my $ti = 0;
	my $ntripl = 0;
	my %tripletsFets = ();
	if ($ferunits eq "all") {
		foreach $tsp1(sort(@t_taxa_tree1)) {
			foreach $tsp2(sort(@t_taxa_tree1)) {
				foreach $tsp3(sort(@t_taxa_tree1)) {
					
					if ( ($tsp1 ne $tsp2) && ($tsp1 ne $tsp3) && ($tsp2 ne $tsp3) && ($tsp1 ne "") && ($tsp2 ne "") && ($tsp3 ne "") ) {
						my @AuxTriplets = ("$tsp1", "$tsp2", "$tsp3");
						@AuxTriplets = sort @AuxTriplets;
						$tspa1 = @AuxTriplets[0];
						$tspa2 = @AuxTriplets[1];
						$tspa3 = @AuxTriplets[2];
						if ($tripletsFets{"$tspa1-$tspa2-$tspa3"} ne "Y") {
							
							my @telimsp = ();
							foreach $te(@t_taxa_tree1) {
								if ( "$tspa1-$tspa2-$tspa3" !~ /$te/) {
									push @telimsp, $te;
								}
							}
							
							my $tpTree1 = &pruneTree($treeTriplet1, @telimsp);
							my $tpTree2 = &pruneTree($treeTriplet2, @telimsp);
							my $resolvedTree1 = "Y"; 
							my $resolvedTree2 = "Y";
							if ($tpTree1 =~ /^\(\w+,\w+,\w+\);$/) {
								$resolvedTree1 = "N";
							}
							if ($tpTree2 =~ /^\(\w+,\w+,\w+\);$/) {
								$resolvedTree2 = "N";
							}
							
			
							if ( ($resolvedTree1 eq 'N') && ($resolvedTree2 eq 'N') ) {
								
								$tUnresolved++;
							}elsif ( ($resolvedTree1 eq 'Y') && ($resolvedTree2 eq 'N') ) {	
								
								$tr1++;
							}elsif ( ($resolvedTree1 eq 'N') && ($resolvedTree2 eq 'Y') ) {
								
								$tr2++;
							}elsif ( ($resolvedTree1 eq 'Y') && ($resolvedTree2 eq 'Y') ) {
								
								$tpTree1 =~ /\((\w+),(\w+)\)/;
								my $tjunts1a = $1;
								my $tjunts1b = $2;
								if ($tjunts1a gt  $tjunts1b) {
									$tjunts1a = "$tjunts1a,$tjunts1b";
								}else {
									$tjunts1a = "$tjunts1b,$tjunts1a";
								}
	
								
								$tpTree2 =~ /\((\w+),(\w+)\)/;
								my $tjunts2a = $1;
								my $tjunts2b = $2;
								if ($tjunts2a gt  $tjunts2b) {
									$tjunts2a = "$tjunts2a,$tjunts2b";
								}else {
									$tjunts2a = "$tjunts2b,$tjunts2a";
								}
	
								if ($tjunts1a eq $tjunts2a) {
									
									$tsolved++;
								}else {
									
									$tdifferent++;
								}
							}
							$tripletsFets{"$tspa1-$tspa2-$tspa3"} = "Y";
							$ntripl++;
						}
					}
	
				}
			}	
		}
	}
	elsif ($ferunits eq "random") {
		for ($ntripl = 0;$ntripl < $numberoftriplets;$ntripl++) {
			
	
			
			my $a = "";
			my $b = "";
			my $c = "";
			my $d = "";
	
			$a = int(rand($t_n_taxa_tree1))+1;
			$tsp1 = $t_taxa_tree1[$a];
	
			while (($b eq "") || ($tsp2 eq $tsp1)) {
				$b = int(rand($t_n_taxa_tree1))+1;
				$tsp2 = $t_taxa_tree1[$b];
			}
			while ( ($c eq "") || ($tsp3 eq $tsp1) || ($tsp3 eq $tsp2) ) {
				$c = int(rand($t_n_taxa_tree1))+1;
				$tsp3 = $t_taxa_tree1[$c];
			}
			my @AuxTriplets = ("$tsp1", "$tsp2", "$tsp3");
			@AuxTriplets = sort @AuxTriplets;
			$tspa1 = @AuxTriplets[0];
			$tspa2 = @AuxTriplets[1];
			$tspa3 = @AuxTriplets[2];

				
				my @telimsp = ();
				foreach $te(@t_taxa_tree1) {
					if ( "$tspa1-$tspa2-$tspa3" !~ /$te/) {
						push @telimsp, $te;
					}
				}
				
				my $tpTree1 = &pruneTree($treeTriplet1, @telimsp);
				my $tpTree2 = &pruneTree($treeTriplet2, @telimsp);
				my $resolvedTree1 = "Y"; 
				my $resolvedTree2 = "Y";
				if ($tpTree1 =~ /^\(\w+,\w+,\w+\);$/) {
					$resolvedTree1 = "N";
				}
				if ($tpTree2 =~ /^\(\w+,\w+,\w+\);$/) {
					$resolvedTree2 = "N";
				}
				
				
				if ( ($resolvedTree1 eq 'N') && ($resolvedTree2 eq 'N') ) {
					
					$tUnresolved++;
				}elsif ( ($resolvedTree1 eq 'Y') && ($resolvedTree2 eq 'N') ) {	
					
					$tr1++;
				}elsif ( ($resolvedTree1 eq 'N') && ($resolvedTree2 eq 'Y') ) {
					
					$tr2++;
				}elsif ( ($resolvedTree1 eq 'Y') && ($resolvedTree2 eq 'Y') ) {
					
					$tpTree1 =~ /\((\w+),(\w+)\)/;
					my $tjunts1a = $1;
					my $tjunts1b = $2;
					if ($tjunts1a gt  $tjunts1b) {
						$tjunts1a = "$tjunts1a,$tjunts1b";
					}else {
						$tjunts1a = "$tjunts1b,$tjunts1a";
					}

					
					$tpTree2 =~ /\((\w+),(\w+)\)/;
					my $tjunts2a = $1;
					my $tjunts2b = $2;
					if ($tjunts2a gt  $tjunts2b) {
						$tjunts2a = "$tjunts2a,$tjunts2b";
					}else {
						$tjunts2a = "$tjunts2b,$tjunts2a";
					}

					if ($tjunts1a eq $tjunts2a) {
						
						$tsolved++;
					}else {
						$tdifferent++;
					}
				}
				$tripletsFets{"$tspa1-$tspa2-$tspa3"} = "Y";
		}
	}


	

	
	my $tDC = 0;
	if ($ntripl ne "") {
		$tDC = $tdifferent / $ntripl;
	}

	
	my $tEA = 0;
	if ($ntripl ne "") {
		$tEA = ( $tdifferent + $tr1 + $tr2 + $tUnresolved ) / $ntripl;
	}

	
	my $tSJA = 0;
	if ( ($tdifferent + $tsolved) != 0) {
		$tSJA = $tdifferent / ($tdifferent + $tsolved);
	}

	
	my $tSSJA = 0; 
	if ( ($tdifferent + $tsolved + $tUnresolved) != 0) {
		$tSSJA = $tdifferent / ($tdifferent + $tsolved + $tUnresolved);
	}

	
	my $tSD = 0;
	if ( ( 2*$tdifferent + 2*$tsolved  + $tr1 + $tr2 ) != 0 ) {
		$tSD = ( 2*$tdifferent + $tr1 + $tr2 ) / ( 2*$tdifferent + 2*$tsolved  + $tr1 + $tr2 );
	}

	

	return ($ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD);
}	

sub quartetsDistance() {
	my ($treeQuartet1, $treeQuartet2, $units) = @_;
	
	
	
	my @q_taxa_tree1 = sort split '\W+', $treeQuartet1;
	my $q_n_taxa_tree1 = (@q_taxa_tree1 - 1);
	

	
	my @q_taxa_tree2 = sort split '\W+', $treeQuartet2;
	my $q_n_taxa_tree2 = (@q_taxa_tree2 - 1);
	

	my $elim_taxa = $q_n_taxa_tree1 - 4;
	my $number_quartets = $q_n_taxa_tree1*($q_n_taxa_tree1-1)*($q_n_taxa_tree1-2)*($q_n_taxa_tree1-3) / 24;
	

	
	my $numberofquartets = 100;
	if ($units eq "all") {
		$ferunits = "all";
	}elsif ($units eq "random") {
		$ferunits = "random";
	}elsif ($units eq "relative") {
		if ($number_quartets <= 100) {
			$ferunits = "all";
		}else {
			$ferunits = "random";
		}
	}
	

	
	my $qsolved = 0;
	my $qdifferent = 0;
	my $qUnresolved = 0;
	my $qr1 = 0;
	my $qr2 = 0;
	my $qi = 0;
	my $nquart = 0;
	my %quartetsFets = ();
	if ($ferunits eq "all") {
		foreach $qsp1(sort(@q_taxa_tree1)) {
			foreach $qsp2(sort(@q_taxa_tree1)) {
				foreach $qsp3(sort(@q_taxa_tree1)) {
					foreach $qsp4(sort(@q_taxa_tree1)) {
						if ( ($qsp1 ne $qsp2) && ($qsp1 ne $qsp3) && ($qsp1 ne $qsp4) && ($qsp2 ne $qsp3) && ($qsp2 ne $qsp4) && ($qsp3 ne $qsp4) && ($qsp1 ne "") && ($qsp2 ne "") && ($qsp3 ne "") && ($qsp4 ne "")  ) {
							
							
	
							
							
	
							my @AuxQuartets = ("$qsp1", "$qsp2", "$qsp3", "$qsp4");
							@AuxQuartets = sort @AuxQuartets;
	
							$qspa1 = @AuxQuartets[0];
							$qspa2 = @AuxQuartets[1];
	
							$qspa3 = @AuxQuartets[2];
							$qspa4 = @AuxQuartets[3];
	
							if ($quartetsFets{"$qspa1-$qspa2-$qspa3-$qspa4"} ne "Y") {
								
								my @qelimsp = ();
								foreach $qe(@q_taxa_tree1) {
									if ( "$qspa1-$qspa2-$qspa3-$qspa4" !~ /$qe/) {
										push @qelimsp, $qe;
									}
								}
								
								my $qpTree1 = &pruneTree($treeQuartet1, @qelimsp);
								my $qpTree2 = &pruneTree($treeQuartet2, @qelimsp);

								my $resolvedTree1 = "N"; 
								my $resolvedTree2 = "N";
								if ($qpTree1 =~ /^\(\(\w+,\w+\),\(\w+,\w+\)\);$/) {
									$resolvedTree1 = "Y";
								}else {
									if ($qpTree1 =~ s/\((\w+),(\w+)\)//) {
										my $t_x1 = $1;
										my $t_x2 = $2;
										$qpTree1 =~ /(\w+)\W+(\w+)/;
										my $t_x3 = $1;
										my $t_x4 = $2;
										$qpTree1 = "(($t_x1,$t_x2),($t_x3,$t_x4));";
										
										$resolvedTree1 = "Y";
									}
								}
	
								if ($qpTree2 =~ /^\(\(\w+,\w+\),\(\w+,\w+\)\);$/) {
									$resolvedTree2 = "Y";
								}else {
									if ($qpTree2 =~ s/\((\w+),(\w+)\)//) {
										my $t_x1 = $1;
										my $t_x2 = $2;
										$qpTree2 =~ /(\w+)\W+(\w+)/;
										my $t_x3 = $1;
										my $t_x4 = $2;
										$qpTree2 = "(($t_x1,$t_x2),($t_x3,$t_x4));";
										
										$resolvedTree2 = "Y";
										
									}
								}
	
								
								
								
								if ( ($resolvedTree1 eq 'N') && ($resolvedTree2 eq 'N') ) {
									
									$qUnresolved++;
								}elsif ( ($resolvedTree1 eq 'Y') && ($resolvedTree2 eq 'N') ) {	
									
									$qr1++;
								}elsif ( ($resolvedTree1 eq 'N') && ($resolvedTree2 eq 'Y') ) {
									
									$qr2++;
								}elsif ( ($resolvedTree1 eq 'Y') && ($resolvedTree2 eq 'Y') ) {
									
									$qpTree1 =~ /^\(\((\w+),(\w+)\),\((\w+),(\w+)\)\);$/;
									my $qjunts1a = $1;
									my $qjunts1b = $2;
									my $qjunts1c = $3;
									my $qjunts1d = $4;
									if ($qjunts1a gt  $qjunts1b) {
										$qjunts1a = "$qjunts1a,$qjunts1b";
									}else {
										$qjunts1a = "$qjunts1b,$qjunts1a";
									}
									if ($qjunts1c gt  $qjunts1d) {
										$qjunts1b = "$qjunts1c,$qjunts1d";
									}else {
										$qjunts1b = "$qjunts1d,$qjunts1c";
									}
	
									
									$qpTree2 =~ /^\(\((\w+),(\w+)\),\((\w+),(\w+)\)\);$/;
									my $qjunts2a = $1;
									my $qjunts2b = $2;
									my $qjunts2c = $3;
									my $qjunts2d = $4;
									if ($qjunts2a gt  $qjunts2b) {
										$qjunts2a = "$qjunts2a,$qjunts2b";
									}else {
										$qjunts2a = "$qjunts2b,$qjunts2a";
									}
									if ($qjunts2c >  $qjunts2d) {
										$qjunts2b = "$qjunts2c,$qjunts2d";
									}else {
										$qjunts2b = "$qjunts2d,$qjunts2c";
									}
	
									if ( ($qjunts1a eq $qjunts2a) || ($qjunts1a eq $qjunts2b) || ($qjunts1b eq $qjunts2a) || ($qjunts1b eq $qjunts2b) ) {
										
										$qsolved++;
									}else {
										
										$qdifferent++;
									}
								}
								$quartetsFets{"$qspa1-$qspa2-$qspa3-$qspa4"} = "Y";
								$nquart++;
							}
						}
					}
				}
			}	
		}
	}elsif ($ferunits eq "random") {
		for ($nquart = 0;$nquart < $numberofquartets;$nquart++) { 
			
	
			
			my $a = "";
			my $b = "";
			my $c = "";
			my $d = "";
	
			$a = int(rand($q_n_taxa_tree1))+1;
			$qsp1 = $q_taxa_tree1[$a];
	
			while (($b eq "") || ($qsp2 eq $qsp1)) {
				$b = int(rand($q_n_taxa_tree1))+1;
				$qsp2 = $q_taxa_tree1[$b];
			}
			while ( ($c eq "") || ($qsp3 eq $qsp1) || ($qsp3 eq $qsp2) ) {
				$c = int(rand($q_n_taxa_tree1))+1;
				$qsp3 = $q_taxa_tree1[$c];
			}
			while ( ($d eq "") || ($qsp4 eq $qsp1) || ($qsp4 eq $qsp2) || ($qsp4 eq $qsp3) ) {
				$d = int(rand($q_n_taxa_tree1))+1;
				$qsp4 = $q_taxa_tree1[$d];
			}
	
			my @AuxQuartets = ("$qsp1", "$qsp2", "$qsp3", "$qsp4");
			@AuxQuartets = sort @AuxQuartets;
	
			$qspa1 = @AuxQuartets[0];
			$qspa2 = @AuxQuartets[1];
	
			$qspa3 = @AuxQuartets[2];
			$qspa4 = @AuxQuartets[3];
	
			
			my @qelimsp = ();
			foreach $qe(@q_taxa_tree1) {
				if ( "$qspa1-$qspa2-$qspa3-$qspa4" !~ /$qe/) {
					push @qelimsp, $qe;
				}
			}
			
			my $qpTree1 = &pruneTree($treeQuartet1, @qelimsp);
			my $qpTree2 = &pruneTree($treeQuartet2, @qelimsp);

			my $resolvedTree1 = "N"; 
			my $resolvedTree2 = "N";

			if ($qpTree1 =~ /^\(\(\w+,\w+\),\(\w+,\w+\)\);$/) {
				$resolvedTree1 = "Y";
			}else {
				if ($qpTree1 =~ s/\((\w+),(\w+)\)//) {
					my $t_x1 = $1;
					my $t_x2 = $2;
					$qpTree1 =~ /(\w+)\W+(\w+)/;

					my $t_x3 = $1;
					my $t_x4 = $2;

					$qpTree1 = "(($t_x1,$t_x2),($t_x3,$t_x4));";

					
					$resolvedTree1 = "Y";
				}
			}
	
			if ($qpTree2 =~ /^\(\(\w+,\w+\),\(\w+,\w+\)\);$/) {
				$resolvedTree2 = "Y";
			}else {
				if ($qpTree2 =~ s/\((\w+),(\w+)\)//) {
					my $t_x1 = $1;
					my $t_x2 = $2;

					$qpTree2 =~ /(\w+)\W+(\w+)/;

					my $t_x3 = $1;
					my $t_x4 = $2;

					$qpTree2 = "(($t_x1,$t_x2),($t_x3,$t_x4));";

					
					$resolvedTree2 = "Y";
					
				}
			}
	
			
			
			
			if ( ($resolvedTree1 eq 'N') && ($resolvedTree2 eq 'N') ) {
				
				$qUnresolved++;
			}elsif ( ($resolvedTree1 eq 'Y') && ($resolvedTree2 eq 'N') ) {	
				
				$qr1++;
			}elsif ( ($resolvedTree1 eq 'N') && ($resolvedTree2 eq 'Y') ) {
				
				$qr2++;
			}elsif ( ($resolvedTree1 eq 'Y') && ($resolvedTree2 eq 'Y') ) {
				
				$qpTree1 =~ /^\(\((\w+),(\w+)\),\((\w+),(\w+)\)\);$/;

				my $qjunts1a = $1;
				my $qjunts1b = $2;
				my $qjunts1c = $3;
				my $qjunts1d = $4;
				if ($qjunts1a gt  $qjunts1b) {
					$qjunts1a = "$qjunts1a,$qjunts1b";
				}else {
					$qjunts1a = "$qjunts1b,$qjunts1a";
				}
				if ($qjunts1c gt  $qjunts1d) {
					$qjunts1b = "$qjunts1c,$qjunts1d";
				}else {
					$qjunts1b = "$qjunts1d,$qjunts1c";
				}
	
				
				$qpTree2 =~ /^\(\((\w+),(\w+)\),\((\w+),(\w+)\)\);$/;
				my $qjunts2a = $1;
				my $qjunts2b = $2;
				my $qjunts2c = $3;
				my $qjunts2d = $4;
				if ($qjunts2a gt  $qjunts2b) {
					$qjunts2a = "$qjunts2a,$qjunts2b";
				}else {
					$qjunts2a = "$qjunts2b,$qjunts2a";
				}
				if ($qjunts2c gt  $qjunts2d) {
					$qjunts2b = "$qjunts2c,$qjunts2d";
				}else {
					$qjunts2b = "$qjunts2d,$qjunts2c";
				}
	
				if ( ($qjunts1a eq $qjunts2a) || ($qjunts1a eq $qjunts2b) || ($qjunts1b eq $qjunts2a) || ($qjunts1b eq $qjunts2b) ) {
					
					$qsolved++;
				}else {
					
					$qdifferent++;
				}
			}
			
			$quartetsFets{"$qspa1-$qspa2-$qspa3-$qspa4"} = "Y";
		}
	}


	
	
	my $qDC = 0;
	if ($nquart ne "") {
		$qDC = $qdifferent / $nquart;
	}

	
	my $qEA = 0;
	if ($nquart ne "") {
		$qEA = ( $qdifferent + $qr1 + $qr2 + $qUnresolved ) / $nquart;
	}

	
	my $qSJA = 0;
	if ($qdifferent + $qsolved != 0) {
		$qSJA = $qdifferent / ($qdifferent + $qsolved);
	}
	
	
	my $qSSJA = 0;
	if ( ($qdifferent + $qsolved + $qUnresolved) != 0 ) {
		$qSSJA = $qdifferent / ($qdifferent + $qsolved + $qUnresolved);
	}

	
	my $qSD = 0;
	if ( ( 2*$qdifferent + 2*$qsolved  + $qr1 + $qr2 ) != 0 )  {
		$qSD = ( 2*$qdifferent + $qr1 + $qr2 ) / ( 2*$qdifferent + 2*$qsolved  + $qr1 + $qr2 );
	}
	

	return ($nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD);
}	


sub splitDistance() {
	my ($treeSplit1, $treeSplit2) = @_;
	
	
	
	my $Treeaux1 = $treeSplit1;
	my @SPLITS1 = ();
	my %SpLiTs1 = ();
	my $q1 = 0;
	while ($Treeaux1 !~ /\(\w+,\w+,\w+\)/) {
		$Treeaux1 =~ s/\((\w+),(\w+)\)/$1_$2/;

		my $leftSplit = "$1_$2";

		my $rightSplit = $Treeaux1;
		$rightSplit =~ s/$leftSplit//;
		$rightSplit =~ s/_/ /g;
		$rightSplit =~ s/\W+/ /g;

		$leftSplit =~ s/_/ /g;

		my $right1 = "";
		my $left1 = "";
		my @Righ = split '\s+', $rightSplit;
		foreach $r(sort @Righ) {
			if ($r ne "") {
				$right1 .= "$r ";
			}
		}
		my @Left =  split '\s+', $leftSplit;
		foreach $l(sort @Left) {
			if ($l ne "") {
				$left1 .= "$l ";
			}	
		}
		
		if ($left1 ge $right1) {
			$split1 = "$left1"."|"." $right1";
		}else {
			$split1 = "$right1"."|"." $left1";
		}

		push @SPLITS1, $split1;
		$SpLiTs1{$split1} = "y";
		
		$q1++;
	}
	
	

	
	my $Treeaux2 = $treeSplit2;
	my @SPLITS2 = ();
	my %SpLiTs2 = ();
	my $q2 = 0;
	while ($Treeaux2 !~ /\(\w+,\w+,\w+\)/) {
		$Treeaux2 =~ s/\((\w+),(\w+)\)/$1_$2/;

		my $leftSplit = "$1_$2";

		my $rightSplit = $Treeaux2;
		$rightSplit =~ s/$leftSplit//;
		$rightSplit =~ s/_/ /g;
		$rightSplit =~ s/\W+/ /g;

		$leftSplit =~ s/_/ /g;

		my $right2 = "";
		my $left2 = "";
		my @Righ = split '\s+', $rightSplit;
		foreach $r(sort @Righ) {
			if ($r ne "") {
				$right2 .= "$r ";
			}
		}
		my @Left =  split '\s+', $leftSplit;
		foreach $l(sort @Left) {
			if ($l ne "") {
				$left2 .= "$l ";
			}	
		}
		if ($left2 ge $right2) {
			$split2 = "$left2"."|"." $right2";
		}else {
			$split2 = "$right2"."|"." $left2";
		}
		
		push @SPLITS2, $split2;
		$SpLiTs2{$split2} = "y";
		$q2++;
	}
	

	
	my $minq1q2;
	if ($q1 >= $q2) {
		$minq1q2 = $q1;
	}else {
		$minq1q2 = $q2;
	}
	
	
	my $p = 0;
	my $m = 0;
	
	foreach $split(@SPLITS2) {
		if ($SpLiTs1{$split}) {
			
			$p++;
		}else {
			
		}
		$m++;
	}
	$m *=2;
	$p *=2;

	
        
        my $SplitDist = 0;
        if ($m == 0) {
                $SplitDist = 1;
        }else {
                $SplitDist = 1 - ( $p / $m );
        }

	

	return ($SplitDist, $p, $m);
}



sub pruneTree() {
	my ($t, @spNoCom) = @_; 
	
	foreach $nospcom(@spNoCom) {
		
		if ($t =~ /,$nospcom,/) {
			$t =~ s/,$nospcom,/,/;
		}
		elsif ($t =~ /\(\w+,$nospcom\)/) {
			$t =~ s/\((\w+),$nospcom\)/$1/;
		}
		elsif ($t =~ /\($nospcom,\w+\)/) {
			$t =~ s/\($nospcom,(\w+)\)/$1/;
		}
		elsif ($t =~ /^\($nospcom,/) {
			$t =~ s/\($nospcom,/\(/;	
		}
		elsif ($t =~ /,$nospcom\)$/) {
			$t =~ s/,$nospcom\)$/\)/;	
		}
		elsif ($t =~ /\($nospcom,\(/) {
			$t =~ /^(.+)(\($nospcom,\(.+\);)/;
			my $ini = $1;
			my $fi = $2;
			my @TEMP = split '', $fi;
			$contCAR=0;
			$segonTrobat = 0;
			$FI = "false";
			foreach $CAR(@TEMP) {
				if ($segonTrobat < 2) {
					if ($CAR eq "(") {
						$segonTrobat++;
						$cont++;
					}
					if ($segonTrobat < 2) {
						$TEMP[$contCAR] = '';
					}
				}else{
					if ($CAR eq "(") {
						$cont++;
					}elsif ($CAR eq ")") {
						$cont--;
					}
				}
				if ($segonTrobat >= 2) {	
					if ($cont == 0) {
						$TEMP[$contCAR] = '';
						last;
					}
				}
				$contCAR++;
			}
			$fi2="";
			foreach $CAR(@TEMP) {
				$fi2 .= $CAR;
			}
			$t = "$ini"."$fi2";
		}
		elsif ($t =~ /\),$nospcom\)/) {
			$t =~ /(.+\),$nospcom\))(.+)/;
			
			my $ini = $1;
			my $fi = $2;
			$ini = reverse $ini;
			my @TEMP = split '', $ini;
			$contCAR=0;
			$segonTrobat = 0;
			$FI = "false";
			foreach $CAR(@TEMP) {
				if ($segonTrobat < 2) {
					if ($CAR eq ")") {
						$segonTrobat++;
						$cont++;
					}
					if ($segonTrobat < 2) {
						$TEMP[$contCAR] = '';
					}
				}else{
					if ($CAR eq ")") {
						$cont++;
					}elsif ($CAR eq "(") {
						$cont--;
					}
				}
				if ($segonTrobat >= 2) {	
					if ($cont == 0) {
						$TEMP[$contCAR] = '';
						last;
					}
				}
				$contCAR++;
			}
			$ini ="";
			foreach $CAR(@TEMP) {
				$ini .= $CAR;
			}
			$ini =reverse $ini;
			$t = "$ini"."$fi";
		}
		elsif ($t =~ /\w+,$nospcom\)/) {
			
			$t =~ s/(\w+)(,$nospcom)\)/$1\)/;
			
			
		}
		elsif ($t =~ /\($nospcom,\w+/) {
			
			$t =~ s/\(($nospcom,)(\w+)/\($2/;
			
			
		}
		elsif ($t =~ /(\($nospcom\),)/) {
			$t =~ s/\($nospcom\),//;
		}
		elsif ($t =~ /(,\($nospcom\))/) {
			$t =~ s/,\($nospcom\)//;
		}
		elsif ($t =~ /\W$nospcom\W/) {
			print "c---> $t\n";
			print "### Error ### $nospcom\n";
			exit;
		} 

		if ($t =~ /^\(.+;$/) {
		}else {
			$t = "("."$t";
			$t =~ s/;/\);/;
		}

		if ($t =~ /^.+\);$/) {
		}else {
			$t = "("."$t";
			$t =~ s/;/\);/;
		}
		
	}

	
	my @TEMP = split '', $t;
	$contCAR = 0;
	
	
	foreach $CAR(@TEMP) {
		if ($CAR eq "(") {
			$contPAR++;
		}elsif ($CAR eq ")") {
			$contPAR--;
		}
		if ($contPAR == 0) {
			$contCAR++;
			last;
		}
		else {
			$contCAR++;
		}
	}
	if ($TEMP[$contCAR] ne ";") {
		$provT = "("."$t";
		$t = $provT;
		$t =~ s/;/\);/;
	}
	
	

	
	my $t_stax = $t;
	$t_stax =~ s/\w+/*/g;
	$t_stax =~ s/;//;
	
	while ($t_stax =~ s/\(\*,\*\)/\*/) {
		
	}	
	
	
	
	$t =~ s/\((\w+)\)/$1/g;
	

	
	if ($t_stax eq "*") {
		
	}elsif ($t_stax eq "(*,*,*)") {
		
	}else {
		
		$tprova = $t;
		my @TEMP = split '', $t;
		$contPar = $t =~ s/\(/\(/g;
		$anterior = "";
		for ($k=1;$k<=$contPar;$k++) {
			$cont_i = 1;
			$contCAR = 0;
			$start_cont = "n";
			$CONTCARMIG = 0;
		
			foreach $CAR(@TEMP) {
				if ($CAR eq "(") {
					if ($k eq $cont_i) {
						
						$start_cont = "y";
						$contCAR++;
					}else {
						$cont_i++;
					}
					
				}elsif ($CAR eq ")") {
					if($start_cont eq "y") {
						$contCAR--;
					}
					if ($contCAR == 0) {
						$start_cont = "n";
					}
				}
				elsif($start_cont eq "y") {
					$CONTCARMIG++;
				}
			}
			
			$start_cont = "n";
			$contCAR=0;
			$cont_i = 1;
			if ($anterior eq $CONTCARMIG) {
				$t_aux = "";
				$primerTrobat = 'n';
				
				foreach $CAR(@TEMP) {
					if ($CAR eq "(") {
						if ($k == $cont_i) {
							$cont_i++;
							$contCAR++;
						}else {
							$cont_i++;
							$t_aux .= $CAR;
							if ($cont_i > $k) {
								$contCAR++;
							}
						}
					}
					elsif ($CAR eq ")") {
						$contCAR--;
						if ($contCAR == 0) {
							if ($primerTrobat eq 'n') {
								$primerTrobat = 's';
							}else {
								$t_aux .= $CAR;
							}
						}else {
							$t_aux .= $CAR;
						}
					}
					else {
						$t_aux .= $CAR;
					}
				}
				$t = $t_aux;
			}
			
			
			
			$anterior = $CONTCARMIG;
		}
		
		
	}

	
	my $t_stax = $t;
	$t_stax =~ s/\w+/*/g;
	$t_stax =~ s/;//;
	while ($t_stax =~ s/\(\*,\*\)/\*/) {
	}
	
	
	if ($t_stax =~ /\(+\*\)+/) {
		$nombreExtra = $t_stax =~ s/\(/\(/g;
		while ($nombreExtra > 0) {
			
			$t =~ s/^\(//;
			$t =~ s/\);$/;/;
			$nombreExtra--;
		}
		
		
	}elsif ($t_stax =~ /\(\(+\*,\*,\*\)\)+/) {
		$nombreExtra = $t_stax =~ s/\(/\(/g;
		while ($nombreExtra > 0) {
			
			$t =~ s/^\(//;
			$t =~ s/\);$/;/;
			$nombreExtra--;
		}
		
		
	}
	

	
	if ($t =~ /^\(.+;$/) {
	}else {
		$t = "("."$t";
		$t =~ s/;/\);/;
	}
	if ($t =~ /^.+\);$/) {
	}else {
		$t = "("."$t";
		$t =~ s/;/\);/;
	}
	

	return $t;	
}


sub unRootTrees() {
	my $t = shift;
	
	
	my $t_stax = $t;
	$t_stax =~ s/\w+/*/g;
	$t_stax =~ s/;//;
	
	
	while ($t_stax =~ s/\(\*,\*\)/\*/) {
		
	}
	
	if ($t_stax eq "*") {
		
		
		$t = &fRootToUnroot($t);
	}elsif ($t_stax eq "(*,*,*)") {
		
		
		
	}else {
		
		
		sleep(3);
		die "# ERROR # (a) Unrooting the trees --> $t_stax \n"
	}

	return $t;
}

sub fRootToUnroot () {
	my $t = shift;
	
	if ($t =~ /^\((\w+),(.+)\);/){
		$ini = $1;
		$fi = $2;
		$fi =~ s/^\(//;
		$fi =~ s/\)$//;
		$t = "("."$ini".","."$fi".");";
		
	}elsif ($t =~ /\((.+),(\w+)\);$/){
		$ini = $1;
		$fi = $2;
		$ini =~ s/^\(//;
		$ini =~ s/\)$//;
		$t = "("."$ini".","."$fi".");";
		
	}
	elsif ($t =~ /^\((\(.+\)),(\(.+\))\);$/) {
		my $t_aux = $t;
		$t_aux =~ s/^\(//;
		$t_aux =~ s/\);$//;
		$t_aux =~ /^(\(.+\)),(\(.+\))$/;
		@TEMP = split //, $t_aux;
		$cont_par = 0;
		$trobat = "false";
		$ini = "";
		$fi="";
		foreach $CAR(@TEMP) {
			
			if ($trobat eq "false") {
				$ini .= $CAR;
			}else {
				$fi .= $CAR;
			}

			if ($CAR eq "(") {
				$cont_par++;
			}
			elsif ($CAR eq ")") {
				$cont_par--;
				if ($cont_par == 0) {
					$trobat = "true";
				}
			}
		}
		$fi =~ s/^,//;
		

		
		$sp_ini = $ini =~ s/(\w+)/$1/g;
		$sp_fi = $fi =~ s/(\w+)/$1/g;
		
		if ($sp_ini >= $sp_fi) {
			$fi =~ s/^\(//;
			$fi =~ s/\)$//;
			
		}else {
			$ini =~ s/^\(//;
			$ini =~ s/\)$//;
			
		}
		$t = "("."$ini".","."$fi".");";
		
	}
	else {
		
		die "# ERROR # (b) Unrooting the trees --> $t_stax \n"
	}
	return $t;
}



sub MoS() {
	my ($name1_IN, $tree1_IN, $name2_IN, $tree2_IN) = @_;
	my %trees = ();	
	$trees{$name1_IN} = $tree1_IN;
	$trees{$name2_IN} = $tree2_IN;

	@treesINput = ($name1_IN, $name2_IN);
	my $Nsp1 = 0;
	my $family1 = 'S';
	my $Nsp2 = 0;
	my $family2 = "S";
	
	
	my %sp1 = ();
	my %sp2=();
	my $i=1;
	foreach $tree(@treesINput) {
		$treeSP = $trees{$tree};
		while ($treeSP =~ /(\w+)/) {
			$sp = $1;
			$treeSP =~ s/$sp//;
			if ($i == 1) {
				$sp1{$sp}++;	
				$Nsp1++;
				if ($sp1{$sp} > 1) {
					$family1 = "M";
				}
			}elsif ($i == 2) {
				$sp2{$sp}++;
				$Nsp2++;
				if ($sp2{$sp} > 1) {
					$family2 = "M";
				}
			}
		}
		$i=2;
	}
	
	
	
	
	
	
	
	
	
	
	
	return ($Nsp1, $family1, $Nsp2, $family2);
}


sub topdSM() {
	my ($name1_INsm, $tree1_INsm, $name2_INsm, $tree2_INsm, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level,$units) = @_;
	my %trees = ();
	my $topdSM = 1;

	
	$trees{$name1_INsm} = $tree1_INsm;

	
	$trees{$name2_INsm} = $tree2_INsm;
	my @SingleGeneTrees = &fMtS($name2_INsm, $trees{$name2_INsm}, $nombreSGTgenerats, $SingleMultiple, $print_prev);

	my $contSingTrees = 0;
	my @topds = ();
	my @topdsUnpruned =();
	my $topdSM = 0;
	my $topdSMunpruned = 0;
	my ($SplitDistSM, $SplitDistSMsd, $qSM, $mSM, $qSMsd, $mSMsd) = ();
	my @SPLITS = ();
	my @DIFFER =();
	my @TOTAL = ();
	my $SplitDistSM = 0;

	
	my $num_sp_elimSM = 0;
	my $e_n_taxa_tree1SM= 0;
	my $minSM = 0;
	my $num_sp_elimSMsd = 0;
	my $e_n_taxa_tree1SMsd= 0;

	my $minSMsd = 0;
	my $dif_elem = "";
	my %DisTaxSM = ();

	
  	my $nquartSM = 0; 
	my $qsolvedSM = 0;
	my $qdifferentSM = 0; 
	my $qr2SM = 0; 
	my $qr1SM = 0; 
	my $qUnresolvedSM = 0;
	my $qDCSM = 0;
	my $qEASM = 0;
	my $qSJASM = 0;
	my $qSSJASM = 0;
	my $qSDSM = 0;

  	my $nquartSMsd = 0; 
	my $qsolvedSMsd = 0;
	my $qdifferentSMsd = 0; 
	my $qr2SMsd = 0; 
	my $qr1SMsd = 0; 
	my $qUnresolvedSMsd = 0;
	my $qDCSMsd = 0;
	my $qEASMsd = 0;
	my $qSJASMsd = 0;
	my $qSSJASMsd = 0;
	my $qSDSMsd = 0;

	
	my $ntriplSM =  0;
	my $tsolvedSM =  0;
	my $tdifferentSM =  0;
	my $tr1SM =  0;
	my $tr2SM =  0;
	my $tUnresolvedSM =  0; 
	my $tDCSM =  0; 
	my $tEASM =  0;
	my $tSJASM =  0;
	my $tSSJASM =  0; 
	my $tSDSM = 0;	

	my $ntriplSMsd =  0;
	my $tsolvedSMsd =  0;
	my $tdifferentSMsd =  0;
	my $tr1SMsd =  0;
	my $tr2SMsd =  0;
	my $tUnresolvedSMsd =  0; 
	my $tDCSMsd =  0; 
	my $tEASMsd =  0;
	my $tSJASMsd =  0;
	my $tSSJASMsd =  0; 
	my $tSDSMsd = 0;	

	
	my @Array_num_sp_elimSM = ();
	my @Array_e_n_taxa_tree1SM = ();
	my @Array_minSM = ();
	
	
  	my @Array_nquartSM = ();
	my @Array_qsolvedSM = ();
	my @Array_qdifferentSM = ();
	my @Array_qr2SM = ();
	my @Array_qr1SM = ();
	my @Array_qUnresolvedSM = ();
	my @Array_qDCSM = ();
	my @Array_qEASM = ();
	my @Array_qSJASM = ();
	my @Array_qSSJASM = ();
	my @Array_qSDSM = ();

	
	my @Array_ntriplSM = ();
	my @Array_tsolvedSM = ();
	my @Array_tdifferentSM = ();
	my @Array_tr1SM = ();
	my @Array_tr2SM = ();
	my @Array_tUnresolvedSM = ();
	my @Array_tDCSM = ();
	my @Array_tEASM = ();
	my @Array_tSJASM = ();
	my @Array_tSSJASM = ();
	my @Array_tSDSM = ();

	foreach $singTree(@SingleGeneTrees) {
		$nomtreesingle = "$name2_INsm"."_$contSingTrees";
		($topd, $topdUnpruned, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($name1_INsm, $tree1_INsm, $nomtreesingle, $singTree, $print_prev, $method, $level, $units);

		
		$topdSM += $topd;

		
		$SplitDistSM += $SplitDist;
		$qSM += $q;
		$mSM += $m;

		
		$num_sp_elimSM += $num_sp_elim;
		$e_n_taxa_tree1SM += $e_n_taxa_tree1;
		$minSM += $min;
		$DisTaxSM{$DisTax}++;
		
		
  		$nquartSM += $nquart; 
		$qsolvedSM += $qsolved;
		$qdifferentSM += $qdifferent; 
		$qr2SM += $qr2; 
		$qr1SM += $qr1; 
		$qUnresolvedSM += $qUnresolved;
		$qDCSM += $qDC;
		$qEASM += $qEA;
		$qSJASM += $qSJA;
		$qSSJASM += $qSSJA;
		$qSDSM += $qSD;

		
		$ntriplSM +=  $ntripl;
		$tsolvedSM +=  $tsolved;
		$tdifferentSM +=  $tdifferent;
		$tr1SM +=  $tr1;
		$tr2SM +=  $tr2;
		$tUnresolvedSM +=  $tUnresolved; 
		$tDCSM +=  $tDC; 
		$tEASM +=  $tEA;
		$tSJASM +=  $tSJA;
		$tSSJASM +=  $tSSJA; 
		$tSDSM += $tSD;	

		
		push @SPLITS, $SplitDist;
		push @DIFFER, $q;
		push @TOTAL, $m;

		
		push @topds, $topd;
		$topdSMunpruned += $topdUnpruned;
		push @topdsUnpruned, $topdUnpruned;

		
		push @Array_num_sp_elimSM, $num_sp_elim;
		push @Array_e_n_taxa_tree1SM, $e_n_taxa_tree1;
		push @Array_minSM, $min;

		
  		push @Array_nquartSM, $nquart; 
		push @Array_qsolvedSM, $qsolved;
		push @Array_qdifferentSM, $qdifferent; 
		push @Array_qr2SM, $qr2; 
		push @Array_qr1SM, $qr1; 
		push @Array_qUnresolvedSM, $qUnresolved;
		push @Array_qDCSM, $qDC;
		push @Array_qEASM, $qEA;
		push @Array_qSJASM, $qSJA;
		push @Array_qSSJASM, $qSSJA;
		push @Array_qSDSM, $qSD;

		
		push @Array_ntriplSM, $ntripl;
		push @Array_tsolvedSM, $tsolved;
		push @Array_tdifferentSM, $tdifferent;
		push @Array_tr1SM, $tr1;
		push @Array_tr2SM, $tr2;
		push @Array_tUnresolvedSM, $tUnresolved; 
		push @Array_tDCSM, $tDC; 
		push @Array_tEASM, $tEA;
		push @Array_tSJASM, $tSJA;
		push @Array_tSSJASM, $tSSJA; 
		push @Array_tSDSM, $tSD;

		$contSingTrees++;
	}

	my $SD = 0;
	my $SDunpruned = 0;
	if ($com <= 3) {
		($topdSM, $SD, $topdSMunpruned, $SDunpruned) = (0,0,0,0);
		if ($print_prev eq 'n') {
		}else{
			if ( ($method eq 'all') || ($method eq 'nodal') ) {
				printf "* Nodal Distance SM (Prunde/Unpruned): (%6.6f +/-%6.6f ) / (%6.6f +/-%6.6f )\n", $topdSM, $SD, $topdSMunpruned, $SDunpruned;
			}
			if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
				print "* Split Distance SM [differents/possibles]: 0 +/- 0.000 [ 0 +/- 0.000 \/ 0 +/- 0.000 ]\n";
			}
			if ( ($method eq 'all') || ($method eq 'disagree') ) {
				print "* Disagreement [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: [none]\n";
			}
			if ( ($method eq 'all') || ($method eq 'quartets') ) {
				print "* Quartets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
			}
			if ( ($method eq 'all') || ($method eq 'triplets') ) {
				print "* Triplets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
			}
		}
	}else {
		
		$topdSM /= $contSingTrees;
		$topdSMunpruned /= $contSingTrees;	
	
		
		$SD = 0;
		foreach $topdsi(@topds) {
			$SD += ($topdsi - $topdSM)**2;
		}
		$SD = sqrt ( $SD / $contSingTrees );

		$SDunpruned = 0;
		foreach $topdsuni(@topdsUnpruned) {
			$SDunpruned += ($topdsuni - $topdSMunpruned)**2;
		}
		$SDunpruned = sqrt ( $SDunpruned / $contSingTrees );

		
		$SplitDistSM /= $contSingTrees;
		$qSM /= $contSingTrees;
		$mSM /= $contSingTrees;

		
		$SplitDistSMsd = 0;
		foreach $kk(@SPLITS) {
			$SplitDistSMsd += ($kk - $SplitDistSM)**2;
		}
		$SplitDistSMsd = sqrt ( $SplitDistSMsd / $contSingTrees );

		$qSMsd = 0;
		foreach $kk(@DIFFER) {
			$qSMsd += ($kk - $qSM)**2;
		}
		$qSMsd = sqrt ( $qSMsd / $contSingTrees );

		$mSMsd = 0;
		foreach $kk(@TOTAL) {
			$mSMsd += ($kk - $mSM)**2;
		}
		$mSMsd = sqrt ( $mSMsd / $contSingTrees );

	
		
		$num_sp_elimSM  /= $contSingTrees;
		$e_n_taxa_tree1SM  /= $contSingTrees;
		$minSM  /= $contSingTrees;

		foreach $kk(@Array_num_sp_elimSM) {
			$num_sp_elimSMsd += ($kk - $num_sp_elimSM)**2;
		}
		$num_sp_elimSMsd = sqrt ($num_sp_elimSMsd/$contSingTrees);

		foreach $kk(@Array_e_n_taxa_tree1SM) {
			$e_n_taxa_tree1SMsd += ($kk - $e_n_taxa_tree1SM)**2;
		}
		$e_n_taxa_tree1SMsd = sqrt ($e_n_taxa_tree1SMsd/$contSingTrees);

		foreach $kk(@Array_minSM) {
			$minSMsd += ($kk - $minSM)**2;
		}
		$minSMsd = sqrt ($minSMsd/$contSingTrees);

		foreach $DisTaxSMelement(sort keys(%DisTaxSM)) {
			$aux_dist__ = $DisTaxSM{$DisTaxSMelement};
			$DisTaxSMelement =~ s/-/ /g;
			$dif_elem .= "[$DisTaxSMelement("."$aux_dist__".")] ";
			
		}

		
  		$nquartSM  /= $contSingTrees;
		$qsolvedSM  /= $contSingTrees;
		$qdifferentSM  /= $contSingTrees;
		$qr2SM  /= $contSingTrees;
		$qr1SM  /= $contSingTrees;
		$qUnresolvedSM  /= $contSingTrees;
		$qDCSM  /= $contSingTrees;
		$qEASM  /= $contSingTrees;
		$qSJASM  /= $contSingTrees;
		$qSSJASM  /= $contSingTrees;
		$qSDSM  /= $contSingTrees;

		foreach $kk(@Array_nquartSM) {
			$nquartSMsd += ($kk - $nquartSM)**2;
		}
		$nquartSMsd = sqrt ($nquartSMsd/$contSingTrees);

		foreach $kk(@Array_qsolvedSM) {
			$qsolvedSMsd += ($kk - $qsolvedSM)**2;
		}
		$qsolvedSMsd = sqrt ($qsolvedSMsd/$contSingTrees);

		foreach $kk(@Array_qdifferentSM) {
			$qdifferentSMsd += ($kk - $qdifferentSM)**2;
		}
		$qdifferentSMsd = sqrt ($qdifferentSMsd/$contSingTrees);

		foreach $kk(@Array_qr2SM) {
			$qr2SMsd += ($kk - $qr2SM)**2;
		}
		$qr2SMsd = sqrt ($qr2SMsd/$contSingTrees);

		foreach $kk(@Array_qr1SM) {
			$qr1SMsd += ($kk - $qr1SM)**2;
		}
		$qr1SMsd = sqrt ($qr1SMsd/$contSingTrees);

		foreach $kk(@Array_qUnresolvedSM) {
			$qUnresolvedSMsd += ($kk - $qUnresolvedSM)**2;
		}
		$qUnresolvedSMsd = sqrt ($qUnresolvedSMsd/$contSingTrees);

		foreach $kk(@Array_qDCSM) {
			$qDCSMsd += ($kk - $qDCSM)**2;
		}	
		$qDCSMsd = sqrt ($qDCSMsd/$contSingTrees);

		foreach $kk(@Array_qEASM) {
			$qEASMsd += ($kk - $qEASM)**2;
		}
		$qEASMsd = sqrt ($qEASMsd/$contSingTrees);

		foreach $kk(@Array_qSJASM) {
			$qSJASMsd += ($kk - $qSJASM)**2;
		}
		$qSJASMsd = sqrt ($qSJASMsd/$contSingTrees);
		
		foreach $kk(@Array_qSSJASM) {
			$qSSJASMsd += ($kk - $qSSJASM)**2;
		}
		$qSSJASMsd = sqrt ($qSSJASMsd/$contSingTrees);

		foreach $kk(@Array_qSDSM) {
			$qSDSMsd += ($kk - $qSDSM)**2;
		}
		$qSDSMsd = sqrt ($qSDSMsd/$contSingTrees);

		
		$ntriplSM  /= $contSingTrees;
		$tsolvedSM  /= $contSingTrees;
		$tdifferentSM  /= $contSingTrees;
		$tr1SM  /= $contSingTrees;
		$tr2SM  /= $contSingTrees;
		$tUnresolvedSM  /= $contSingTrees;
		$tDCSM  /= $contSingTrees;
		$tEASM  /= $contSingTrees;
		$tSJASM  /= $contSingTrees;
		$tSSJASM  /= $contSingTrees;
		$tSDSM  /= $contSingTrees;

		foreach $kk(@Array_ntriplSM) {
			$ntriplSMsd += ($kk - $ntriplSM)**2;
		}
		$ntriplSMsd = sqrt ($ntriplSMsd/$contSingTrees);

		foreach $kk(@Array_tsolvedSM) {
			$tsolvedSMsd += ($kk - $tsolvedSM)**2;
		}
		$tsolvedSMsd = sqrt ($tsolvedSMsd/$contSingTrees);

		foreach $kk(@Array_tdifferentSM) {
			$tdifferentSMsd += ($kk - $tdifferentSM)**2;
		}
		$tdifferentSMsd = sqrt ($tdifferentSMsd/$contSingTrees);

		foreach $kk(@Array_tr1SM) {
			$tr1SMsd += ($kk - $tr1SM)**2;
		}
		$tr1SMsd = sqrt ($tr1SMsd/$contSingTrees);

		foreach $kk(@Array_tr2SM) {
			$tr2SMsd += ($kk - $tr2SM)**2;
		}
		$tr2SMsd = sqrt ($tr2SMsd/$contSingTrees);

		foreach $kk(@Array_tUnresolvedSM) {
			$tUnresolvedSMsd += ($kk - $tUnresolvedSM)**2;
		}
		$tUnresolvedSMsd = sqrt ($tUnresolvedSMsd/$contSingTrees);

		foreach $kk(@Array_tDCSM) {
			$tDCSMsd += ($kk - $tDCSM)**2;
		}
		$tDCSMsd = sqrt ($tDCSMsd/$contSingTrees);

		foreach $kk(@Array_tEASM) {
			$tEASMsd += ($kk - $tEASM)**2;
		}
		$tEASMsd = sqrt ($tEASMsd/$contSingTrees);

		foreach $kk(@Array_tSJASM) {
			$tSJASMsd += ($kk - $tSJASM)**2;
		}
		$tSJASMsd = sqrt ($tSJASMsd/$contSingTrees);

		foreach $kk(@Array_tSSJASM) {
			$tSSJASMsd += ($kk - $tSSJASM)**2;
		}
		$tSSJASMsd = sqrt ($tSSJASMsd/$contSingTrees);

		foreach $kk(@Array_tSDSM) {
			$tSDSMsd += ($kk - $tSDSM)**2;
		}
		$tSDSMsd = sqrt ($tSDSMsd/$contSingTrees);
		
		if ($print_prev eq 'n') {
		}else{
			if ( ($method eq 'all') || ($method eq 'nodal') ) {
				printf "* Nodal Distance SM (Prunde/Unpruned): (%6.6f +/- %6.6f ) / (%6.6f +/- %6.6f )\n", $topdSM, $SD, $topdSMunpruned, $SDunpruned;
			}
			if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
				printf "* Split Distance SM [differents/possibles]: $SplitDistSM +/- %5.3f [ $qSM +/- %5.3f \/ $mSM +/- %5.3f ]\n", $SplitDistSMsd, $qSMsd, $mSMsd ;
			}
			if ( ($method eq 'all') || ($method eq 'disagree') ) {
				printf "* Disagreement [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd ;
			}
			if ( ($method eq 'all') || ($method eq 'quartets') ) {
				printf "* Quartets: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd;
			}
			if ( ($method eq 'all') || ($method eq 'triplets') ) {
				printf "* Triplets: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd;
			}
		}
	}
	
	
	return ($topdSM, $SD, $topdSMunpruned, $SDunpruned, $com, $pcCom2, $SplitDistSM, $SplitDistSMsd, $qSM, $qSMsd, $mSM, $mSMsd, $dif_elem, $num_sp_elimSM, $num_sp_elimSMsd, $e_n_taxa_tree1SM, $e_n_taxa_tree1SMsd, $minSM, $minSMsd, $nquartSM, $nquartSMsd, $qsolvedSM, $qsolvedSMsd, $qdifferentSM, $qdifferentSMsd, $qr2SM, $qr2SMsd, $qr1SM, $qr1SMsd, $qUnresolvedSM, $qUnresolvedSMsd, $qDCSM, $qDCSMsd, $qEASM, $qEASMsd, $qSJASM, $qSJASMsd, $qSSJASM, $qSSJASMsd, $qSDSM, $qSDSMsd, $ntriplSM, $ntriplSMsd, $tsolvedSM, $tsolvedSMsd, $tdifferentSM, $tdifferentSMsd,$tr2SM, $tr2SMsd, $tr1SM, $tr1SMsd, $tUnresolvedSM, $tUnresolvedSMsd, $tDCSM, $tDCSMsd, $tEASM, $tEASMsd, $tSJASM, $tSJASMsd, $tSSJASM, $tSSJASMsd, $tSDSM, $tSDSMsd);
}


sub topdMM() {
	my ($name1_INmm, $tree1_INmm, $name2_INmm, $tree2_INmm, $nombreSGTgenerats, $SingleMultiple, $print_prev, $method, $level,$units) = @_;
	my %trees = ();	
	$trees{$name1_INmm} = $tree1_INmm;
	$trees{$name2_INmm} = $tree2_INmm;

	my $topdMM = 1;
	
	
	
	$trees{$name1_INmm} = $tree1_INmm;
	my @SingleGeneTrees1 = &fMtS($name1_INmm, $trees{$name1_INmm}, $nombreSGTgenerats, $SingleMultiple, $print_prev);
	
	
	$trees{$name2_INmm} = $tree2_INmm;
	my @SingleGeneTrees2 = &fMtS($name2_INmm, $trees{$name2_INmm}, $nombreSGTgenerats, $SingleMultiple, $print_prev);

	
	my $contSingTrees1=0;
	my @topds = ();
	my @topdsUnpruned = ();
	my $topdMMunpruned =0;
	my $topdMM = 0;
	my ($SplitDistMM, $SplitDistMMsd, $qMM, $mMM, $qMMsd, $mMMsd) = (0,0,0);
	my @SPLITS = ();
	my @DIFFER =();
	my @TOTAL = ();

	
	my $num_sp_elimMM = 0;
	my $e_n_taxa_tree1MM= 0;
	my $minMM = 0;
	my $num_sp_elimMMsd = 0;
	my $e_n_taxa_tree1MMsd= 0;

	my $minMMsd = 0;
	my $dif_elem = "";
	my %DisTaxMM = ();

	
  	my $nquartMM = 0; 
	my $qsolvedMM = 0;
	my $qdifferentMM = 0; 
	my $qr2MM = 0; 
	my $qr1MM = 0; 
	my $qUnresolvedMM = 0;
	my $qDCMM = 0;
	my $qEAMM = 0;
	my $qSJAMM = 0;
	my $qSSJAMM = 0;
	my $qSDMM = 0;

  	my $nquartMMsd = 0; 
	my $qsolvedMMsd = 0;
	my $qdifferentMMsd = 0; 
	my $qr2MMsd = 0; 
	my $qr1MMsd = 0; 
	my $qUnresolvedMMsd = 0;
	my $qDCMMsd = 0;
	my $qEAMMsd = 0;
	my $qSJAMMsd = 0;
	my $qSSJAMMsd = 0;
	my $qSDMMsd = 0;

	
	my $ntriplMM =  0;
	my $tsolvedMM =  0;
	my $tdifferentMM =  0;
	my $tr1MM =  0;
	my $tr2MM =  0;
	my $tUnresolvedMM =  0; 
	my $tDCMM =  0; 
	my $tEAMM =  0;
	my $tSJAMM =  0;
	my $tSSJAMM =  0; 
	my $tSDMM = 0;	

	my $ntriplMMsd =  0;
	my $tsolvedMMsd =  0;
	my $tdifferentMMsd =  0;
	my $tr1MMsd =  0;
	my $tr2MMsd =  0;
	my $tUnresolvedMMsd =  0; 
	my $tDCMMsd =  0; 
	my $tEAMMsd =  0;
	my $tSJAMMsd =  0;
	my $tSSJAMMsd =  0; 
	my $tSDMMsd = 0;	

	
	my @Array_num_sp_elimMM = ();
	my @Array_e_n_taxa_tree1MM = ();
	my @Array_minMM = ();
	
	
  	my @Array_nquartMM = ();
	my @Array_qsolvedMM = ();
	my @Array_qdifferentMM = ();
	my @Array_qr2MM = ();
	my @Array_qr1MM = ();
	my @Array_qUnresolvedMM = ();
	my @Array_qDCMM = ();
	my @Array_qEAMM = ();
	my @Array_qSJAMM = ();
	my @Array_qSSJAMM = ();
	my @Array_qSDMM = ();

	
	my @Array_ntriplMM = ();
	my @Array_tsolvedMM = ();
	my @Array_tdifferentMM = ();
	my @Array_tr1MM = ();
	my @Array_tr2MM = ();
	my @Array_tUnresolvedMM = ();
	my @Array_tDCMM = ();
	my @Array_tEAMM = ();
	my @Array_tSJAMM = ();
	my @Array_tSSJAMM = ();
	my @Array_tSDMM = ();

	foreach $singTree1(@SingleGeneTrees1) {
		$contSingTrees2 = 0;
		foreach $singTree2(@SingleGeneTrees2) {
			$nomtreesingle1 = "$name1_INmm"."_$contSingTrees1";
			$nomtreesingle2 = "$name2_INmm"."_$contSingTrees2";
			
			($topd, $topdUnpruned, $com, $pcCom2, $SplitDist, $q, $m, $num_sp_elim, $e_n_taxa_tree1, $min, $DisTax, $nquart, $qsolved, $qdifferent, $qr2, $qr1, $qUnresolved, $qDC, $qEA, $qSJA, $qSSJA, $qSD, $ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &topd($nomtreesingle1, $singTree1, $nomtreesingle2, $singTree2, $print_prev, $method, $level, $units);
			
			
			$SplitDistMM += $SplitDist;
			$qMM += $q;
			$mMM += $m;

			push @SPLITS, $SplitDist;
			push @DIFFER, $q;
			push @TOTAL, $m;

			
			$topdMM += $topd;
			push @topds, $topd;
			$topdMMunpruned += $topdUnpruned;
			push @topdsUnpruned, $topdUnpruned;

			
			$num_sp_elimMM += $num_sp_elim;
			$e_n_taxa_tree1MM += $e_n_taxa_tree1;
			$minMM += $min;
			$DisTaxMM{$DisTax}++;
			
			
			$nquartMM += $nquart; 
			$qsolvedMM += $qsolved;
			$qdifferentMM += $qdifferent; 
			$qr2MM += $qr2; 
			$qr1MM += $qr1; 
			$qUnresolvedMM += $qUnresolved;
			$qDCMM += $qDC;
			$qEAMM += $qEA;
			$qSJAMM += $qSJA;
			$qSSJAMM += $qSSJA;
			$qSDMM += $qSD;
	
			
			$ntriplMM +=  $ntripl;
			$tsolvedMM +=  $tsolved;
			$tdifferentMM +=  $tdifferent;
			$tr1MM +=  $tr1;
			$tr2MM +=  $tr2;
			$tUnresolvedMM +=  $tUnresolved; 
			$tDCMM +=  $tDC; 
			$tEAMM +=  $tEA;
			$tSJAMM +=  $tSJA;
			$tSSJAMM +=  $tSSJA; 
			$tSDMM += $tSD;	
	
			
			push @Array_num_sp_elimMM, $num_sp_elim;
			push @Array_e_n_taxa_tree1MM, $e_n_taxa_tree1;
			push @Array_minMM, $min;
	
			
			push @Array_nquartMM, $nquart; 
			push @Array_qsolvedMM, $qsolved;
			push @Array_qdifferentMM, $qdifferent;
			push @Array_qr2MM, $qr2; 
			push @Array_qr1MM, $qr1; 
			push @Array_qUnresolvedMM, $qUnresolved;
			push @Array_qDCMM, $qDC;
			push @Array_qEAMM, $qEA;
			push @Array_qSJAMM, $qSJA;
			push @Array_qSSJAMM, $qSSJA;
			push @Array_qSDMM, $qSD;
	
			
			push @Array_ntriplMM, $ntripl;
			push @Array_tsolvedMM, $tsolved;
			push @Array_tdifferentMM, $tdifferent;
			push @Array_tr1MM, $tr1;
			push @Array_tr2MM, $tr2;
			push @Array_tUnresolvedMM, $tUnresolved; 
			push @Array_tDCMM, $tDC; 
			push @Array_tEAMM, $tEA;
			push @Array_tSJAMM, $tSJA;
			push @Array_tSSJAMM, $tSSJA; 
			push @Array_tSDMM, $tSD;

			$contSingTrees2++;
		}
		$contSingTrees1++;
	}
	
	my $SD = 0;
	my $SDunpruned = 0;
	if ($com <= 3) {
		($topdMM, $SD, $topdMMunpruned, $SDunpruned) = (0,0,0,0);
		if ($print_prev eq 'n') {
		}else{
			if ( ($method eq 'all') || ($method eq 'nodal') ) {
				printf "* Nodal Distance MM (Prunde/Unpruned): (%6.6f +/-%6.6f ) / (%6.6f +/-%6.6f )\n", $topdMM, $SD, $topdMMunpruned, $SDunpruned;
			}
			if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
				print "* Split Distance MM [differents/possibles]: 0 +/- 0.000 [ 0 +/- 0.000 \/ 0 +/- 0.000 ]\n";
			}
			if ( ($method eq 'all') || ($method eq 'disagree') ) {
				print "* Disagreement [ taxa disagree / all taxa ]: [0.000 +/- 0.000 / 0.000 +/- 0.000], New Split Distance: 0.000 +/- 0.000, Taxa disagree: [none]\n";
			}
			if ( ($method eq 'all') || ($method eq 'quartets') ) {
				print "* Quartets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
			}
			if ( ($method eq 'all') || ($method eq 'triplets') ) {
				print "* Triplets: 0.000 +/- 0.000, s: 0.000 +/- 0.000, d: 0.000  +/- 0.000, r1: 0.000  +/- 0.000, r2: 0.000 +/- 0.000, u: 0.000 +/- 0.000, DC: 0.000 +/- 0.000, EA: 0.000 +/- 0.000, SJA: 0.000 +/- 0.000, SSJA: 0.000 +/- 0.000, SD: 0.000 +/- 0.000\n";
			}
		}
	}else {
		
		$topdMM /= ( $contSingTrees1 * $contSingTrees2);
		$topdMMunpruned /= ( $contSingTrees1 * $contSingTrees2);

		
		foreach $topdsi(@topds) {
			$SD += ($topdsi - $topdMM)**2;
		}
		$SD = sqrt ( $SD / ( $contSingTrees1 * $contSingTrees2) );
	
		
		foreach $topdsuni(@topdsUnpruned) {
			$SDunpruned += ($topdsuni - $topdMMunpruned)**2;
		}
		$SDunpruned = sqrt ( $SDunpruned / ( $contSingTrees1 * $contSingTrees2) );

		
		$SplitDistMM /= ( $contSingTrees1 * $contSingTrees2);
		$qMM /= ( $contSingTrees1 * $contSingTrees2);
		$mMM /= ( $contSingTrees1 * $contSingTrees2);

		
		$SplitDistMMsd = 0;
		foreach $kk(@SPLITS) {
			$SplitDistMMsd += ($kk - $SplitDistMM)**2;
		}
		$SplitDistMMsd = sqrt ( $SplitDistMMsd / ( $contSingTrees1 * $contSingTrees2) );

		$qSMsd = 0;
		foreach $kk(@DIFFER) {
			$qMMsd += ($kk - $qMM)**2;
		}
		$qMMsd = sqrt ( $qMMsd / ( $contSingTrees1 * $contSingTrees2) );

		$mMMsd = 0;
		foreach $kk(@TOTAL) {
			$mMMsd += ($kk - $mMM)**2;
		}
		$mMMsd = sqrt ( $mMMsd / ( $contSingTrees1 * $contSingTrees2) );

		

		
		$num_sp_elimMM  /= ( $contSingTrees1 * $contSingTrees2);
		$e_n_taxa_tree1MM  /= ( $contSingTrees1 * $contSingTrees2);
		$minMM  /= ( $contSingTrees1 * $contSingTrees2);

		foreach $kk(@Array_num_sp_elimMM) {
			$num_sp_elimMMsd += ($kk - $num_sp_elimSM)**2;
		}
		$num_sp_elimMMsd = sqrt ($num_sp_elimMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_e_n_taxa_tree1MM) {
			$e_n_taxa_tree1MMsd += ($kk - $e_n_taxa_tree1MM)**2;
		}
		$e_n_taxa_tree1MMsd = sqrt ($e_n_taxa_tree1MMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_minMM) {
			$minMMsd += ($kk - $minMM)**2;
		}
		$minMMsd = sqrt ($minMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $DisTaxMMelement(sort keys(%DisTaxMM)) {
			$aux_dist__ = $DisTaxMM{$DisTaxMMelement};
			$DisTaxMMelement =~ s/-/ /g;
			$dif_elem .= "[$DisTaxMMelement("."$aux_dist__".")] ";
			
		}

		
  		$nquartMM  /= ( $contSingTrees1 * $contSingTrees2);
		$qsolvedMM  /= ( $contSingTrees1 * $contSingTrees2);
		$qdifferentMM  /= ( $contSingTrees1 * $contSingTrees2);
		$qr2MM  /= ( $contSingTrees1 * $contSingTrees2);
		$qr1MM  /= ( $contSingTrees1 * $contSingTrees2);
		$qUnresolvedMM  /= ( $contSingTrees1 * $contSingTrees2);
		$qDCMM  /= ( $contSingTrees1 * $contSingTrees2);
		$qEAMM  /= ( $contSingTrees1 * $contSingTrees2);
		$qSJAMM  /= ( $contSingTrees1 * $contSingTrees2);
		$qSSJAMM  /= ( $contSingTrees1 * $contSingTrees2);
		$qSDMM  /= ( $contSingTrees1 * $contSingTrees2);

		foreach $kk(@Array_nquartMM) {
			$nquartMMsd += ($kk - $nquartMM)**2;
		}
		$nquartMMsd = sqrt ($nquartMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_qsolvedMM) {
			$qsolvedMMsd += ($kk - $qsolvedSM)**2;
		}
		$qsolvedMMsd = sqrt ($qsolvedMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_qdifferentMM) {
			$qdifferentMMsd += ($kk - $qdifferentMM)**2;
		}
		$qdifferentMMsd = sqrt ($qdifferentMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_qr2MM) {
			$qr2MMsd += ($kk - $qr2MM)**2;
		}
		$qr2MMsd = sqrt ($qr2MMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_qr1MM) {
			$qr1MMsd += ($kk - $qr1MM)**2;
		}
		$qr1MMsd = sqrt ($qr1MMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_qUnresolvedMM) {
			$qUnresolvedMMsd += ($kk - $qUnresolvedMM)**2;
		}
		$qUnresolvedMMsd = sqrt ($qUnresolvedMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_qDCMM) {
			$qDCmMsd += ($kk - $qDCMM)**2;
		}	
		$qDCmMsd = sqrt ($qDCMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_qEAMM) {
			$qEAMMsd += ($kk - $qEASM)**2;
		}
		$qEAMMsd = sqrt ($qEAMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_qSJAMM) {
			$qSJAMMsd += ($kk - $qSJASM)**2;
		}
		$qSJAMMsd = sqrt ($qSJAMMsd/( $contSingTrees1 * $contSingTrees2));
		
		foreach $kk(@Array_qSSJAMM) {
			$qSSJAMMsd += ($kk - $qSSJAMM)**2;
		}
		$qSSJAMMsd = sqrt ($qSSJAMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_qSDMM) {
			$qSDMMsd += ($kk - $qSDMM)**2;
		}
		$qSDMMsd = sqrt ($qSDMMsd/( $contSingTrees1 * $contSingTrees2));

		
		$ntriplMM  /= ( $contSingTrees1 * $contSingTrees2);
		$tsolvedMM  /= ( $contSingTrees1 * $contSingTrees2);
		$tdifferentMM  /= ( $contSingTrees1 * $contSingTrees2);
		$tr1MM  /= ( $contSingTrees1 * $contSingTrees2);
		$tr2MM  /= ( $contSingTrees1 * $contSingTrees2);
		$tUnresolvedMM  /= ( $contSingTrees1 * $contSingTrees2);
		$tDCMM  /= ( $contSingTrees1 * $contSingTrees2);
		$tEAMM  /= ( $contSingTrees1 * $contSingTrees2);
		$tSJAMM  /= ( $contSingTrees1 * $contSingTrees2);
		$tSSJAMM  /= ( $contSingTrees1 * $contSingTrees2);
		$tSDMM  /= ( $contSingTrees1 * $contSingTrees2);

		foreach $kk(@Array_ntriplMM) {
			$ntriplMMsd += ($kk - $ntriplMM)**2;
		}
		$ntriplMMsd = sqrt ($ntriplMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_tsolvedMM) {
			$tsolvedMMsd += ($kk - $tsolvedMM)**2;
		}
		$tsolvedMMsd = sqrt ($tsolvedMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_tdifferentMM) {
			$tdifferentMMsd += ($kk - $tdifferentMM)**2;
		}
		$tdifferentMMsd = sqrt ($tdifferentMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_tr1MM) {
			$tr1MMsd += ($kk - $tr1MM)**2;
		}
		$tr1MMsd = sqrt ($tr1MMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_tr2MM) {
			$tr2MMsd += ($kk - $tr2MM)**2;
		}
		$tr2MMsd = sqrt ($tr2MMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_tUnresolvedMM) {
			$tUnresolvedMMsd += ($kk - $tUnresolvedMM)**2;
		}
		$tUnresolvedMMsd = sqrt ($tUnresolvedMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_tDCMM) {
			$tDCMMsd += ($kk - $tDCMM)**2;
		}
		$tDCMMsd = sqrt ($tDCMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_tEAMM) {
			$tEAMMsd += ($kk - $tEAMM)**2;
		}
		$tEAMMsd = sqrt ($tEAMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_tSJAMM) {
			$tSJAMMsd += ($kk - $tSJAMM)**2;
		}
		$tSJAMMsd = sqrt ($tSJAMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_tSSJAMM) {
			$tSSJAMMsd += ($kk - $tSSJAMM)**2;
		}
		$tSSJAMMsd = sqrt ($tSSJAMMsd/( $contSingTrees1 * $contSingTrees2));

		foreach $kk(@Array_tSDMM) {
			$tSDMMsd += ($kk - $tSDMM)**2;
		}
		$tSDMMsd = sqrt ($tSDMMsd/( $contSingTrees1 * $contSingTrees2));

		

		if ($print_prev eq 'n') {
		}else{
			if ( ($method eq 'all') || ($method eq 'nodal') ) {
				printf "* Nodal Distance MM (Prunde/Unpruned): (%6.6f +/- %6.6f ) / (%6.6f +/- %6.6f )\n", $topdMM, $SD, $topdMMunpruned, $SDunpruned;
			}
			if ( ($method eq 'all') || ($method eq 'split') || ($method eq 'disagree')) {
				printf "* Split Distance MM [differents/possibles]: $SplitDistMM +/- %5.3f [ $qMM +/- %5.3f \/ $mMM +/- %5.3f ]\n", $SplitDistMMsd, $qMMsd, $mMMsd ;
			}
			if ( ($method eq 'all') || ($method eq 'disagree') ) {
				printf "* Disagreement [ taxa disagree / all taxa ]: [%6.6f +/- %6.6f / %6.6f +/- %6.6f], New Split Distance: %6.6f +/- %6.6f, Taxa disagree: $dif_elem\n", $num_sp_elimMM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1MMsd, $minMM, $minMMsd ;
			}
			if ( ($method eq 'all') || ($method eq 'quartets') ) {
				printf "* Quartets: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $nquartMM, $nquartMMsd, $qsolvedMM, $qsolvedMMsd, $qdifferentMM, $qdifferentMMsd, $qr2MM, $qr2MMsd, $qr1MM, $qr1MMsd, $qUnresolvedMM, $qUnresolvedMMsd, $qDCMM, $qDCMMsd, $qEAMM, $qEAMMsd, $qSJAMM, $qSJAMMsd, $qSSJAMM, $qSSJAMMsd, $qSDMM, $qSDMMsd;
			}
			if ( ($method eq 'all') || ($method eq 'triplets') ) {
				printf "* Triplets: %6.6f +/- %6.6f, s: %6.6f +/- %6.6f, d: %6.6f  +/- %6.6f, r1: %6.6f  +/- %6.6f, r2: %6.6f +/- %6.6f, u: %6.6f +/- %6.6f, DC: %6.6f +/- %6.6f, EA: %6.6f +/- %6.6f, SJA: %6.6f +/- %6.6f, SSJA: %6.6f +/- %6.6f, SD: %6.6f +/- %6.6f\n", $ntriplMM, $ntriplMMsd, $tsolvedMM, $tsolvedMMsd, $tdifferentMM, $tdifferentMMsd,$tr2MM, $tr2MMsd, $tr1MM, $tr1MMsd, $tUnresolvedMM, $tUnresolvedMMsd, $tDCMM, $tDCMMsd, $tEAMM, $tEAMMsd, $tSJAMM, $tSJAMMsd, $tSSJAMM, $tSSJAMMsd, $tSDMM, $tSDMMsd;
			}
		}
	}
	
	
	return ($topdMM, $SD, $topdMMunpruned, $SDunpruned, $com, $pcCom2, $SplitDistMM, $SplitDistMMsd, $qMM, $qMMsd, $mMM, $mMMsd,$dif_elem, $num_sp_elimMM, $num_sp_elimMMsd, $e_n_taxa_tree1MM, $e_n_taxa_tree1MMsd, $minMM, $minMMsd, $nquartMM, $nquartMMsd, $qsolvedMM, $qsolvedMMsd, $qdifferentMM, $qdifferentMMsd, $qr2MM, $qr2MMsd, $qr1MM, $qr1MMsd, $qUnresolvedMM, $qUnresolvedMMsd, $qDCMM, $qDCMMsd, $qEAMM, $qEAMMsd, $qSJAMM, $qSJAMMsd, $qSSJAMM, $qSSJAMMsd, $qSDMM, $qSDMMsd, $ntriplMM, $ntriplMMsd, $tsolvedMM, $tsolvedMMsd, $tdifferentMM, $tdifferentMMsd,$tr2MM, $tr2MMsd, $tr1MM, $tr1MMsd, $tUnresolvedMM, $tUnresolvedMMsd, $tDCMM, $tDCMMsd, $tEAMM, $tEAMMsd, $tSJAMM, $tSJAMMsd, $tSSJAMM, $tSSJAMMsd, $tSDMM, $tSDMMsd);
}


sub fMtS() {
	my ($mgfNameTree1, $mgfTree1, $nombreSGTgenerats, $SingleMultiple, $print_prev) = @_;
	if ($print_prev eq 'n') {
	}else{
		print "\n\nFrom multiple to single gene family .... \n\n";
	}
	my @treesSingle = ();

	my ($i_v, $i_v, $i_v2, $i_v3, $contsp) = (0, 0,0,0,0);

	
	my $espMultiG = $mgfTree1;
	$espMultiG =~ s/\W/ /g;
	$espMultiG =~ s/\s+/-/g;

	my @espMultiGarray = split '-', $espMultiG;
	my %espMultiGhash = ();
	foreach $espmi(@espMultiGarray) {
		$espMultiGhash{$espmi}++;
	}	
	my %espMultiGhash_original = %espMultiGhash;

	
	%spambrepetit = ();
	my %numberreps = ();
	my $VARIANTS = 1;
	foreach $espmi2(@espMultiGarray) {
		if ($espMultiGhash{$espmi2} > 1) {
			$numberreps{$espmi2} = $espMultiGhash{$espmi2};;
			$VARIANTS *= $espMultiGhash{$espmi2};
			while ($espMultiGhash{$espmi2} > 0) {
				$mgfTree1 =~ s/(\W)($espmi2)(\W)/$1 $2 _ $espMultiGhash{$espmi2} $3/;
				$mgfTree1 =~ s/\s//g;	
				$espMultiGhash{$espmi2}--;
			}
			$spambrepetit{$espmi2} = 's';
		}
	}

	
	my $fileFMTS eq "";
	if ($mgfNameTree1 eq "") {
		$fileFMTS = "fmts";
		if ($SingleMultiple eq "single") {
			open (OUT, ">$fileFMTS");
		}
	}elsif ($mgfNameTree1 =~ /random/) {
	}
	else {
		$fileFMTS = "fmts"."_"."$mgfNameTree1";
		if ($SingleMultiple eq "single") {
			open (OUT, ">$fileFMTS");
		}
	}

	
	my %Species_rel = ();
	if ($nombreSGTgenerats eq "relative") {
		my $mgfTree1_rel = $mgfTree1;
		$mgfTree1_rel =~ s/_|\d+//g;

		my @sp_rel = split /\W+/, $mgfTree1_rel;
		my $NUM_TREES_rel = 1;
		foreach my $i_rel(@sp_rel) {
			if ($i_rel ne "") {
				$Species_rel{$i_rel}++;
			}
		}
		foreach my $i_rel(sort keys(%Species_rel)) {
			$NUM_TREES_rel *= $Species_rel{$i_rel};
		}

		if ($NUM_TREES_rel <= 100) {
			$nombreSGTgenerats = "all"; 
		}elsif ($NUM_TREES_rel > 100) {
			$nombreSGTgenerats = 100;
		}
	}

	if ($nombreSGTgenerats eq "all") {
		if ($print_prev eq 'n') {
		}else{
			print "\nThe number of single gene trees from $mgfNameTree1 is $VARIANTS.\n";
		}
		
		
		my %sp_variants = ();
		my %contsp = 0;
		my %numsp_var = 0;
		my %total_Variants = ();
		foreach $sp_rep(keys (%spambrepetit)) {
			
			my $cont_rp=0;
			$contsp++;
			while ($cont_rp <  $numberreps{$sp_rep}) {
				$cont_rp++;
				$afegir = "$sp_rep"."_$cont_rp";
				$sp_variants{$contsp} .= "$afegir"."#";
				$total_Variants{$afegir} = "v";
				$numsp_var{$contsp}++;
			}
			
		}
		
		
		$i_v=1;
		%numsp_var2 = %numsp_var;
		my $arbresGenerats = 0;
		while ($i_v <= $VARIANTS) {
			my @spapodarMM = ();
			my %spPodar = ();
			$i_v2=1;
			while ($i_v2 <= $contsp) {
				@VAriANTS = split /\#/, $sp_variants{$i_v2};
				
				$aux = @VAriANTS[$numsp_var2{$i_v2}-1];
				$spPodar{@VAriANTS[$numsp_var2{$i_v2}-1]} = "n";
				$i_v2++;
			}
			
			$i_v3=$contsp;
			while ($i_v3 >= 1) {
				if ($numsp_var2{$i_v3} == 1) {
					$numsp_var2{$i_v3} = $numsp_var{$i_v3};
				}else {
					$numsp_var2{$i_v3}--;
					last;
				}
				$i_v3--;
			}
			
			
			foreach $v(sort keys(%total_Variants)) {
				if ($spPodar{$v} eq "n") {
					
				}else {
					push @spapodarMM, $v;
				}
			}

			
			$podatTree = &pruneTree($mgfTree1, @spapodarMM);
			$podatTree = &unRootTrees($podatTree);
			

			
			$podatTree =~ s/_\d+//g;
			if ($SingleMultiple eq "single") {
				print OUT "$podatTree\n";
			}
			
			
			push @treesSingle, $podatTree;
			$arbresGenerats++;

			$i_v++;
		}
		
	}else {
		if ($print_prev eq 'n') {
		}else{
			print "\nThe fmts algorithm generates $nombreSGTgenerats single gene trees. The number of single gene trees from $mgfNameTree1 is $VARIANTS.\n";
		}
		my $arbresGenerats = 0;
		while ($arbresGenerats < $nombreSGTgenerats) {
			@spapodar=();
			$podatTree = ();
			
			
			my @spapodarMM = ();
			foreach $espmi2(keys (%spambrepetit)) {
				$aleat = int(rand($espMultiGhash_original{$espmi2})+1);
				$NOelim = "$espmi2"."_"."$aleat";
				$i__e=0;
				foreach $espmi3(@espMultiGarray) {
					if ($espmi3 eq $espmi2) {
						$i__e++;
						$elim = "$espmi3"."_"."$i__e";
						if ($elim eq $NOelim)  {
						}else {
							push @spapodarMM, $elim;
							
						}
					}
				}
			}
			
			
			
			$podatTree = &pruneTree($mgfTree1, @spapodarMM);
			$podatTree = &unRootTrees($podatTree);
			
			
			$podatTree =~ s/_\d+//g;
			if ($SingleMultiple eq "single") {
				print OUT "$podatTree\n";
			}
	
			
			push @treesSingle, $podatTree;
			$arbresGenerats++;
		}
	}

	if ($SingleMultiple eq "single") {
		close OUT;
	}
	
	return (@treesSingle);
}



sub fMtS2() {
	my ($mgfNameTree1, $mgfTree1, $nombreSGTgenerats, $SingleMultiple, $print_prev) = @_;
	if ($print_prev eq 'n') {
	}else{
		print "\n\nFrom multiple to single gene family .... \n\n";
	}
	my @treesSingle = ();

	my ($i_v, $i_v, $i_v2, $i_v3, $contsp) = (0, 0,0,0,0);

	
	my $espMultiG = $mgfTree1;
	$espMultiG =~ s/\W/ /g;
	$espMultiG =~ s/\s+/-/g;

	my @espMultiGarray = split '-', $espMultiG;
	my %espMultiGhash = ();
	foreach $espmi(@espMultiGarray) {
		$espMultiGhash{$espmi}++;
	}	
	my %espMultiGhash_original = %espMultiGhash;

	
	%spambrepetit = ();
	my %numberreps = ();
	my $VARIANTS = 1;
	foreach $espmi2(@espMultiGarray) {
		if ($espMultiGhash{$espmi2} > 1) {
			$numberreps{$espmi2} = $espMultiGhash{$espmi2};;
			$VARIANTS *= $espMultiGhash{$espmi2};
			while ($espMultiGhash{$espmi2} > 0) {
				$mgfTree1 =~ s/(\W)($espmi2)(\W)/$1 $2 _ $espMultiGhash{$espmi2} $3/;
				$mgfTree1 =~ s/\s//g;	
				$espMultiGhash{$espmi2}--;
			}
			$spambrepetit{$espmi2} = 's';
		}
	}

	
	my $fileFMTS eq "";
	if ($mgfNameTree1 eq "") {
		$fileFMTS = "fmts";
		if ($SingleMultiple eq "single") {
			open (OUT, ">$fileFMTS");
		}
	}elsif ($mgfNameTree1 =~ /random/) {
	}
	else {
		$fileFMTS = "fmts"."_"."$mgfNameTree1";
		if ($SingleMultiple eq "single") {
			open (OUT, ">$fileFMTS");
		}
	}

	
	my %Species_rel = ();
	if ($nombreSGTgenerats eq "relative") {
		my $mgfTree1_rel = $mgfTree1;
		$mgfTree1_rel =~ s/_|\d+//g;

		my @sp_rel = split /\W+/, $mgfTree1_rel;
		my $NUM_TREES_rel = 1;
		foreach my $i_rel(@sp_rel) {
			if ($i_rel ne "") {
				$Species_rel{$i_rel}++;
			}
		}
		foreach my $i_rel(sort keys(%Species_rel)) {
			$NUM_TREES_rel *= $Species_rel{$i_rel};
		}

		if ($NUM_TREES_rel <= 100) {
			$nombreSGTgenerats = "all"; 
		}elsif ($NUM_TREES_rel > 100) {
			$nombreSGTgenerats = 100;
		}
	}

	if ($nombreSGTgenerats eq "all") {
		if ($print_prev eq 'n') {
		}else{
			print "\nThe number of single gene trees from $mgfNameTree1 is $VARIANTS.\n";
		}
		
		
		my %sp_variants = ();
		my %contsp = 0;
		my %numsp_var = 0;
		my %total_Variants = ();
		foreach $sp_rep(keys (%spambrepetit)) {
			
			my $cont_rp=0;
			$contsp++;
			while ($cont_rp <  $numberreps{$sp_rep}) {
				$cont_rp++;
				$afegir = "$sp_rep"."_$cont_rp";
				$sp_variants{$contsp} .= "$afegir"."#";
				$total_Variants{$afegir} = "v";
				$numsp_var{$contsp}++;
			}
			
		}
		
		
		$i_v=1;
		%numsp_var2 = %numsp_var;
		my $arbresGenerats = 0;
		while ($i_v <= $VARIANTS) {
			my @spapodarMM = ();
			my %spPodar = ();
			$i_v2=1;
			while ($i_v2 <= $contsp) {
				@VAriANTS = split /\#/, $sp_variants{$i_v2};
				
				$aux = @VAriANTS[$numsp_var2{$i_v2}-1];
				$spPodar{@VAriANTS[$numsp_var2{$i_v2}-1]} = "n";
				$i_v2++;
			}
			
			$i_v3=$contsp;
			while ($i_v3 >= 1) {
				if ($numsp_var2{$i_v3} == 1) {
					$numsp_var2{$i_v3} = $numsp_var{$i_v3};
				}else {
					$numsp_var2{$i_v3}--;
					last;
				}
				$i_v3--;
			}
			
			
			foreach $v(sort keys(%total_Variants)) {
				if ($spPodar{$v} eq "n") {
					
				}else {
					push @spapodarMM, $v;
				}
			}

			
			$podatTree = &pruneTree($mgfTree1, @spapodarMM);
			
			

			
			$podatTree =~ s/_\d+//g;
			if ($SingleMultiple eq "single") {
				print OUT "$podatTree\n";
			}
			
			
			push @treesSingle, $podatTree;
			$arbresGenerats++;

			$i_v++;
		}
		
	}else {
		if ($print_prev eq 'n') {
		}else{
			print "\nThe fmts algorithm generates $nombreSGTgenerats single gene trees. The number of single gene trees from $mgfNameTree1 is $VARIANTS.\n";
		}
		my $arbresGenerats = 0;
		while ($arbresGenerats < $nombreSGTgenerats) {
			@spapodar=();
			$podatTree = ();
			
			
			my @spapodarMM = ();
			foreach $espmi2(keys (%spambrepetit)) {
				$aleat = int(rand($espMultiGhash_original{$espmi2})+1);
				$NOelim = "$espmi2"."_"."$aleat";
				$i__e=0;
				foreach $espmi3(@espMultiGarray) {
					if ($espmi3 eq $espmi2) {
						$i__e++;
						$elim = "$espmi3"."_"."$i__e";
						if ($elim eq $NOelim)  {
						}else {
							push @spapodarMM, $elim;
							
						}
					}
				}
			}
			
			
			
			$podatTree = &pruneTree($mgfTree1, @spapodarMM);
			
			
			
			$podatTree =~ s/_\d+//g;
			if ($SingleMultiple eq "single") {
				print OUT "$podatTree\n";
			}
	
			
			push @treesSingle, $podatTree;
			$arbresGenerats++;
		}
	}

	if ($SingleMultiple eq "single") {
		close OUT;
	}
	
	return (@treesSingle);
}



sub random() {
	my ($name1_IN, $tree1_IN, $name2_IN, $tree2_IN) = @_;
	my ($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom) = ();

	$name1_INrandom = "$name1_IN"."_random";
	$tree1_INrandom = $tree1_IN;
	$name2_INrandom = "$name2_IN"."_random";
	$tree2_INrandom = $tree2_IN;

	
	$tree1_INrandom =~ s/\(/ /g;
	$tree1_INrandom =~ s/\)/ /g;
	$tree1_INrandom =~ s/,/ /g;
	$tree1_INrandom =~ s/;//g;
	$tree1_INrandom =~ s/\s+/\t/g;
	my @sprandom1 = split '\t', $tree1_INrandom;
	
	
	my $sp_rand1 = 0;
	foreach my $i_rand(@sprandom1) {
		if ($i_rand ne '') {
			$sp_rand1++;
		}
	}
	
	
	$tree2_INrandom =~ s/\(/ /g;
	$tree2_INrandom =~ s/\)/ /g;
	$tree2_INrandom =~ s/,/ /g;
	$tree2_INrandom =~ s/;//g;
	$tree2_INrandom =~ s/\s+/\t/g;

	my @sprandom2 = split '\t', $tree2_INrandom;

	
	my $sp_rand2 = 0;
	foreach my $i_rand(@sprandom2) {
		if ($i_rand ne '') {
			$sp_rand2++;
		}
	}

	
	$tree1_INrandom = $tree1_IN;
	$tree2_INrandom = $tree2_IN;
	$tree1_INrandom =~ s/\w+/\*/g;
	$tree2_INrandom =~ s/\w+/\*/g;

	
	my @sprandom1_col = @sprandom1;
	while ($tree1_INrandom =~ /\*/) {
		$random_rand = int(rand($sp_rand1)+1);
		if (@sprandom1_col[$random_rand] ne "") {
			$tree1_INrandom =~ s/\*/@sprandom1_col[$random_rand]/;
			@sprandom1_col[$random_rand] = "";
		}
	}

	
	my @sprandom2_col = @sprandom2;
	while ($tree2_INrandom =~ /\*/) {
		$random_rand = int(rand($sp_rand2)+1);
		if (@sprandom2_col[$random_rand] ne "") {
			$tree2_INrandom =~ s/\*/@sprandom2_col[$random_rand]/;
			@sprandom2_col[$random_rand] = "";
		}
	}
	return ($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom);
}


sub random2() {
	my ($name1_IN, $tree1_IN, $name2_IN, $tree2_IN) = @_;
	$name1_INrandom = "$name1_IN"."_random";
	$tree1_INrandom = $tree1_IN;
	$name2_INrandom = "$name2_IN"."_random";
	$tree2_INrandom = $tree2_IN;

	
	$tree1_INrandom =~ s/\(/ /g;
	$tree1_INrandom =~ s/\)/ /g;
	$tree1_INrandom =~ s/,/ /g;
	$tree1_INrandom =~ s/;//g;
	$tree1_INrandom =~ s/\s+/\t/g;
	my @sprandom1 = split '\t', $tree1_INrandom;
	
	
	my $sp_rand1 = 0;
	foreach my $i_rand(@sprandom1) {
		if ($i_rand ne '') {
			$sp_rand1++;
		}
	}
	
	
	$tree2_INrandom =~ s/\(/ /g;
	$tree2_INrandom =~ s/\)/ /g;
	$tree2_INrandom =~ s/,/ /g;
	$tree2_INrandom =~ s/;//g;
	$tree2_INrandom =~ s/\s+/\t/g;

	my @sprandom2 = split '\t', $tree2_INrandom;

	
	my $sp_rand2 = 0;
	foreach my $i_rand(@sprandom2) {
		if ($i_rand ne '') {
			$sp_rand2++;
		}
	}
	
	
	$tree1_INrandom = &findRandomTopology($sp_rand1);
	
	

	
	$tree2_INrandom = &findRandomTopology($sp_rand2);
	
	

	
	$tree1_INrandom =~ s/\w+/\*/g;
	$tree2_INrandom =~ s/\w+/\*/g;
	
	

	
	my @sprandom1_col = @sprandom1;
	while ($tree1_INrandom =~ /\*/) {
		$random_rand = int(rand($sp_rand1)+1);
		if (@sprandom1_col[$random_rand] ne "") {
			$tree1_INrandom =~ s/\*/@sprandom1_col[$random_rand]/;
			@sprandom1_col[$random_rand] = "";
		}
	}

	
	my @sprandom2_col = @sprandom2;
	while ($tree2_INrandom =~ /\*/) {
		$random_rand = int(rand($sp_rand2)+1);
		if (@sprandom2_col[$random_rand] ne "") {
			$tree2_INrandom =~ s/\*/@sprandom2_col[$random_rand]/;
			@sprandom2_col[$random_rand] = "";
		}
	}
	
	
	
	return ($name1_INrandom, $tree1_INrandom, $name2_INrandom, $tree2_INrandom);
}
sub findRandomTopology {
	my $numSP = shift;
	my $randomTopology = "";

	
	my %sprtop = ();
	for (my $irtop=0; $irtop < $numSP; $irtop++) {
		$sprtop{$irtop} = $irtop;
		
	}

	
	my $spdelhash = $numSP;
	my $snrt = 0;
	while ($snrt < $numSP) {
		if ($sprtop{'3'} eq '') {
			$randomTopology = "($sprtop{'0'},$sprtop{'1'},$sprtop{'2'})";
		}else {
			
			my $ale1 = "";
			my $ale2 = "";
			while ($ale1 == $ale2) {
				$ale1 = int(rand($spdelhash));
				$ale2 = int(rand($spdelhash));
			}
	
			

			
			
			
			
	
			$randomTopology = "($sprtop{$ale1},$sprtop{$ale2})";

			
			$sprtop{$ale1} = "";
			$sprtop{$ale2} = "";
		
			%sprtopAUX = ();
			$spdelhash=0;
			foreach $sp_hash(keys(%sprtop)) {
				if ($sprtop{$sp_hash} eq '') {
				}else {
					$sprtopAUX{$spdelhash} = $sprtop{$sp_hash};
					
					$spdelhash++;
				}
			}
			$sprtopAUX{$spdelhash} = $randomTopology;
			
			
			%sprtop = %sprtopAUX;
		}
		
		
		$snrt = $randomTopology =~ s/(\d+)/$1/g;
		
		
	}
	$randomTopology .= ";";
	return $randomTopology;
}

sub meanSD() {
	my @data = @_;
	my ($sub_mean,$sub_sd) = ();
	my $cont_data =0;
	my $sum_data = 0;
	my $sum_dif = 0; 
	
	foreach my $datai(@data) {
		$sum_data += $datai;
		$cont_data++;
	}
	$sub_mean = $sum_data / $cont_data;

	
	foreach my $dataisd(@data) {
		$sum_dif += ($dataisd - $sub_mean)**2 ;
	}
	$sub_sd = sqrt ( $sum_dif / $cont_data );
	return ($sub_mean,$sub_sd);
}






sub SP {
	my $LocTrees=shift;
	my %name_tree = &getTrees($LocTrees);
	
	open (OUTsp, ">species");
	my $name = "";
	my $tree = "";
	foreach $name(sort keys(%name_tree)) {
		my @species = ();
		my %Species = ();
		$tree = $name_tree{$name};
	
		@sp = split /\W+/, $tree;
		foreach my $i(@sp) {
			if ($i ne "") {
				$Species{$i}++;
			}
		}

		print OUTsp ">$name:$tree\n";
		
	
		foreach my $i(sort keys(%Species)) {
			print OUTsp "$i ($Species{$i})\t";
		}
		print OUTsp "\n";
	}
	close OUTsp;
	
	print "\nThanks for using TOPD-fMtS .................................. Bye.\n\n";
}






sub prunne_trees {
	my $LocTrees=shift;
	my $LocSpecies=shift;
	my %name_tree = &getTrees($LocTrees);
	
	open (INspecies, "<$LocSpecies") || die "cannot open $LocSpecies\n";
	my @speciestoprunne = ();
	while (<INspecies>) {
		chomp;
		if ( ($_ ne "") && ($_ ne '#') ) {
			push @speciestoprunne, $_;
		}
	}	
	close INspecies;

	open (OUTprunned, ">pruned_trees");
	my $name = "";
	my $tree = "";
	foreach $name(sort keys(%name_tree)) {
		$tree = $name_tree{$name};
		
		
		$pTree = &pruneTree($tree, @speciestoprunne);

		
		$pTree = &unRootTrees($pTree);

		print OUTprunned "$name\t$pTree\n";
	}
	close OUTprunned;

	print "\nThanks for using TOPD-fMtS .................................. Bye.\n\n";
}





sub fmts_trees {
	my $LocTrees = shift;
	my $type = shift;
	my %name_tree = &getTrees($LocTrees);

	my $name = "";
	open (OUTfmts, ">fmts");
	foreach $name(sort keys(%name_tree)) {
		$tree = $name_tree{$name};
		my @species = ();
		my %Species = ();
		$tree = $name_tree{$name};

		if ($type eq "relative") {
			@sp = split /\W+/, $tree;
			$NUM_TREES = 1;
			foreach my $i(@sp) {
				if ($i ne "") {
					$Species{$i}++;
				}
			}
			foreach my $i(sort keys(%Species)) {
				$NUM_TREES *= $Species{$i};
			}
	
			if ($NUM_TREES < 100) {
				my @SingleGeneTrees1 = &fMtS('random', $tree, 'all', 'multiple' , 'n');
				print OUTfmts ">all_$NUM_TREES|$name:$tree\n";
				foreach my $single (@SingleGeneTrees1) {
					print OUTfmts "$single\n";
				}
			}elsif ($NUM_TREES > 100) {
				my @SingleGeneTrees1 = &fMtS('random', $tree, '100', 'multiple' , 'n');
				print OUTfmts ">random_100|$name:$tree\n";
				foreach my $single (@SingleGeneTrees1) {
					print OUTfmts "$single\n";
				}
			}	
		}elsif ($type eq "random") {
			my @SingleGeneTrees1 = &fMtS('random', $tree, '100', 'multiple' , 'n');
			print OUTfmts ">$name:$tree\n";
			foreach my $single (@SingleGeneTrees1) {
				print OUTfmts "$single\n";
			}
		}elsif ($type eq "all")  {
			my @SingleGeneTrees1 = &fMtS('random', $tree, 'all', 'multiple' , 'n');
			print OUTfmts ">$name:$tree\n";
			foreach my $single (@SingleGeneTrees1) {
				print OUTfmts "$single\n";
			}
		}else {
			my @SingleGeneTrees1 = &fMtS('random', $tree, 'all', 'multiple' , 'n');
			print OUTfmts ">$name:$tree\n";
			foreach my $single (@SingleGeneTrees1) {
				print OUTfmts "$single\n";
			}
		}
	}

	print "\nThanks for using TOPD-fMtS .................................. Bye.\n\n";
	
	close OUTfmts;
}







sub getQuartets() {
	my $LocTrees = shift;
	my %name_tree = &getTrees($LocTrees);	

	open (OUTquart, ">quartets");
	foreach $nametrees(sort keys(%name_tree)) {
		$treeQuartet1 = $name_tree{$nametrees};

		my @species = split /\W+/, $treeQuartet1;
		my %SPECIES = ();
		foreach my $aux(@species) {
			$SPECIES{$aux} = "";
		}
		
		foreach $a(keys(%SPECIES)) {
			while ($treeQuartet1 =~ s/\($a,$a\)/$a/) {
				
			}
		}

		print OUTquart ">$nametrees|$treeQuartet1|";
		
		

		
		my @q_taxa_tree1 = sort split '\W+', $treeQuartet1;
		my $q_n_taxa_tree1 = (@q_taxa_tree1 - 1);
	
		my $number_quartets = $q_n_taxa_tree1*($q_n_taxa_tree1-1)*($q_n_taxa_tree1-2)*($q_n_taxa_tree1-3) / 24;
		print OUTquart "$number_quartets\n";

		
		my $qi = 0;
		my $nquart = 0;
		my %quartetsFets = ();
		foreach $qsp1(sort(@q_taxa_tree1)) {
			foreach $qsp2(sort(@q_taxa_tree1)) {
				foreach $qsp3(sort(@q_taxa_tree1)) {
					foreach $qsp4(sort(@q_taxa_tree1)) {
						if ( ($qsp1 ne $qsp2) && ($qsp1 ne $qsp3) && ($qsp1 ne $qsp4) && ($qsp2 ne $qsp3) && ($qsp2 ne $qsp4) && ($qsp3 ne $qsp4) && ($qsp1 ne "") && ($qsp2 ne "") && ($qsp3 ne "") && ($qsp4 ne "")  ) {
							my @AuxQuartets = ("$qsp1", "$qsp2", "$qsp3", "$qsp4");
							@AuxQuartets = sort @AuxQuartets;
	
							$qspa1 = @AuxQuartets[0];
							$qspa2 = @AuxQuartets[1];
	
							$qspa3 = @AuxQuartets[2];
							$qspa4 = @AuxQuartets[3];
	
							if ($quartetsFets{"$qspa1-$qspa2-$qspa3-$qspa4"} ne "Y") {
								print OUTquart "$qspa1-$qspa2-$qspa3-$qspa4\t";
								my @qelimsp = ();
								foreach $qe(@q_taxa_tree1) {
									if ( "$qspa1-$qspa2-$qspa3-$qspa4" !~ /$qe/) {
										push @qelimsp, $qe;
									}
								}
								
								my $qpTree1 = &pruneTree($treeQuartet1, @qelimsp);
								
								
								
								my @quartetselim = &remove_quartets($qpTree1);
								if ($quartetselim[0] ne "") {
									$qpTree1 = $quartetselim[0];
									if ($quartetselim[1] ne "") {
										$qpTree1 .= "\t$quartetselim[1]";
										if ($quartetselim[2] ne "") {
											$qpTree1 .= "\t$quartetselim[2]";
											if ($quartetselim[3] ne "") {

											}
										}
									}
								}

								
								
								print OUTquart "$qpTree1\n";
							
								$quartetsFets{"$qspa1-$qspa2-$qspa3-$qspa4"} = "Y";
								$nquart++;
							}
						}
					}
				}
			}	
		}
	}

	print "\nThanks for using TOPD-fMtS .................................. Bye.\n\n";
	
	close OUTquart;
}	



sub remove_quartets {
	my $treerem = shift;

	
	
	my @SingleGeneTreesrem = &fMtS('random', $treerem, 'relative', 'multiple' , 'n');

	
	my %AUX = ();
	foreach my $aux(@SingleGeneTreesrem) {
		$AUX{$aux} = "";
	}
	

	my %quartetsrem = ();
	my %fetsrem = ();
	foreach my $tree1rem(sort keys %AUX) {
		foreach my $tree2rem(sort keys %AUX) {
			
			if ( ($tree1rem ne $tree2rem) && ($fetsrem{$tree1rem}{$tree2rem} eq "") && ($fetsrem{$tree2rem}{$tree1rem}) eq "") {
				my ($spdist, $q, $m) = &splitDistance($tree1rem,$tree2rem);
				if ( $spdist == 0) {
				
				
					if ( ($quartetsrem{$tree1rem} eq "y") || ($quartetsrem{$tree1rem} eq "") ) {
						$quartetsrem{$tree1rem} = "y";
						$quartetsrem{$tree2rem} = "n";
					}elsif ( ($quartetsrem{$tree2rem} eq "y") || ($quartetsrem{$tree2rem} eq "") ) {
						$quartetsrem{$tree1rem} = "n";
						$quartetsrem{$tree2rem} = "y";
					}
				}else {
					if ($quartetsrem{$tree1rem} ne "n")  {
						$quartetsrem{$tree1rem} = "y";
					}
					if ($quartetsrem{$tree2rem} ne "n")  {
						$quartetsrem{$tree2rem} = "y";
					}
				}
				
				$fetsrem{$tree1rem}{$tree2rem} = "y";
			}else {
				if ($quartetsrem{$tree1rem} ne "n")  {
					$quartetsrem{$tree1rem} = "y";
				}
			}
		}
	}

	my @QUARTETSREM = ();
	foreach $aaa(keys (%quartetsrem)) {
		if ($quartetsrem{$aaa} eq "y") {
			push @QUARTETSREM, $aaa;
			
		}
	}
	
	
	return @QUARTETSREM;
}





sub getTriplets() {

	my $LocTrees = shift;
	my %name_tree = &getTrees($LocTrees);
	
	open (OUTtriplets, ">triplets");
	foreach $nametrees(sort keys(%name_tree)) {
		$treeTriplet1 = $name_tree{$nametrees};



		my @species = split /\W+/, $treeTriplet1;
		my %SPECIES = ();
		foreach my $aux(@species) {
			$SPECIES{$aux} = "";
		}
		
		foreach $a(keys(%SPECIES)) {
			while ($treeTriplet1 =~ s/\($a,$a\)/$a/) {
				
			}
		}

		print OUTtriplets ">$nametrees|$treeTriplet1|";

		
		my @t_taxa_tree1 = sort split '\W+', $treeTriplet1;
		my $t_n_taxa_tree1 = (@t_taxa_tree1 - 1);


		my $number_triplets = $t_n_taxa_tree1*($t_n_taxa_tree1-1)*($t_n_taxa_tree1-2) / 6;
		print OUTtriplets "$number_triplets\n";

		
		my $ti = 0;
		my $ntripl = 0;
		my %tripletsFets = ();
		foreach $tsp1(sort(@t_taxa_tree1)) {
			foreach $tsp2(sort(@t_taxa_tree1)) {
				foreach $tsp3(sort(@t_taxa_tree1)) {
					
					if ( ($tsp1 ne $tsp2) && ($tsp1 ne $tsp3) && ($tsp2 ne $tsp3) && ($tsp1 ne "") && ($tsp2 ne "") && ($tsp3 ne "") ) {
						my @AuxTriplets = ("$tsp1", "$tsp2", "$tsp3");
						@AuxTriplets = sort @AuxTriplets;
						$tspa1 = @AuxTriplets[0];
						$tspa2 = @AuxTriplets[1];
						$tspa3 = @AuxTriplets[2];
						if ($tripletsFets{"$tspa1-$tspa2-$tspa3"} ne "Y") {
							print OUTtriplets "$tspa1-$tspa2-$tspa3\t";
							
							my @telimsp = ();
							foreach $te(@t_taxa_tree1) {
								if ( "$tspa1-$tspa2-$tspa3" !~ /$te/) {
									push @telimsp, $te;
								}
							}
							
							my $tpTree1 = &pruneTree($treeTriplet1, @telimsp);
	
							
							
							my @tripletselim = &remove_triplets($tpTree1);
							
							if ($tripletselim[0] ne "") {
								$tpTree1 = $tripletselim[0];
								if ($tripletselim[1] ne "") {
									$tpTree1 .= "\t$tripletselim[1]";
									if ($tripletselim[2] ne "") {
										$tpTree1 .= "\t$tripletselim[2]";
										if ($tripletselim[3] ne "") {
											$tpTree1 .= "\t$tripletselim[3]";
										}
									}
								}
							}
							


							print OUTtriplets "$tpTree1\n";
							
							$tripletsFets{"$tspa1-$tspa2-$tspa3"} = "Y";
							$ntripl++;
						}
					}
				}
			}
		}
	}

	print "\nThanks for using TOPD-fMtS .................................. Bye.\n\n";
}

sub remove_triplets {
	my $treerem = shift;

	
	
	my @SingleGeneTreesrem = &fMtS2('random', $treerem, 'relative', 'multiple' , 'n');

	
	my %AUX = ();
	foreach my $aux(@SingleGeneTreesrem) {
		
		$AUX{$aux} = "";
	}
	

	my %tripletsrem = ();
	my %fetsrem = ();
	foreach my $tree1rem(sort keys %AUX) {
		foreach my $tree2rem(sort keys %AUX) {
			
			if ( ($tree1rem ne $tree2rem) && ($fetsrem{$tree1rem}{$tree2rem} eq "") && ($fetsrem{$tree2rem}{$tree1rem}) eq "") {
				
				my ($ntripl, $tsolved, $tdifferent, $tr1, $tr2, $tUnresolved, $tDC, $tEA, $tSJA, $tSSJA, $tSD) = &tripletsDistance($tree1rem, $tree2rem, 'all');
				if ( $tsolved == 1) {
					if ( ($tripletsrem{$tree1rem} eq "y") || ($tripletsrem{$tree1rem} eq "") ) {
						$tripletsrem{$tree1rem} = "y";
						$tripletsrem{$tree2rem} = "n";
					}elsif ( ($tripletsrem{$tree2rem} eq "y") || ($tripletsrem{$tree2rem} eq "") ) {
						$tripletsrem{$tree1rem} = "n";
						$tripletsrem{$tree2rem} = "y";
					}
				}else {
					if ($tripletsrem{$tree1rem} ne "n")  {
						$tripletsrem{$tree1rem} = "y";
					}
					if ($tripletsrem{$tree2rem} ne "n")  {
						$tripletsrem{$tree2rem} = "y";
					}
				}
				
				$fetsrem{$tree1rem}{$tree2rem} = "y";
			}else {
				if ($tripletsrem{$tree1rem} ne "n")  {
					$tripletsrem{$tree1rem} = "y";
				}
			}
		}
	}

	my @TRIPLETSREM= ();
	foreach $aaa(keys (%tripletsrem)) {
		if ($tripletsrem{$aaa} eq "y") {
			push @TRIPLETSREM, $aaa;
			
		}
	}
	
	
	return @TRIPLETSREM;
}






sub quickquartets {
	my $LocTrees = shift;
	my $numberofquartets = shift;

	
	my %name_tree = &getTrees($LocTrees);

	if ( ($numberofquartets eq "") || ($numberofquartets =~ /\D/) ) {
		$numberofquartets = 100;
	}

	open (OUTqq, ">qq_result");


	my %fetsqq = ();
	foreach $nametrees1(sort keys(%name_tree)) {
		foreach $nametrees2(sort keys(%name_tree)) {
			if ( ($nametrees1 ne $nametrees2) && ($fetsqq{$nametrees1}{$nametrees2} ne "y") && ($fetsqq{$nametrees2}{$nametrees1} ne "y") ) {
				$fetsqq{$nametrees1}{$nametrees2} = "y";
				
				
				$treeQuartet1 = $name_tree{$nametrees1};
				$treeQuartet2 = $name_tree{$nametrees2};

				my %SPECIEStot = ();
				
				my @species1 = split /\W+/, $treeQuartet1;
				my %SPECIES1 = ();
				foreach my $aux(@species1) {
					$SPECIEStot{$aux}++;
					$SPECIES1{$aux} = "y";
				}
				foreach my $a(keys(%SPECIES1)) {
					while ($treeQuartet1 =~ s/\($a,$a\)/$a/) {
						
					}
				}

				
				my @species2 = split /\W+/, $treeQuartet2;
				my %SPECIES2 = ();
				foreach my $aux(@species2) {
					$SPECIEStot{$aux}++;
					$SPECIES2{$aux} = "y";
				}
				foreach my $a(keys(%SPECIES2)) {
					while ($treeQuartet2 =~ s/\($a,$a\)/$a/) {
						
					}
				}
				
				
				my @nocomuns1 = ();
				my @nocomuns2 = ();
				my @comuns = ();
				foreach $aux(keys %SPECIEStot) {
					if ($aux ne "") {
						if ( ($SPECIES1{$aux} ne "y") && ($SPECIES2{$aux} ne "n") ) {
							for (my $qqq = 0;$qqq < $SPECIEStot{$aux};$qqq++) {
								push @nocomuns1, $aux;
							}
						}
						elsif ( ($SPECIES1{$aux} ne "n") && ($SPECIES2{$aux} ne "y") ) {
							for (my $qqq = 0;$qqq < $SPECIEStot{$aux};$qqq++) {
								push @nocomuns2, $aux;
							}
						}elsif( ($SPECIES1{$aux} eq "y") && ($SPECIES2{$aux} eq "y") ) {
							push @comuns, $aux;
						}
					}
				}

				my $spscomuns = @comuns;
				if ($spscomuns >= 4) {
					
					$treeQuartet1 = &pruneTree($treeQuartet1, @nocomuns1);
					$treeQuartet2 = &pruneTree($treeQuartet2, @nocomuns2);
					
					
					
	
					
	
					
					my @q_taxa_tree1 = sort split '\W+', $treeQuartet1;
					my $q_n_taxa_tree1 = (@q_taxa_tree1 - 1);
	
					
					my @q_taxa_tree2 = sort split '\W+', $treeQuartet2;
					my $q_n_taxa_tree2 = (@q_taxa_tree2 - 1);
	
					my $number_quartets = $q_n_taxa_tree1*($q_n_taxa_tree1-1)*($q_n_taxa_tree1-2)*($q_n_taxa_tree1-3) / 24;
	
					
					my $qi = 0;
					my $qqdist = "";
					my $ndist = "";
	
					my $qqdist2 = "";
					my $ndist2 = "";
	
					my $totalsp = @comuns;
					for (my $nquart = 0;$nquart < $numberofquartets;$nquart++) {
						
	
						
						my $a = "";
						my $b = "";
						my $c = "";
						my $d = "";
	
						$a = int(rand($totalsp));
						$qsp1 = $comuns[$a];
	
						while (($b eq "") || ($qsp2 eq $qsp1)) {
							$b = int(rand($totalsp));
							$qsp2 = $comuns[$b];
						}
						while ( ($c eq "") || ($qsp3 eq $qsp1) || ($qsp3 eq $qsp2) ) {
							$c = int(rand($totalsp));
							$qsp3 = $comuns[$c];
						}
						while ( ($d eq "") || ($qsp4 eq $qsp1) || ($qsp4 eq $qsp2) || ($qsp4 eq $qsp3) ) {
							$d = int(rand($totalsp));
							$qsp4 = $comuns[$d];
						}
						
	
						
						
						my @qelimsp = ();
						foreach $qe(@q_taxa_tree1) {
							if ( "$qsp1-$qsp2-$qsp3-$qsp4" !~ /$qe/) {
								push @qelimsp, $qe;
							}
						}
						
						my $qpTree1 = "";
						$qpTree1 = &pruneTree($treeQuartet1, @qelimsp);
						
						
						my @quartetselim = &remove_quartets($qpTree1);
						my $nqqq1 = 0;
						if ($quartetselim[0] ne "") {
							$nqqq1++;
							$qpTree1 = $quartetselim[0];
							if ($quartetselim[1] ne "") {
								$qpTree1 .= "-$quartetselim[1]";
								$nqqq1++;
							}
							if ($quartetselim[2] ne "") {
								$qpTree1 .= "-$quartetselim[2]";
								$nqqq1++;
							}
						}
						
						
	
						
						my @qelimsp = ();
						foreach $qe(@q_taxa_tree2) {
							if ( "$qsp1-$qsp2-$qsp3-$qsp4" !~ /$qe/) {
								push @qelimsp, $qe;
							}
						}
						
						my $qpTree2 = "";
						$qpTree2 = &pruneTree($treeQuartet2, @qelimsp);
						
						
						my @quartetselim = &remove_quartets($qpTree2);
						my $nqqq2 = 0;
						if ($quartetselim[0] ne "") {
							$nqqq2++;
							$qpTree2 = $quartetselim[0];
							if ($quartetselim[1] ne "") {
								$qpTree2 .= "-$quartetselim[1]";
								$nqqq2++;
							}
							if ($quartetselim[2] ne "") {
								$qpTree2 .= "-$quartetselim[2]";
								$nqqq2++;
							}
						}
						
						
						
						
						
						
						
						if ( ($nqqq1 == 3) && ($nqqq2 == 3) ) {
							$qqdist += 0;
						}elsif ( ( ($nqqq1 == 3) && ($nqqq2 == 2) ) || ( ($nqqq1 == 2) && ($nqqq2 == 3) ) ) {
							$qqdist += 1/3;
						}
						elsif ( ( ($nqqq1 == 3) && ($nqqq2 == 1) ) || ( ($nqqq1 == 1) && ($nqqq2 == 3) ) ) {
							$qqdist += 2/3;
						}
						elsif ( ($nqqq1 == 2) && ($nqqq2 ==2)) {
							$qpTree1=~ /(.+)\-(.+)/;
							my $qa1 = $1;
							my $qa2 = $2;
							$qpTree2 =~ /(.+)\-(.+)/;
							my $qb1 = $1;
							my $qb2 = $2;
						
							my ($spdist1, $q1, $m1) = &splitDistance($qa1,$qb1);
							my ($spdist2, $q2, $m2) = &splitDistance($qa1,$qb2);
							my ($spdist3, $q1, $m1) = &splitDistance($qa2,$qb1);
							my ($spdist4, $q2, $m2) = &splitDistance($qa2,$qb2);
						
							$dot = $spdist1 + $spdist2 + $spdist3 + $spdist4;
		
							if ($dot == 2) {
								$qqdist += 0;
							}elsif ($dot == 3) {
								$qqdist += 2/3;
							}
						}
						elsif ( ($nqqq1 == 1) && ($nqqq2 ==2)) {
							$qpTree2 =~ /(.+)\-(.+)/;
							my $qa1 = $1;
							my $qa2 = $2;
							my ($spdist1, $q1, $m1) = &splitDistance($qpTree1,$qa1);
							my ($spdist2, $q2, $m2) = &splitDistance($qpTree1,$qa2);
							$qqdist += ($spdist1 + $spdist2) / 2;
							
						}
						elsif ( ($nqqq1 == 2) && ($nqqq2 == 1)) {
							$qpTree1=~ /(.+)\-(.+)/;
							my $qa1 = $1;
							my $qa2 = $2;
							my ($spdist1, $q1, $m1) = &splitDistance($qpTree2,$qa1);
							my ($spdist2, $q2, $m2) = &splitDistance($qpTree2,$qa2);
							$qqdist += ($spdist1 + $spdist2) / 2;
							
						}
						else {
							my ($spdist, $q, $m) = &splitDistance($qpTree1,$qpTree2);
							$qqdist += $spdist;
							
						}
						$ndist++;
	
	
	
						
						if ( ($nqqq1 == 3) && ($nqqq2 == 3) ) {
							$qqdist2 += 0;
						}elsif ( ( ($nqqq1 == 3) && ($nqqq2 == 2) ) || ( ($nqqq1 == 2) && ($nqqq2 == 3) ) ) {
							$qqdist2 += 0;
						}
						elsif ( ( ($nqqq1 == 3) && ($nqqq2 == 1) ) || ( ($nqqq1 == 1) && ($nqqq2 == 3) ) ) {
							$qqdist2 += 0;
						}
						elsif ( ($nqqq1 == 2) && ($nqqq2 ==2)) {
							$qqdist2 += 0;
						}
						elsif ( ($nqqq1 == 1) && ($nqqq2 ==2)) {
							$qpTree2 =~ /(.+)\-(.+)/;
							my $qa1 = $1;
							my $qa2 = $2;
							my ($spdist1, $q1, $m1) = &splitDistance($qpTree1,$qa1);
							my ($spdist2, $q2, $m2) = &splitDistance($qpTree1,$qa2);
							if ( ($spdist1 == 0) || ($spdist2 == 0)  ) {
								$qqdist2 += 0;
							}else {
								$qqdist2 += 1;
							}
							
						}
						elsif ( ($nqqq1 == 2) && ($nqqq2 == 1)) {
							$qpTree1=~ /(.+)\-(.+)/;
							my $qa1 = $1;
							my $qa2 = $2;
							my ($spdist1, $q1, $m1) = &splitDistance($qpTree2,$qa1);
							my ($spdist2, $q2, $m2) = &splitDistance($qpTree2,$qa2);
							if ( ($spdist1 == 0) || ($spdist2 == 0)  ) {
								$qqdist2 += 0;
							}else {
								$qqdist2 += 1;
							}
						}
						else {	
							my ($spdist, $q, $m) = &splitDistance($qpTree1,$qpTree2);
							$qqdist2 += $spdist;
						}
						$ndist2++;
	
					}
					
					
					$qqdist /= $ndist;
					
	
					
					$qqdist2 /= $ndist2;
					
	
					printf "$nametrees1\t$nametrees2\t$ndist\t%6.3f\t%6.3f\n", $qqdist, $qqdist2;	
					print OUTqq "######################### topd $nametrees1 - $nametrees2 #########################\n";
					print OUTqq "* Quartets: $ndist\n";
					printf OUTqq "* qqdist: %6.3f\n", $qqdist;
					printf OUTqq "* qqdist_min: %6.3f\n", $qqdist2;
					
				}else {
					printf "$nametrees1\t$nametrees2\t-\t-\t-\n", $qqdist, $qqdist2;
					print OUTqq "######################### topd $nametrees1 - $nametrees2 #########################\n";
					print OUTqq "* Quartets: $ndist\n";
					printf OUTqq "* qqdist: %6.3f\n", $qqdist;
					printf OUTqq "* qqdist_min: %6.3f\n", $qqdist2;
					
				}
			}
		}
	}

	print "\nThanks for using TOPD-fMtS .................................. Bye.\n\n";
	
	close OUTqq;
}






sub duplidist {
	my $LocTrees = shift;
	my %name_tree = &getTrees($LocTrees);

	
	&SP($LocTrees);

	
	my %tree_species = ();
	my $nametree = "";
	my $speciestree = "";
	open (IN, "<species") || die "Cannot open species\n";
	while (<IN>) {
		chomp;
		if (/^>(.+)\:/) {
			$nametree = $1; print "$nametree\n";
		}
		else {
			$speciestree = $_;
			$tree_species{$nametree} = $speciestree;
		}	
	}
	close IN;

	
	open (OUTdd, ">dd_result");
	foreach $nametrees1(sort keys(%tree_species)) {
		my %ALLSP = ();
		my @speciestree1 = split /\t/, $tree_species{$nametrees1};
		my $nsp1 = @speciestree1;
		my %sptree1 = ();
		foreach $aux(@speciestree1) {
			if ($aux ne "") {
				$aux =~ /(.+) \((.+)\)/;
				$sptree1{$1} = $2;
				$ALLSP{$1} = "y";
				
			}
		}

		foreach $nametrees2(sort keys(%tree_species)) {
			if ($nametrees1 ne $nametrees2) {

				my @speciestree2 = split /\t/,  $tree_species{$nametrees2};
				my $nsp2 = @speciestree2;
				my %sptree2 = ();
				foreach $aux(@speciestree2) {
					if ($aux ne "") {
						$aux =~ /(.+) \((.+)\)/;
						$sptree2{$1} = $2;
						$ALLSP{$1} = "y";
						
					}
				}

				print "######################### topd $nametrees1 - $nametrees2 #########################\n";
				print "* taxa (ntree1) (ntree2): ";
				print OUTdd "######################### topd $nametrees1 - $nametrees2 #########################\n";
				print OUTdd "* Taxa (ntree1) (ntree2): ";
				my $dist12 = "";
				my $DIST12 = "";
				my $contCOMP = "";
				foreach $t(sort keys(%ALLSP)) {
					if ( ($sptree1{$t} ne "") && ($sptree2{$t} ne "") ) {
						print "$t ($sptree1{$t}) ($sptree2{$t}) | ";
						print OUTdd "$t ($sptree1{$t}) ($sptree2{$t}) | ";

						$dist12 = ($sptree1{$t} - $sptree2{$t})**2;
						$DIST12 += $dist12;
						$contCOMP++;
					}elsif ( ($sptree1{$t} ne "") && ($sptree2{$t} eq "") ) {
						print "$t ($sptree1{$t}) (0) | ";
						print OUTdd "$t ($sptree1{$t}) (0) | ";

						$dist12 = ($sptree1{$t} - $sptree2{$t})**2;
						$DIST12 += $dist12;
						$contCOMP++;
					}elsif ( ($sptree1{$t} eq "") && ($sptree2{$t} ne "") ) {
						print "$t (0) ($sptree2{$t}) | ";
						print OUTdd "$t (0) ($sptree2{$t}) | ";

						$dist12 = ($sptree1{$t} - $sptree2{$t})**2;
						$DIST12 += $dist12;
						$contCOMP++;
					}
				}
	
				
				my $duplidist = sqrt($DIST12 / $contCOMP);
				printf "\n* Duplidist: %6.3f\n", $duplidist;
				printf OUTdd"\n* Duplidist: %6.3f\n", $duplidist;
			}
		}
	}

	print "\nThanks for using TOPD-fMtS .................................. Bye.\n\n";
	
	close OUTdd;
}






sub credits() {
	print "\n";
	print "##########################################\n";
	print "#                                        #\n";
	print "#      TOPD-fMtS version 3.3             #\n";
	print "#      (August 2007)                     #\n";
	print "#                                        #\n";
	print "#      Author: Pere Puigbo               #\n";
	print "#      E-mail: ppuigbo\@urv.net           #\n";
	print "#      http://genomes.urv.es/topd        #\n";
	print "#                                        #\n";
	print "##########################################\n\n";
}

