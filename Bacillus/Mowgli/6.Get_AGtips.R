#!/usr/local/bin/Rscript

### For each transfer prediction into Geobacillus, we identify the Geobacillus tips involved in that transfer. This is done for each group, at each Mowgli penalty ##

## 05/2017 - Completely reworked the script to use lapply. Importantly, in some rare cases HGTs into Anoxy/Geobacillus were nested within other transfers. This meant that for a certain gene family (e.g. 1872, T=5) the same tips were counted as part of more than one transfer ##

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load(tidyverse, fs, parallel, ape, gdata, phylobase, phangorn, gtools)

GetTerminalTip	<- function(node, tree) {
	pattern	<- paste0("^(", node, "_)+")
	terminal_tip	<- str_subset(as.character(tipLabels(tree)), pattern)
	if (length(terminal_tip) > 1) stop("Should not have more than one terminal tip")
	return(terminal_tip)
}

penalty_list	<- c(3, 4)
include_core	<- TRUE

master		<- "/Users/aesin/Desktop/Bacillus"
mowOut_dir	<- file.path(master, "Mowgli/Mowgli_output")
if (include_core == TRUE) {
	output_dir <- file.path(mowOut_dir, "Per_penalty_tips", "Core_incl")	
} else {
	output_dir <- file.path(mowOut_dir, "Per_penalty_tips", "Core_excl")	
}
dir_create(output_dir)


coreKeep_file	<- file.path(master, "Core_genomes", "Genome_lists", "coreToKeep.tsv")
coreKeep_tbl	<- read_tsv(coreKeep_file)

##

NodeLabPattern = "[0-9]+$"
tip_number_tbl	<- data.frame(Penalty = numeric(), Number.Of.Tips = numeric())

for (penalty in penalty_list) {

	message(paste0("Identifying Core tips per group for penalty: ", penalty, "..."), appendLF = FALSE)

	directories <- mixedsort(dir(file.path(mowOut_dir, paste0("Output_", penalty))))
	  
	## WARNINGS ABOUT NODES NOT FOUND REFER TO INABILITY TO FIND NODES WHICH ARE IN FACT TIPS (SOLO TRANSFER TO ONE SPECIES) ##
	per_group_HGT_tips	<- mclapply(directories, function(directory) {
		# message(paste0("Penalty: ", penalty, " == Directory: ", directory))
		dir <- file.path(mowOut_dir, paste0("Output_", penalty), directory)

		# Read in transfer data from the Parsed_events.tsv file #
		if (!file.exists(file.path(dir, "Parsed_events.tsv"))) {
			return()
		}

		transfer_data	<- read.table(file.path(dir, "Parsed_events.tsv"), header = TRUE, sep = "\t")
		external_trans	<- transfer_data[transfer_data$HGT.Type == "OutAG", ]
		internal_trans	<- transfer_data[transfer_data$HGT.Type == "AG2AG", ]

		OutAG_nodes		<- external_trans$GeneT.Child

		# If there are no external transfers detected, then the fate of all core species is either Vertical or Root
		if (nrow(external_trans) == 0) {
			# message("No external transfers == SKIP")
			return()
		}

		# Read in the reconciled gene tree #
		full_gt4 <- phylo4(read.tree(file.path(dir, "FullGeneTree.mpr")))

		HGT_transfer_tips	<- lapply(1:nrow(external_trans), function(OutAG_index) {
			OutAG			<- external_trans[OutAG_index, ]
			index			<- as.character(OutAG$Transfer.Index)
			genet_receiver	<- as.character(OutAG$GeneT.Child)

			## Get all the nodes descendent of this OutAG receiver node
			child_nodes_raw	<- names(descendants(full_gt4, node = genet_receiver, type = "all"))

			## Check whether another OutAG receiver node is 'internalised' in this transfer
			## New approach to find nested nodes because terminal nodes are now given as
			## binomial names instead of node IDs
			find_nested_nodes	<- lapply(OutAG_nodes, function(OutAG_node) {

				## If we find that an OutAG node contains an underscore, it must be a
				## terminal tip so we use a pattern allowing us to catch the
				## tip together with it's trailing ID. E.g.:
				## 471223_1 + _1905

				if (length(grep("[_]+", OutAG_node)) != 0) {
					pattern	<- paste0("(^", OutAG_node, "_)+")

				## If there is no underscore, the OutAG node is an internal branch and we want to 
				## find an exact copy (i.e. we don't want to find "198" when the node)
				## we are looking for is "19"
				} else {
					pattern = paste0("(^", OutAG_node, "$)+")
				}
				
				if (length(grep(pattern, child_nodes_raw, value = TRUE)) != 0) {
					return(as.character(OutAG_node))
				}
			})
			nested_transfers	<- unlist(Filter(length, find_nested_nodes))


			if (length(nested_transfers) == 0) {
				nested_transfers 	<- NA
			} else {
				nested_transfers	<- paste0(nested_transfers, collapse = " ")
			}
			
			# Check for any AG2AG events derived from this OutAG
			extra_tips_list		<- vector(mode = "list")
			derived_AG2AG_list	<- internal_trans[internal_trans$Heritage == index,]
			if (nrow(derived_AG2AG_list) != 0) {
				for (index in 1:nrow(derived_AG2AG_list)) {
					derived_AG2AG		<- derived_AG2AG_list[index, ]
					derived_receiver	<- as.character(derived_AG2AG$GeneT.Child)
					
					# If the AG2AG node already lies within the OutAG clade...
					if (length(grep(derived_receiver, child_nodes_raw, value = TRUE)) != 0) {
						# If the derived AG2AG transfer nodes are already within the OutAG transfer nodes and we want them included, continue
						if (include_core == TRUE) {
							next
						# If we don't want to include tips that derive from secondary AG2AG transfers, remove them from the tree
						} else {
							excluded_tips	<- names(descendants(full_gt4, node = derived_receiver, type = "tips"))
						}
						
					# If the AG2AG node derived from an OutAG transfer does not lie within the OutAG gene tree clade...
					} else {
						# If we want them included..
						if (include_core == TRUE) {
							extra_tips		<-  names(descendants(full_gt4, node = derived_receiver, type = "tips"))
							if (length(extra_tips) == 0) {extra_tips	<- GetTerminalTip(derived_receiver, full_gt4)}
							extra_tips_list	<- c(extra_tips_list, extra_tips)
						# Otherwise just continue
						} else {
							next
						}
					}
				}
			}
			

			core_tips  <- names(descendants(full_gt4, node = genet_receiver, type = "tips"))
			if (length(core_tips) == 0) {
				core_tips	<- GetTerminalTip(genet_receiver, full_gt4)
			}
			all_tips	<- c(core_tips, extra_tips_list)

			anoxygeo_tips	<- grep(pattern = paste(paste0("^", coreKeep_tbl$Taxid, "_"), collapse = "|"), all_tips, value = TRUE) 
			anoxygeo_tips_trim	<- str_replace(anoxygeo_tips, pattern = "(_[0-9]+$)", replacement = "")

			## A hack: nesting can be extensive. E.g. check 1262 (19-set) @ T5. We have an OutAG -> AG2AG -> AG2AG
			## but the second AG2AG is nested within first, leading to double naming. Instead of trying to deal with
			## every eventuality, just take the unique names.
			anoxygeo_tips_trim_unique	<- unique(anoxygeo_tips_trim)
			if (length(anoxygeo_tips_trim_unique) != length(anoxygeo_tips_trim)) {
				message(paste0("Penalty: ", penalty, " == Directory: ", directory))
			}

			AGtips_char	<- paste0(" ", paste(anoxygeo_tips_trim_unique, collapse = " "), " ")
			out_df		<- data.frame(Directory = directory, Donor_edge = OutAG$Donor.Edge, Receiver_edge = OutAG$Receiver.Edge, GeneT_child_node = genet_receiver, Tip_number = length(anoxygeo_tips_trim_unique), Species = AGtips_char, Nested_T = nested_transfers, stringsAsFactors = FALSE)

			return(out_df)

		})


		## In some cases, we have nested nodes with secondary AG2AG transfers (that are also the nested) but are not removed in the above steps. We corrected for that here
		corrected_tips_out	<- lapply(HGT_transfer_tips, function(element) {
			if (is.na(element$Nested_T) == FALSE) {

				tip_list		<- unlist(str_split(str_trim(element$Species), " "))
				nested_nodes	<- unlist(str_split(element$Nested_T, " "))

				# For each nest node, exclude the tips belonging to the nested node from the nestING node tips
				for (node in nested_nodes) {
					tips_to_exclude	<- unlist(lapply(HGT_transfer_tips, function(element2) {
						if (element2$GeneT_child_node == node) {
							tips	<- unlist(str_split(str_trim(element2$Species), " "))
							return(tips)
						}
					}))
					tip_list		<- tip_list[!tip_list %in% tips_to_exclude]
				}
				# Replace original tip list with the derived tip list
				element$Species		<- paste0(" ", paste(tip_list, collapse = " "), " ")
				element$Tip_number	<- length(tip_list)
			}
			return(element)
		})
		full_table	<- do.call(rbind.data.frame, corrected_tips_out)
		trunc_table	<- full_table[,-ncol(full_table)]
		return(trunc_table)
	}, mc.cores = 10)

	per_group_HGT_tips	<- per_group_HGT_tips[lapply(per_group_HGT_tips, length) > 0]
	per_penalty_df		<- do.call(rbind.data.frame, per_group_HGT_tips)
	# Filter out phantom transfers
	numPhantom			<- nrow(per_penalty_df %>% filter(Tip_number == 0))
	per_penalty_df		<- per_penalty_df %>% filter(Tip_number != 0)

	total_num_tips		<- sum(per_penalty_df$Tip_number)
	tip_number_tbl[nrow(tip_number_tbl)+1,]	<- c(penalty, total_num_tips)

	write.table(per_penalty_df, file = file.path(output_dir, paste0("Per_penalty_tips_t", penalty, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

	message(paste0("\rIdentifying core tips per group for penalty: ", penalty, "... done"))
	message(paste0("Number of phantom transfers for penalty: ", penalty, " == ", numPhantom))
}

write.table(tip_number_tbl, file = file.path(output_dir, "Stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)