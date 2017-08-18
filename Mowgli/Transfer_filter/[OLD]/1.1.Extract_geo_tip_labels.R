#!/usr/local/bin/Rscript

### For each transfer prediction into Geobacillus, we identify the Geobacillus tips involved in that transfer. This is done for each group, at each Mowgli penalty ##

## 05/2017 - Completely reworked the script to use lapply. Importantly, in some rare cases HGTs into Anoxy/Geobacillus were nested within other transfers. This meant that for a certain gene family (e.g. 1872, T=5) the same tips were counted as part of more than one transfer ##

library(ape)
library(phylobase)
library(phangorn)
library(stringr)
library(gtools)

GetTerminalTip	<- function(node, tree) {
	 pattern = paste0("(_", node, "$)+")
	 terminal_tip <- grep(pattern, as.character(tipLabels(tree)), value = TRUE)
	 return(terminal_tip)
}

penalty_list	<- c(3, 4, 5, 6)
include_AG2AG	<- TRUE

master_dir <- "/Users/aesin/Desktop/Mowgli/Mowgli_outputs"
if (include_AG2AG == TRUE) {
	output_dir <- paste0(master_dir, "/Per_penalty_tips/AG2AG_incl/")	
} else {
	output_dir <- paste0(master_dir, "/Per_penalty_tips/AG2AG_excl/")	
}
dir.create(output_dir, showWarnings = FALSE)

##

NodeLabPattern = "[0-9]+$"
tip_number_tbl	<- data.frame(Penalty = numeric(), Number.Of.Tips = numeric())

for (penalty in penalty_list) {
	directories <- mixedsort(dir(paste0(master_dir, "/Mow_test_out_t", penalty)))
	  
	## WARNINGS ABOUT NODES NOT FOUND REFER TO INABILITY TO FIND NODES WHICH ARE IN FACT TIPS (SOLO TRANSFER TO ONE SPECIES) ##
	per_group_HGT_tips	<- lapply(directories, function(directory) {
		message(paste0("Penalty: ", penalty, " == Directory: ", directory))
		dir <- paste0(master_dir, "/Mow_test_out_t", penalty, "/", directory)

		# Read in transfer data from the Parsed_events.tsv file #
		transfer_data	<- read.table(paste0(dir, "/Parsed_events.tsv"), header = TRUE, sep = "\t")
		external_trans	<- transfer_data[transfer_data$HGT.Type == "OutAG", ]
		internal_trans	<- transfer_data[transfer_data$HGT.Type == "AG2AG", ]

		OutAG_nodes		<- external_trans$GeneT.Child

		# If there are no external transfers detected, then the fate of all AG species is either Vertical or Root
		if (nrow(external_trans) == 0) {
			# message("No external transfers == SKIP")
			return()
		}

		# Read in the reconciled gene tree #
		full_gt4 <- phylo4(read.tree(paste0(dir, "/FullGeneTree.mpr")))

		HGT_transfer_tips	<- lapply(1:nrow(external_trans), function(OutAG_index) {
			OutAG			<- external_trans[OutAG_index, ]
			index			<- as.character(OutAG$Transfer.Index)
			genet_receiver	<- as.character(OutAG$GeneT.Child)

			# Get all the nodes descendent of this OutAG receiver node
			child_nodes_raw	<- names(descendants(full_gt4, node = genet_receiver, type = "all"))

			# Check whether another OutAG receiver node is 'internalised' in this transfer
			internal_xfer	<- child_nodes_raw[str_extract(child_nodes_raw, NodeLabPattern) %in% OutAG_nodes]
			if (length(internal_xfer) > 0) {
				# Exclude all tips belonging to the secondary transfer(s) from the full gene tree
				excluded_tips	<- unlist(lapply(internal_xfer, function(int_node) names(descendants(full_gt4, node = int_node, type = "tips"))))
				final_tree		<- subset(full_gt4, tips.exclude = excluded_tips)
			} else {
				final_tree		<- full_gt4
			}

			# Check for any AG2AG events derived from this OutAG
			extra_tips_list		<- vector(mode = "list")
			derived_AG2AG_list	<- internal_trans[internal_trans$Heritage == index,]
			if (length(derived_AG2AG_list) != 0) {
				for (index in 1:nrow(derived_AG2AG_list)) {
					derived_AG2AG		<- derived_AG2AG_list[index, ]
					derived_receiver	<- as.character(derived_AG2AG$GeneT.Child)
					
					# If the AG2AG node already lies within the OutAG clade...
					if (length(child_nodes_raw[str_extract(child_nodes_raw, NodeLabPattern) %in% derived_receiver]) != 0) {
						# If the derived AG2AG transfer nodes are already within the OutAG transfer nodes and we want them included, continue
						if (include_AG2AG == TRUE) {
							next
						# If we don't want to include tips that derive from secondary AG2AG transfers, remove them from the tree
						} else {
							excluded_tips	<- names(descendants(full_gt4, node = derived_receiver, type = "tips"))
						}
						
					# If the AG2AG node derived from an OutAG transfer does not lie within the OutAG gene tree clade...
					} else {
						# If we want them included..
						if (include_AG2AG == TRUE) {
							extra_tips		<-  names(descendants(final_tree, node = derived_receiver, type = "tips"))
							if (length(extra_tips) == 0) {extra_tips	<- GetTerminalTip(derived_receiver, final_tree)}
							extra_tips_list	<- c(extra_tips_list, extra_tips)
						# Otherwise just continue
						} else {
							next
						}
					}
				}
			}
			

			core_tips  <- names(descendants(final_tree, node = genet_receiver, type = "tips"))
			if (length(core_tips) == 0) {
				core_tips	<- GetTerminalTip(genet_receiver, final_tree)
			}
			all_tips	<- c(core_tips, extra_tips_list)

			toMatch <- c("Geobacillus", "Anoxybacillus")
			anoxygeo_tips	<- grep(all_tips, pattern = paste(toMatch, collapse = "|"), value = TRUE)
			anoxygeo_tips_trim	<- str_replace(anoxygeo_tips, pattern = "(_[0-9]+$)", replacement = "")

			return(data.frame(Directory = directory, Donor_edge = OutAG$Donor.Edge, Receiver_edge = OutAG$Receiver.Edge, Tip_number = length(anoxygeo_tips_trim), Species = paste(anoxygeo_tips_trim, collapse = " ")))

		})
		return(do.call(rbind.data.frame, HGT_transfer_tips))
	})

	per_group_HGT_tips	<- per_group_HGT_tips[lapply(per_group_HGT_tips, length) > 0]
	per_penalty_df		<- do.call(rbind.data.frame, per_group_HGT_tips)

	total_num_tips		<- sum(per_penalty_df$Tip_number)
	tip_number_tbl[nrow(tip_number_tbl)+1,]	<- c(penalty, total_num_tips)

	write.table(per_penalty_df, file = paste0(output_dir, "Per_penalty_tips_t", penalty, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

write.table(tip_number_tbl, file = paste0(output_dir, "Stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)