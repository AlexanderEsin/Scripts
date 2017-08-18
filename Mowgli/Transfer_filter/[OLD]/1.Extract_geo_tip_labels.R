#!/usr/local/bin/Rscript

### For each transfer prediction into Geobacillus, we identify the Geobacillus tips involved in that transfer. This is done for each group, at each Mowgli penalty ##

## 05/2017 - Completely reworked the script to use lapply. Importantly, in some rare cases HGTs into Anoxy/Geobacillus were nested within other transfers. This meant that for a certain gene family (e.g. 1872, T=5) the same tips were counted as part of more than one transfer ##

library(ape)
library(phylobase)
library(phangorn)
library(stringr)
library(gtools)

###########################################################################

penalty_list <- c(3, 4, 5, 6)

master_dir <- "/Users/aesin/Desktop/Mowgli/Mowgli_outputs"
output_dir <- paste0(master_dir, "/Per_penalty_tips/")
dir.create(output_dir, showWarnings = FALSE)

###########################################################################
tmp <- lapply(penalty_list, function(penalty) {

  directories <- mixedsort(dir(paste0(master_dir, "/Mow_test_out_t", penalty)))
    
  ## WARNINGS ABOUT NODES NOT FOUND REFER TO INABILITY TO FIND NODES WHICH ARE IN FACT TIPS (SOLO TRANSFER TO ONE SPECIES) ##
  per_group_HGT_tips	<- lapply(directories, function(directory) {
    dir <- paste0(master_dir, "/Mow_test_out_t", penalty, "/", directory)
  	setwd(dir = dir)

  	# Read in transfer data from the Transfers_in.tsv file #
  	transfer_data	<- read.table("Transfers_in.tsv", header = FALSE, sep = "\n")[-1,]
  	num_events		<- length(transfer_data) / 2

  	# Extract the donor/receptor node data #
  	hgt_lines <- grep("HGT", transfer_data, value = TRUE)
  	hgt_lines <- str_replace(hgt_lines, pattern = "\t", replacement = " ")

  	# If there are no HGTs - continue #
  	if (length(hgt_lines) == 0) {
  		return()
  	}

  	# Extract the gene tree transfer node #
  	node_data_split <- unlist(strsplit(grep("node", transfer_data, value = TRUE), " "))
  	nodes <- grep("([0-9]+)", node_data_split, value = TRUE)

  	# Read in the reconciled gene tree #
  	full_gt4 <- phylo4(read.tree("FullGeneTree.mpr"))

  	HGT_transfer_tips	<- lapply(seq(num_events), function(index) {
  		hgt_line	<- hgt_lines[index]
  		node		<- nodes[index]

      pattern = "[0-9]+$"
  		all_child_nodes	<- names(descendants(full_gt4, node = node, type = "all"))
  		internal_xfer <- all_child_nodes[str_extract(all_child_nodes, pattern) %in% nodes]

  		if (length(internal_xfer) > 0) {
  			print(paste0("A transfer is nested in another transfer. Gene family: ", directory))
  			excluded_tips <- unlist(lapply(internal_xfer, function(int_node) names(descendants(full_gt4, node = int_node, type = "tips"))))
  			subset_tree   <- subset(full_gt4, tips.exclude = excluded_tips)
  			all_child_tips  <- names(descendants(subset_tree, node = node, type = "tips")) 
  		} else {
  			all_child_tips  <- names(descendants(full_gt4, node = node, type = "tips"))

  			# If we can't find any tips, that means our transfer was into a terminal edge - i.e. a tip, so we need to find the tip corresponding to that node value #
  			if (length(all_child_tips) == 0) {
  			  pattern = paste0("(_", node, "$)+")
  			  all_child_tips <- grep(pattern, as.character(tipLabels(full_gt4)), value = TRUE)
  			}
  		}

  		toMatch <- c("Geobacillus", "Anoxybacillus")
  		anoxygeo_child_tips	<- grep(all_child_tips, pattern = paste(toMatch, collapse = "|"), value = TRUE)
  		anoxygeo_child_tips_trim	<- str_replace(anoxygeo_child_tips, pattern = "(_[0-9]+$)", replacement = "")

      message(paste0("Penalty: ", penalty, "-- Directory: ", directory))

  		# Write out the list of species for each transfer event for each group #
  		return(data.frame(Directory = directory, HGT = hgt_line, Gene_tree_node = node, Species = paste(anoxygeo_child_tips_trim, collapse = " ")))
  	})

  	return(do.call(rbind.data.frame, HGT_transfer_tips))
  })

  per_group_HGT_tips <- per_group_HGT_tips[lapply(per_group_HGT_tips, length) > 0]
  per_penalty_df  <- do.call(rbind.data.frame, per_group_HGT_tips)
  write.table(per_penalty_df, file = paste0(output_dir, "Per_penalty_tips_t", penalty, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  return()
})



###########################################################################






