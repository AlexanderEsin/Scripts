#!/usr/bin/Rscript
### For each transfer prediction into Geobacillus, we identify the Geobacillus tips involved in that transfer. This is done for each group, at each Mowgli penalty ##

library(ape)
library(phylobase)
library(phangorn)
library(stringr)
library(gtools)

###########################################################################

penalty_list <- c(4, 5, 6, 8, 10, 20)
master_dir <- "/users/aesin/Desktop/Mowgli/Mowgli_outputs"
outpit_dir <- paste0(master_dir, "/Per_penalty_tips")
dir.create(outpit_dir, showWarnings = F)

###########################################################################

for (penalty in penalty_list) {
  penalty_dir <- (paste0(master_dir, "/Mow_test_out_t", penalty))
  setwd(penalty_dir)
  directories <- mixedsort(dir())
  
  data_table <- matrix(nrow = length(directories)*3, ncol = 4)
  table_entry_counter = 1
  directory_counter = 1
  
  for (directory in directories) {
    folder_dir <- as.character(paste0(penalty_dir, "/", directory))
    setwd(folder_dir)

    # Get the value of the nodes from the Transfers_in.tsv file #
    transfer_data <- as.character(read.delim("Transfers_in.tsv", header = F, sep = "\n")[,1])
    hgt_lines <- grep("HGT", transfer_data, value = TRUE)
    hgt_lines <- str_replace(hgt_lines, pattern = "\t", replacement = " ")
    
    # If there are no HGTs - continue #
    if (length(hgt_lines) == 0) {
      next
    }
    
    # Process the lines defining at which nodes the transfers occurred; extract just the numbers #
    node_lines <- grep("node", transfer_data, value = TRUE)
    node_lines_split <- unlist(strsplit(node_lines, " "))
    nodes <- grep("([0-9]+)", node_lines_split, value = TRUE)
    
    # Get the subset of Geobacillus tips given the node at which transfer occurred #
    full_gene_tree <- read.tree("FullGeneTree.mpr")
    full_gt4 <- phylo4(full_gene_tree)
    
    tips_needed <- list(); node_counter = 1
    for (node in nodes) {
      child_tips <- names(descendants(full_gt4, node = node))
      # If we can't find any tips, that means our transfer was into a terminal edge - i.e. a tip, so we need to find the tip corresponding to that node value #
      if (length(child_tips) == 0) {
        pattern = paste0("(_", node, "$)+")
        child_tips <- grep(pattern, as.character(tipLabels(full_gt4)), value = TRUE)
      }
      
      tips_needed[[node_counter]] <- child_tips
      node_counter = node_counter + 1
    }
    
    # Keep only the Geobacillus tips and strip away the mowgli number at the end of the tip names #
    toMatch <- c("Geobacillus", "Anoxybacillus")
    anoxy_geo_tips <- lapply(tips_needed, grep, pattern = paste(toMatch, collapse = "|"), value = TRUE)
    anoxy_geo_tips <- anoxy_geo_tips[lapply(anoxy_geo_tips, length)>0]
    ageo_tips_trim <- lapply(anoxy_geo_tips, str_replace, pattern = "(_[0-9]+$)", replacement = "")
    
    # Write out the list of species for each transfer event for each group #
    for (n in 1:length(ageo_tips_trim)) {
      data_table[table_entry_counter,] <- c(directory, hgt_lines[n], nodes[n], paste(ageo_tips_trim[[n]], collapse = " ")) 
      table_entry_counter = table_entry_counter + 1
    }

    # Counter #
    message(directory_counter)
    directory_counter = directory_counter + 1
    
  }

  data_df <- as.data.frame(data_table)
  data_df_natrim <- data_df[which(is.na(data_df$V1) == FALSE),]
  colnames(data_df_natrim) <- c("Group", "HGT", "Gene_tree_node", "Species")

  setwd(outpit_dir)
  write.table(data_df_natrim, file = paste0("Per_penalty_tips_t", penalty, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

}



###########################################################################






