library(ape)
library(phylobase)
library(phytools)
library(stringr)
library(geiger)
library(RSvgDevice)
library(reshape2)
library(ggplot2)

###########################################################################
## Annotate branches ##
annotate_transfer_num <- function (edge_entry, index) {
    edge <- edge_entry[1:2]
    edge_label <- as.character(round(edge_entry[3], digits = 2))
    # print(edge_label)

    edgelabels(edge_label, index, adj = c(0.5, -0.25), bg = "white", frame = "none", cex = 0.7)
}

###########################################################################
## Variables ##
penalties = c(4:6)
min_taxa = 0; # 0 or 50

single_pen_plot = TRUE
single_pen_val = 5

input_edge_dir <- "/Users/aesin/Desktop/Mowgli/Long_distance_HGT/Full/Scenarios_1_2/Receptor_edges/"

if (min_taxa > 4) {
	file_name_add <- paste0("_min_", min_taxa)
} else {
	file_name_add <- ""
}

###########################################################################

##############################
# Prepare the cladogram tree #
##############################

species_tree_dir = "/Users/aesin/Desktop/Mowgli/Species_tree/Reconciled_species_tree"
setwd(species_tree_dir)

species_tree <- read.tree("outputSpeciesTree.mpr")

geo_tips <- grep("Geobacillus_", species_tree$tip.label, value = TRUE)
anoxy_tips <- grep("Anoxybacillus_", species_tree$tip.label, value = TRUE)
anoxy_geo_species_tips <- c(geo_tips, anoxy_tips)

anoxy_geo_species <- drop.tip(species_tree, setdiff(species_tree$tip.label, anoxy_geo_species_tips))
anoxy_geo_4 <- phylo4(anoxy_geo_species)

##################################
## Get the node labels as a data frame ##
node_lbs_df <- data.frame(labels(anoxy_geo_4, "all"))
internal_nodes_species <- nodeId(anoxy_geo_4, "internal")
names(node_lbs_df)[1] <- "Labels"

## Prepare a data frame with the edges so we can count ##
anoxy_geo_edge_df <- data.frame(anoxy_geo_species$edge)
names(anoxy_geo_edge_df) <- c("E1", "E2")

## Make a column for each penalty ##
for (penalty in penalties) {
    cname <- paste0("T", penalty)
    anoxy_geo_edge_df$colblank <- 0
    names(anoxy_geo_edge_df)[names(anoxy_geo_edge_df)=="colblank"] <- cname
}

###########################################################################

##################################
# Prepare the time-resolved tree #
##################################

## Read in the time-resolved tree ##
bacillaceae <- read.tree(file = "/Users/aesin/Desktop/Consensus_trees/Consensus_finals/Bacillaceae_N1000_newick.txt")

anoxy_tips <- grep("Anoxybacillus_", bacillaceae$tip.label, value = TRUE)
geo_tips <- grep("Geobacillus_", bacillaceae$tip.label, value = TRUE)
anoxy_geo_time_tips <- c(geo_tips, anoxy_tips)

anoxy_geo_time <- drop.tip(bacillaceae, setdiff(bacillaceae$tip.label, anoxy_geo_time_tips))

##################################
## Annotate the time-resolved tree with the correct node labels imported from the species tree-derived tree ##
anoxy_geo_time_relab <- anoxy_geo_time
for (tip in anoxy_geo_time_tips) {
    species_tree_tip_name <- grep(tip, anoxy_geo_species_tips, value = T)
    original_tip_position <- match(tip, anoxy_geo_time_relab$tip.label)
    # Replace with new name #
    anoxy_geo_time_relab$tip.label[original_tip_position] <- species_tree_tip_name
}

## Remove node labels (BS-value) and convert to phylo4##
anoxy_geo_time_relab$node.label <- NULL
anoxy_geo_time_relab_4 <- phylo4(anoxy_geo_time_relab)

## Get all the internal nodes
internal_nodes_time <- nodeId(anoxy_geo_time_relab_4, "internal")
unmapped_edges <- vector(mode = "list")
unmapped_count = 1

for (internal_time in internal_nodes_time) {
    matched = FALSE
    for (internal_species in internal_nodes_species) {
        if (setequal(tips(anoxy_geo_time_relab, internal_time), tips(anoxy_geo_species, internal_species)) == TRUE) {
            ## Get the label of the node from the species-cladorgram tree ##
            node_label <- anoxy_geo_4@label[internal_species]

            ## Index of the internal node label in the time-resolved tree ##
            node_index <- which(names(nodeLabels(anoxy_geo_time_relab_4)) == internal_time)
            
            ## Assign the name of the node to the time-resolved tree ##
            nodeLabels(anoxy_geo_time_relab_4)[node_index] <- node_label

            print(paste0("Time tree node: ", internal_time, " maps to: ", internal_species))
            matched = TRUE
            break
        }
    }
    if (matched == FALSE) {
        print(paste0("Time tree node: ", internal_time, " did not map to any corresponding species tree node"))
        edge_index <- which(anoxy_geo_time_relab_4@edge[,2] == as.character(internal_time))
        unmapped_edges[[unmapped_count]] <- as.vector(anoxy_geo_time_relab_4@edge[edge_index,])
        unmapped_count = unmapped_count + 1
    }
}


##################################
## Prepare the count table for the time-resolved tree as for the species-derived prune above ##
node_lbs_time_df <- data.frame(labels(anoxy_geo_time_relab_4, "all"))
names(node_lbs_time_df)[1] <- "Labels"

## Prepare a data frame with the edges so we can count ##
anoxy_geo_time_edge_df <- data.frame(anoxy_geo_time$edge)
names(anoxy_geo_time_edge_df) <- c("E1", "E2")

## Make a column for each penalty ##
for (penalty in penalties) {
    cname <- paste0("T", penalty)
    anoxy_geo_time_edge_df$colblank <- 0
    names(anoxy_geo_time_edge_df)[names(anoxy_geo_time_edge_df)=="colblank"] <- cname
}


###########################################################################

for (penalty in penalties) {

    ## Read in the refined receptor edge list for a particular penalty ##
    receptor_df <- read.table(file = paste0(input_edge_dir, "T", penalty, "_refined_receptor_edges", file_name_add, ".txt"), header = FALSE, sep = "\t")

    ## Column name for this penalty ##
    cname <- paste0("T", penalty)

    ## For each transfer, we can identify the edge ##
    at_root = 0
    for(i in 1:nrow(receptor_df)) {
        row <- receptor_df[i,]

        if (row[1] == 3728 || row[2] == 3728) {
            at_root = at_root + 1
        } else {
            edge_2_species <- grep(row[2], node_lbs_df$Labels)

            edge_2_time <- grep(row[2], node_lbs_time_df$Labels)

            ## Fill in the count for the species-derived tree ##
            anoxy_geo_edge_df[,cname][which(anoxy_geo_edge_df$E2 == edge_2_species)] = anoxy_geo_edge_df[,cname][which(anoxy_geo_edge_df$E2 == edge_2_species)] + 1

            ## Fill in the count for the time-derived tree ##
            anoxy_geo_time_edge_df[,cname][which(anoxy_geo_time_edge_df$E2 == edge_2_time)] = anoxy_geo_time_edge_df[,cname][which(anoxy_geo_time_edge_df$E2 == edge_2_time)] + 1
        }
    }
}

###########################################################################
## Save a copy of each df with the raw numbers ##
anoxy_geo_edge_raw_df <- anoxy_geo_edge_df
anoxy_geo_time_edge_raw_df <- anoxy_geo_time_edge_df

## Make the transfers in at each branch a fraction of the total transfers ##
anoxy_geo_edge_df <- cbind(anoxy_geo_edge_df[,1:2],sweep(anoxy_geo_edge_df[,-1:-2],2,colSums(anoxy_geo_edge_df[,-1:-2]), "/")*100)
anoxy_geo_time_edge_df <- cbind(anoxy_geo_time_edge_df[,1:2],sweep(anoxy_geo_time_edge_df[,-1:-2],2,colSums(anoxy_geo_time_edge_df[,-1:-2]), "/")*100)

## Add a mean column ##
anoxy_geo_edge_df$mean <- rowMeans(anoxy_geo_edge_df[,-1:-2])
anoxy_geo_time_edge_df$mean <- rowMeans(anoxy_geo_time_edge_df[,-1:-2])

###########################################################################
## Perform correlations ##
x_corr <- cbind(anoxy_geo_time_edge_df[,1:2],anoxy_geo_time_edge_df$mean, anoxy_geo_time$edge.length)
overall_pearson <- round(cor(x_corr[,3:4], method = "pearson"), digits = 4)[2]
overall_spearman <- round(cor(x_corr[,3:4], method = "spearman"), digits = 4)[2]
writeLines(paste0("Pearson correlation of transfers over all branch lengths: ", overall_pearson, "\n", "Spearman correlation of transfers over all branch lengths: ", overall_spearman))


## Remove those taxa that are based on incomplete chromosome assembly ##
refined_x_corr <- x_corr
poor_assembly_taxa <- c("Geobacillus_caldoxylosilyticus", "Geobacillus_sp_G11MC16", "Geobacillus_sp_G1w1", "Geobacillus_thermocatenulatus", "Geobacillus_sp_WSUCF1", "Geobacillus_sp_MAS1")

for (taxon in poor_assembly_taxa) {
    node <- match(taxon, anoxy_geo_time$tip)
    print(node)
    refined_x_corr <- refined_x_corr[-which(refined_x_corr$E2 == node),]
}

contig_removed_pearson <- round(cor(refined_x_corr[,3:4], method = "pearson"), digits = 4)[2]
contig_removed_spearman <- round(cor(refined_x_corr[,3:4], method = "spearman"), digits = 4)[2]

writeLines(paste0("Pearson correlation of transfers over branches with contig-derived assembly branches removed: ", contig_removed_pearson, "\nSpearman correlation of transfers over branches with contig-derived assembly branches removed: ", contig_removed_spearman))

## Remove anoxybacillus taxa and supporting nodes ##
geo_only_with_edge <- x_corr[c(-1:-3,-5),]
geo_only <- geo_only_with_edge[,3:4]

geo_only_pearson <- round(cor(geo_only, method = "pearson"), digits = 4)[2]
geo_only_spearman <- round(cor(geo_only, method = "spearman"), digits = 4)[2]
writeLines(paste0("Pearson correlation of transfers over Geobacillus only branches: ", geo_only_pearson, "\nSpearman correlation of transers over Geobacillus only branches: ", geo_only_spearman))


## Keep only the branches that are terminal ##
terminal_branches_df <- data.frame(E1 = numeric(), E2 = numeric(), count = numeric(), edge_length = numeric())

no_anoxy_geo_corr <- x_corr[c(-1:-3,-5),]

for (taxon in geo_tips) {

    node <- match(taxon, anoxy_geo_time$tip)
    print(taxon)
    print(node)

    terminal_branches_df[nrow(terminal_branches_df)+1,] <- no_anoxy_geo_corr[which(no_anoxy_geo_corr$E2 == node),]
}

terminal_only_pearson <- round(cor(terminal_branches_df[,3:4], method = "pearson"), digits = 4)[2]
terminal_only_spearman <- round(cor(terminal_branches_df[,3:4], method = "spearman"), digits = 4)[2]

writeLines(paste0("Pearson correlation of transfers over terminal Geobacillus branches only: ", terminal_only_pearson, "\nSpearman correlation of transers over terminal Geobacillus branches only: ", terminal_only_spearman))

## Keep only branches that are NON-terminal ##
no_anoxy_geo_corr <- x_corr[c(-1:-3,-5),]
non_terminal_branches_df <- no_anoxy_geo_corr

for (taxon in geo_tips) {

    node <- match(taxon, anoxy_geo_time$tip)
    print(taxon)
    print(node)

    non_terminal_branches_df <- non_terminal_branches_df[-which(non_terminal_branches_df$E2 == node),]
}

non_terminal_only_pearson <- round(cor(non_terminal_branches_df[,3:4], method = "pearson"), digits = 4)[2]
non_terminal_only_spearman <- round(cor(non_terminal_branches_df[,3:4], method = "spearman"), digits = 4)[2]
writeLines(paste0("Pearson correlation of transfers over non-terminal Geobacillus branches only: ", non_terminal_only_pearson, "\nSpearman correlation of transers over non-terminal Geobacillus branches only: ", non_terminal_only_spearman))

###########################################################################
## Plot the species-derived tree ##
# devSVG(file = "/users/aesin/Dropbox/LSR/Presentation/Recipient_figures/Clado_T5_min0_2.svg")

plotBranchbyTrait(anoxy_geo_species, anoxy_geo_edge_df$mean, method="edges", legend = 5, title = "Fraction of transfers")

for (i in 1:nrow(anoxy_geo_edge_df)) {
    entry <- anoxy_geo_edge_df[i,]

    annotate_transfer_num(entry, i)
}

# dev.off()

####################################
## Plot on the time-resolved tree ##

# devSVG(file = "/users/aesin/Dropbox/LSR/Presentation/Recipient_figures/Time_T5_min0.svg")

plotBranchbyTrait(anoxy_geo_time, anoxy_geo_time_edge_df$mean, method="edges", legend = .1, title = "Fraction of transfers")

for (i in 1:nrow(anoxy_geo_time_edge_df)) {
    entry <- anoxy_geo_time_edge_df[i,]

    annotate_transfer_num(entry, i)
}

## Colour in the "missing" incongruent edges ##

for (i in length(unmapped_edges)) {
    edge <- unmapped_edges[[i]]
    ape::edges(edge[2], edge[1], col = "red", lty = 2, lwd = 2)
}

total_transfers <- nrow(receptor_df)
#legend("bottomright", legend = paste0("Min. sequences = ", min_taxa, ". At root: N = ", at_root, "\nBranch length vs HGT correlation: Overall = ", overall_pearson, "(P) & ", overall_spearman, "(S). Geo only = ", geo_only_pearson, "(P) & ", geo_only_spearman, "(S)."))

#dev.off()

####################################
## Simple correlation ##

# test <- as.data.frame(cbind(refined_x_corr[,1], refined_x_corr[,2]*100))
# names(test) <- c("Relative_transfer", "Edge_length")

# #devSVG(file = "/users/aesin/Dropbox/LSR/Presentation/Figures/Correlation/Correlation_t5_min50.svg", height = 15, width = 20)
# ggplot(test, aes(x = Edge_length, y = Relative_transfer)) + geom_point(color = "red", size = 5) + scale_y_continuous(limits = c(0,10), breaks = seq(from = 0, to = 10, by = 1)) + ggtitle("stuff") + theme(panel.grid.minor = element_blank())
# #dev.off()

###########################################################################
###########################################################################
## Prepare the plots with raw numbers of transfers per branch, at a given penalty ##

  # devSVG(file = "/users/aesin/Dropbox/LSR/Report/Figures/Figure_8/Clado_t5_min0.svg", height = 15, width = 20) #
      # devSVG(file = "/users/aesin/Dropbox/LSR/Report/Figures/Figure_8/Time_t5_min0.svg", height = 15, width = 20) #

if (single_pen_plot == TRUE) {

    colname_pen <- paste0("T", single_pen_val)
    single_pen_df <- anoxy_geo_edge_raw_df[,c(1:2, which(colnames(anoxy_geo_edge_raw_df) == colname_pen))]
    single_pen_time_df <- anoxy_geo_time_edge_raw_df[,c(1:2, which(colnames(anoxy_geo_time_edge_raw_df) == colname_pen))]

    ## Clado ##
    plotBranchbyTrait(anoxy_geo_species, single_pen_df[,3], method="edges", legend = 10, title = "Number of transfers")

    for (i in 1:nrow(single_pen_df)) {
        entry <- single_pen_df[i,]

        annotate_transfer_num(entry, i)
    }

    legend("bottomright", legend = paste0("Total transfers = ", sum(single_pen_df[,3]), ". Min. sequences = ", min_taxa))
    ## Time ##
    plotBranchbyTrait(anoxy_geo_time, single_pen_time_df[,3], method="edges", legend = .1, title = "Number of transfers")
    for (i in 1:nrow(single_pen_time_df)) {
        entry <- single_pen_time_df[i,]

        annotate_transfer_num(entry, i)
    }
    legend("bottomright", legend = paste0("Total transfers = ", sum(single_pen_time_df[,3]), ". Min. sequences = ", min_taxa))
    
}

###########################################################################
###########################################################################
## Prepare a correlation with differentially colored points depending on whether they are terminal / non-terminal. Label each point with an index that can be compared to a phylogeny to identify which point corresponds to which branch ##

x1<-cbind(single_pen_time_df, anoxy_geo_time$edge.length)
colnames(x1)[4] <- "branch_length"
geo_only_for_plot <- x1[c(-1:-3,-5,-15),]
no_outlier <- x1[c(-1:-3,-5, -7),]
devSVG(file = "/users/aesin/Dropbox/T5_min0_geo_only_correlation.svg", height = 15, width = 20) 
ggplot(geo_only_for_plot, aes(x = branch_length, y = T5)) + geom_point(size = 4, color = "red") + geom_smooth(method = "lm", se = F, formula = y~x-1)
spearman_corr_here <- round(cor(geo_only_for_plot[,3:4], method = "spearman"), digits = 4)[2]
pearson_corr_here <- round(cor(geo_only_for_plot[,3:4], method = "pearson"), digits = 4)[2]
dev.off()
#and - outlier #

####################################
## Phylogeny with the branches indexed ##
plotBranchbyTrait(anoxy_geo_time, anoxy_geo_time_edge_df$mean, method="edges", legend = .1, title = "Fraction of transfers")

for (i in 1:nrow(anoxy_geo_time_edge_df)) {
    # entry <- anoxy_geo_time_edge_df[i,]

    # edge <- edge_entry[1:2]
    # print(edge_label)

    edgelabels(i, i, adj = c(0.5, -0.25), bg = "white", frame = "none", cex = 0.7)
}

####################################
## Correlation with differentially labelled points ##

geo_only_with_index <- cbind(geo_only_with_edge, rownames(geo_only))
names(geo_only_with_index)[3:5] <- c("mean_fraction_transfers", "branch_length", "point_index")
geo_only_with_index$terminal <- 0

for (i in 1:nrow(geo_only_with_index)) {
    edge_1 <- geo_only_with_index[i,1]
    edge_2 <- geo_only_with_index[i,2]

    x <- terminal_branches_df[which(terminal_branches_df$E1 == edge_1 & terminal_branches_df$E2 == edge_2),]$E1
    if (length(x) == 1) {
        geo_only_with_index[i,6] <- "Terminal"
    } else {
        geo_only_with_index[i,6] <- "Not_terminal"
    }

}

ggplot(geo_only_with_index, aes(x = branch_length, y = mean_fraction_transfers, label = point_index, color = terminal)) + geom_point(size = 3) + geom_text(aes(label = point_index), hjust = 0, vjust = -1)


###########################################################################
###########################################################################

test_1 <- cbind(anoxy_geo_time_edge_raw_df, anoxy_geo_time$edge.length)[c(-1:-3,-5,-15),]
test_2 <- test_1[,-1:-2]
colnames(test_2)[4] <- "branch_length"
test_3 <- merge(geo_only_with_index, test_2, by.x = "branch_length")

ggplot(test_3, aes(x = branch_length, y = T5, label = point_index, color = terminal)) + geom_point(size = 4)