library(ape)
library(phylobase)
library(phytools)
library(stringr)
library(geiger)
library(RSvgDevice)
library(reshape2)
library(ggplot2)

################################
## Annotate branches function ##
################################

annotate_transfer_num <- function (edge_entry, index) {
    edge <- edge_entry[1:2]
    edge_label <- as.character(round(edge_entry[3], digits = 2))
    # print(edge_label)

    edgelabels(edge_label, index, adj = c(0.5, -0.25), bg = "white", frame = "none", cex = 0.7)
}

###############
## Variables ##
###############

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

##############################
# Prepare the cladogram tree #
##############################

species_tree_dir = "/Users/aesin/Desktop/Mowgli/Species_tree/Reconciled_species_tree"
setwd(species_tree_dir)

species_tree_clado <- read.tree("outputSpeciesTree.mpr")

## Get the set of Geobacillus and Anoxybacillus tips and combine into a single set ##
geo_tips <- grep("Geobacillus_", species_tree_clado$tip.label, value = TRUE)
anoxy_tips <- grep("Anoxybacillus_", species_tree_clado$tip.label, value = TRUE)
anoxy_geo_clado_tips <- c(geo_tips, anoxy_tips)

## Use the Anoxy/Geobacillus tip set to extract the Anoxy/Geobacillus sublclade. Transform to phylo4 object ##
anoxy_geo_clado <- drop.tip(species_tree_clado, setdiff(species_tree_clado$tip.label, anoxy_geo_species_tips))
anoxy_geo_clado4 <- phylo4(anoxy_geo_clado)

## Get the node labels as a data frame ##
node_lbs_df <- data.frame(labels(anoxy_geo_clado4, "all"))
internal_nodes_species <- nodeId(anoxy_geo_clado4, "internal")
names(node_lbs_df)[1] <- "Labels"

## Prepare a data frame with the edges so we can count ##
anoxy_geo_edge_df <- data.frame(anoxy_geo_clado$edge)
names(anoxy_geo_edge_df) <- c("E1", "E2")

## Make a column for each penalty ##
for (penalty in penalties) {
    cname <- paste0("T", penalty)
    anoxy_geo_edge_df$colblank <- 0
    names(anoxy_geo_edge_df)[names(anoxy_geo_edge_df)=="colblank"] <- cname
}

##################################
# Prepare the time-resolved tree #
##################################

## Read in the time-resolved tree ##
bacillaceae_tree <- read.tree(file = "/Users/aesin/Desktop/Consensus_trees/Consensus_finals/Bacillaceae_N1000_newick.txt")

geo_tips <- grep("Geobacillus_", bacillaceae_tree$tip.label, value = TRUE)
anoxy_tips <- grep("Anoxybacillus_", bacillaceae_tree$tip.label, value = TRUE)
anoxy_geo_time_tips <- c(geo_tips, anoxy_tips)

anoxy_geo_time <- drop.tip(bacillaceae_tree, setdiff(bacillaceae_tree$tip.label, anoxy_geo_time_tips))

## Annotate the time-resolved tree with the correct node labels imported from the species tree-derived tree ##
anoxy_geo_time_relab <- anoxy_geo_time
for (tip in anoxy_geo_time_tips) {
    species_tree_tip_name <- grep(tip, anoxy_geo_species_tips, value = T)
    original_tip_position <- match(tip, anoxy_geo_time_relab$tip.label)
    # Replace with new name #
    anoxy_geo_time_relab$tip.label[original_tip_position] <- species_tree_tip_name
}

## Remove node labels (BS-value) and convert to phylo4 ##
anoxy_geo_time_relab$node.label <- NULL
anoxy_geo_time_relab_4 <- phylo4(anoxy_geo_time_relab)

## Get all the internal nodes
internal_nodes_time <- nodeId(anoxy_geo_time_relab_4, "internal")
unmapped_edges <- vector(mode = "list")
unmapped_count = 1

## For each internal node, match it to the nodes from the species tree ##
for (internal_time in internal_nodes_time) {
    matched = FALSE
    for (internal_species in internal_nodes_species) {
        if (setequal(tips(anoxy_geo_time_relab, internal_time), tips(anoxy_geo_clado, internal_species)) == TRUE) {

            ## Get the label of the node from the species-cladogram tree ##
            node_label <- anoxy_geo_clado4@label[internal_species]

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


# Prepare the count table for the time-resolved tree as for the species-derived prune above ##
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


#####################################################################
## Read in the number of transfers at each penalty for each branch ##
#####################################################################

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


##########################################################################
## Process a copy of the above df to be an average across the penalties ##
##########################################################################

## Save a copy of each df with the raw numbers ##
anoxy_geo_edge_raw_df <- anoxy_geo_edge_df
anoxy_geo_time_edge_raw_df <- anoxy_geo_time_edge_df

## Make the transfers in at each branch a fraction of the total transfers ##
anoxy_geo_edge_df <- cbind(anoxy_geo_edge_df[,1:2],sweep(anoxy_geo_edge_df[,-1:-2],2,colSums(anoxy_geo_edge_df[,-1:-2]), "/")*100)
anoxy_geo_time_edge_df <- cbind(anoxy_geo_time_edge_df[,1:2],sweep(anoxy_geo_time_edge_df[,-1:-2],2,colSums(anoxy_geo_time_edge_df[,-1:-2]), "/")*100)

## Add a mean column ##
anoxy_geo_edge_df$mean <- rowMeans(anoxy_geo_edge_df[,-1:-2])
anoxy_geo_time_edge_df$mean <- rowMeans(anoxy_geo_time_edge_df[,-1:-2])


#############################################################
## Process the correlation and graphs for a single penalty ##
#############################################################

## Restrict the transfer number to a specific penalty, and combine it with the branch length data from the time-resolved tree ##
colname_pen <- paste0("T", single_pen_val)
anoxy_geo_brnch_len <- cbind(anoxy_geo_time_edge_raw_df[,c(1:2, which(colnames(anoxy_geo_time_edge_raw_df) == colname_pen))], anoxy_geo_time$edge.length)

## Restrict the data to only Geobacillus branches, including 
geo_only_brnch_len <- anoxy_geo_brnch_len[c(-1:-5),]
geo_only_correlation <- geo_only_brnch_len[,3:4]

geo_only_pearson <- round(cor(geo_only_correlation, method = "pearson"), digits = 4)[2]
geo_only_spearman <- round(cor(geo_only_correlation, method = "spearman"), digits = 4)[2]

writeLines(paste0("Pearson correlation of transfers over all Geobacillus branches: ", overall_pearson, "\n", "Spearman correlation of transfers over all Geobacillus branches: ", overall_spearman))


###################################################################
## Plot the transfer data onto the clado and time-resolved trees ##
###################################################################

## Plot the species-derived (clado) tree using the mean fraction of transfers across penalties ##

# devSVG(file = "/users/aesin/Dropbox/LSR/Presentation/Recipient_figures/Clado_T5_min0_2.svg")
plotBranchbyTrait(anoxy_geo_clado, anoxy_geo_edge_df$mean, method="edges", legend = 5, title = "Fraction of transfers")

for (i in 1:nrow(anoxy_geo_edge_df)) {
    entry <- anoxy_geo_edge_df[i,]

    annotate_transfer_num(entry, i)
}
# dev.off()

## Plot the time-resolved tree using the mean fraction of transfers across penalties. NB one incongruent edge ##

# devSVG(file = "/users/aesin/Dropbox/LSR/Presentation/Recipient_figures/Time_T5_min0.svg")
plotBranchbyTrait(anoxy_geo_time, anoxy_geo_time_edge_df$mean, method="edges", legend = .1, title = "Fraction of transfers")

for (i in 1:nrow(anoxy_geo_time_edge_df)) {
    entry <- anoxy_geo_time_edge_df[i,]

    annotate_transfer_num(entry, i)
}

# Colour in the "missing" incongruent edges #
for (i in length(unmapped_edges)) {
    edge <- unmapped_edges[[i]]
    ape::edges(edge[2], edge[1], col = "red", lty = 2, lwd = 2)
}

total_transfers <- nrow(receptor_df)
legend("bottomright", legend = paste0("Min. sequences = ", min_taxa, ". At root: N = ", at_root, "\nBranch length vs HGT correlation: Overall = ", overall_pearson, "(P) & ", overall_spearman, "(S). Geo only = ", geo_only_pearson, "(P) & ", geo_only_spearman, "(S)."))
#dev.off()


####################################################################################
## Prepare the plots with raw numbers of transfers per branch, at a given penalty ##
####################################################################################

# devSVG(file = "/users/aesin/Dropbox/LSR/Report/Figures/Figure_8/Clado_t5_min0.svg", height = 15, width = 20) #
# devSVG(file = "/users/aesin/Dropbox/LSR/Report/Figures/Figure_8/Time_t5_min0.svg", height = 15, width = 20) #

if (single_pen_plot == TRUE) {

    colname_pen <- paste0("T", single_pen_val)
    single_pen_df <- anoxy_geo_edge_raw_df[,c(1:2, which(colnames(anoxy_geo_edge_raw_df) == colname_pen))]
    single_pen_time_df <- anoxy_geo_time_edge_raw_df[,c(1:2, which(colnames(anoxy_geo_time_edge_raw_df) == colname_pen))]

    ## Clado ##
    plotBranchbyTrait(anoxy_geo_clado, single_pen_df[,3], method="edges", legend = 10, title = "Number of transfers")

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


#########################################
## Phylogeny with the branches indexed ##
#########################################

plotBranchbyTrait(anoxy_geo_time, anoxy_geo_time_edge_df$mean, method="edges", legend = .1, title = "Fraction of transfers")

for (i in 1:nrow(anoxy_geo_time_edge_df)) {
    # entry <- anoxy_geo_time_edge_df[i,]

    # edge <- edge_entry[1:2]
    # print(edge_label)

    edgelabels(i, i, adj = c(0.5, -0.25), bg = "white", frame = "none", cex = 0.7)
}

##########################
## Subspecies groupings ##
##########################

## Aliyu et al 2016 divides the Geobacillus phylogeny into a number of species subgroups (fig 3). I applied that classification to the genomes used in this analysis. Branches separating species/tips belonging to the same subgroup are considered to not have undergone substantial differentiation and purifying selection. ##

## Read in the table constructed based on the publication ##
subspecies_groupings <- read.csv(file = "/users/aesin/Desktop/Geo_analysis/Geo_omes/Subspecies_geo_groups.csv", header = TRUE, sep = ",")

## Get those tips that have a fellow tip in the same subgroup. Isolate all the subgroups that have more than one member ##
subspecies_only <- subset(subspecies_groupings, duplicated(subspecies_groupings$Group) | duplicated(subspecies_groupings$Group, fromLast = TRUE))
unique_subgroups <- unique(subspecies_only$Group)

## Import the table containing anoxy/geo edges, associated branch lengths, and the num. transfers in for a given penalty. Create new column for whether a branch is considered to be subspecies or inter-species ##
sinpen_time_sub_df <- anoxy_geo_brnch_len
names(sinpen_time_sub_df)[ncol(sinpen_time_sub_df)] <- "branch_length"
sinpen_time_sub_df[,"Subspecies_branch"] <- NA

## For each group with more than one member, identify the species. Find all the descendant nodes from the common ancestor of the subspecies (they are monophyletic). Using those nodes, identify all the edges (branches) that we will classify as being a subspecies branch (1). All other branches are interspecies (0) ##
for (group in unique_subgroups) {
    tips_in_subgroup <- as.vector(subspecies_only$Species[which(subspecies_only$Group == group)])
    descend_nodes <- as.vector(descendants(t1, MRCA(t1, tips_in_subgroup), "all"))

    for (node in descend_nodes) {
        print(node)
        sinpen_time_sub_df$Subspecies_branch[which(sinpen_time_sub_df$E2 == node)] <- 1
    }
}

sinpen_time_sub_df$Subspecies_branch[is.na(sinpen_time_sub_df$Subspecies_branch)] <- 0

## Add a point index column, corresponding to the index used in the tree - so that each point in the num. transfer vs branch length can be linked to a specific branch ##
sinpen_time_sub_index_df <- cbind(sinpen_time_sub_df, rownames(sinpen_time_sub_df))
names(sinpen_time_sub_index_df)[ncol(sinpen_time_sub_index_df)] <- "point_index"


## For the correlation plot, restrict the data to only Geobacillus branches. Convert 0/1 to FALSE/TRUE so that the plot is in two colors (discrete). We also need to remove one branch point - because it's incongruent between the species tree used for the reconciliations and the time-resolved tree from which we get the branch lengths ##
geo_only_full_df <- sinpen_time_sub_index_df[c(-1:-5,-15),]
geo_only_full_logical <- geo_only_full_df
geo_only_full_logical$Subspecies_branch <- as.logical(geo_only_full_logical$Subspecies_branch)

## Plot ##
ggplot(geo_only_full_logical, aes(x = branch_length, y = T5, label = point_index, color = Subspecies_branch)) + geom_point(size = 3) + geom_text(aes(label = point_index), hjust = 0, vjust = -1)


## Plot with regression lines, but exclude the outlier (point-index = 7) ##
geo_only_full_logical_no_out_df <- geo_only_full_logical[-which(geo_only_full_logical$point_index == 7),]
ggplot(geo_only_full_logical_no_out_df, aes(x = branch_length, y = T5, label = point_index, color = Subspecies_branch)) + geom_point(size = 3) + geom_text(aes(label = point_index), hjust = 0, vjust = -1) + geom_smooth(method = "lm", se = FALSE)


## Plot the tree highlighted based on subspecies/inter-species branches. Label the edges based on point index ##
plotBranchbyTrait(anoxy_geo_time, sinpen_time_sub_index_df$Subspecies_branch, method="edges", title = "Subspecies branches")
for (i in 1:nrow(anoxy_geo_time_edge_df)) {
    edgelabels(i, i, adj = c(0.5, -0.25), bg = "white", frame = "none", cex = 0.7)
}

##################
## Correlations ##
##################

## Correlate all Geobacillus branches ##
cor_all <- geo_only_full_logical[,3:4]
cor(cor_not_subsp, method = "pearson")

# Pearson : 0.3313325
# Spearman : 0.5943488


## Just those branches which are not "sub-species" ## 
cor_not_subsp <- geo_only_full_logical[which(geo_only_full_logical$Subspecies_branch == "FALSE"),][,c(3,4)]
cor(cor_not_subsp, method = "spearman")

# Pearson : 0.3038134
# Spearman : 0.5633697

## Outlier corrected - branch 7 ##
cor_not_subsp_outlier <- cor_not_subsp[-2,]
cor(cor_not_subsp_outlier, method = "spearman")

# Pearson : 0.7872012
# Spearman : 0.7766565

## Subspecies only branches ##
cor_subsp <- geo_only_full_logical[which(geo_only_full_logical$Subspecies_branch == "TRUE"),][,c(3,4)]
cor(cor_subsp, method = "spearman")

# Pearson : 0.4032259
# Spearman : 0.4341429



########################
## Simple correlation ##
########################

# ggplot(test, aes(x = Edge_length, y = Relative_transfer)) + geom_point(color = "red", size = 5) + scale_y_continuous(limits = c(0,10), breaks = seq(from = 0, to = 10, by = 1)) + ggtitle("stuff") + theme(panel.grid.minor = element_blank())

##############################
### OLD CORRELATIONS BELOW ###
##############################



## Perform correlations ##
# x_corr <- cbind(anoxy_geo_time_edge_df[,1:2],anoxy_geo_time_edge_df$mean, anoxy_geo_time$edge.length)
# overall_pearson <- round(cor(x_corr[,3:4], method = "pearson"), digits = 4)[2]
# overall_spearman <- round(cor(x_corr[,3:4], method = "spearman"), digits = 4)[2]
# writeLines(paste0("Pearson correlation of transfers over all branch lengths: ", overall_pearson, "\n", "Spearman correlation of transfers over all branch lengths: ", overall_spearman))


# ## Remove those taxa that are based on incomplete chromosome assembly ##
# refined_x_corr <- x_corr
# poor_assembly_taxa <- c("Geobacillus_caldoxylosilyticus", "Geobacillus_sp_G11MC16", "Geobacillus_sp_G1w1", "Geobacillus_thermocatenulatus", "Geobacillus_sp_WSUCF1", "Geobacillus_sp_MAS1")

# for (taxon in poor_assembly_taxa) {
#     node <- match(taxon, anoxy_geo_time$tip)
#     print(node)
#     refined_x_corr <- refined_x_corr[-which(refined_x_corr$E2 == node),]
# }

# contig_removed_pearson <- round(cor(refined_x_corr[,3:4], method = "pearson"), digits = 4)[2]
# contig_removed_spearman <- round(cor(refined_x_corr[,3:4], method = "spearman"), digits = 4)[2]

# writeLines(paste0("Pearson correlation of transfers over branches with contig-derived assembly branches removed: ", contig_removed_pearson, "\nSpearman correlation of transfers over branches with contig-derived assembly branches removed: ", contig_removed_spearman))

# ## Remove anoxybacillus taxa and supporting nodes ##
# geo_only_with_edge <- x_corr[c(-1:-3,-5),]
# geo_only <- geo_only_with_edge[,3:4]

# geo_only_pearson <- round(cor(geo_only, method = "pearson"), digits = 4)[2]
# geo_only_spearman <- round(cor(geo_only, method = "spearman"), digits = 4)[2]
# writeLines(paste0("Pearson correlation of transfers over Geobacillus only branches: ", geo_only_pearson, "\nSpearman correlation of transers over Geobacillus only branches: ", geo_only_spearman))


# ## Keep only the branches that are terminal ##
# terminal_branches_df <- data.frame(E1 = numeric(), E2 = numeric(), count = numeric(), edge_length = numeric())

# no_anoxy_geo_corr <- x_corr[c(-1:-3,-5),]

# for (taxon in geo_tips) {

#     node <- match(taxon, anoxy_geo_time$tip)
#     print(taxon)
#     print(node)

#     terminal_branches_df[nrow(terminal_branches_df)+1,] <- no_anoxy_geo_corr[which(no_anoxy_geo_corr$E2 == node),]
# }

# terminal_only_pearson <- round(cor(terminal_branches_df[,3:4], method = "pearson"), digits = 4)[2]
# terminal_only_spearman <- round(cor(terminal_branches_df[,3:4], method = "spearman"), digits = 4)[2]

# writeLines(paste0("Pearson correlation of transfers over terminal Geobacillus branches only: ", terminal_only_pearson, "\nSpearman correlation of transers over terminal Geobacillus branches only: ", terminal_only_spearman))

# ## Keep only branches that are NON-terminal ##
# no_anoxy_geo_corr <- x_corr[c(-1:-3,-5),]
# non_terminal_branches_df <- no_anoxy_geo_corr

# for (taxon in geo_tips) {

#     node <- match(taxon, anoxy_geo_time$tip)
#     print(taxon)
#     print(node)

#     non_terminal_branches_df <- non_terminal_branches_df[-which(non_terminal_branches_df$E2 == node),]
# }

# non_terminal_only_pearson <- round(cor(non_terminal_branches_df[,3:4], method = "pearson"), digits = 4)[2]
# non_terminal_only_spearman <- round(cor(non_terminal_branches_df[,3:4], method = "spearman"), digits = 4)[2]
# writeLines(paste0("Pearson correlation of transfers over non-terminal Geobacillus branches only: ", non_terminal_only_pearson, "\nSpearman correlation of transers over non-terminal Geobacillus branches only: ", non_terminal_only_spearman))