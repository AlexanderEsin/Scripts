### Identify the members of the Inside_group taxa which are outliers, i.e. did not cluster with the rest of the taxa during the initial Astral run with the gene trees.
### Produce a pruned tree missing those taxa so that they can be readded in a subsequent Astral run using the the 1:1 orthologous gene trees from the Inside group to place the taxa.


library(ape)
library(geiger)
library(phangorn)

working_dir <- "/users/aesin/desktop/Geo_final_trees"
reference_dir <- paste0(working_dir, "/Reference")
astral_output_dir <- paste0(working_dir, "/Astral_output")

## Open the Inside Group reference tree and extract the tip set ##

setwd(reference_dir)
inside_reference_tree <- "Inside_group_N100.txt"
inside_reference_data <- read.tree(inside_reference_tree)
inside_reference_tips <- inside_reference_data$tip.label

## Open the Astral inital output and read in the tree ##

setwd(astral_output_dir)
astral_output_tree <- "Astral_output_all.txt"
astral_output_data <- read.tree(astral_output_tree)

## Root the tree based on an arbitary distantly related group ##

astral_output_data_rooted <- root(astral_output_data, "Homo_sapiens")

## Assign arbitrary edge lengths ##

astral_output_data_rooted$edge.length <- rep_len(1, (length(astral_output_data$edge)/2)-1)

## Build matrix of distances between all tips in tree ##

dist_matrix_full <- cophenetic.phylo(astral_output_data_rooted)

## Extract those for the shared tips only ##

dist_matrix_sub <- dist_matrix_full[inside_reference_tips, inside_reference_tips]

## Hierarchically cluster the distances, and isolate the two most basal groups ##

d <- dist(dist_matrix_sub, method = "euclidean")
fit <- hclust(d, method = "ward.D2")
groups <- cutree(fit, k = 2)

## Optionally plot the cluster dendrogram ##

par(cex=0.3, mar=c(5, 8, 4, 1))
plot(fit, xlab="", ylab="", main="", sub="", axes=FALSE)
par(cex=1)
title(xlab="Taxa", ylab="Relative distance", main="Cluster Dendrogram")
rect.hclust(fit, k=2, border="red")

## Identify which group is the smaller, they would be the outliers ##

if (sum(length(which(groups == 1))) > sum(length(which(groups == 2)))) {
	outliers <- c(names(which(groups == 2)))
} else {
	outliers <- c(names(which(groups == 1)))
}

## Prune the outliers from the group ## 

pruned_astral_data <- drop.tip(astral_output_data, outliers)

## Write the pruned tree out ##

write.tree(pruned_astral_data, file = "Inside_group_outliers_pruned.txt")