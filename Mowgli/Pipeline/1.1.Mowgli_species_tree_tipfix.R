#!/usr/bin/Rscript
## This script changes any tip labels ending in "_[number]" to "_[number]A" so that the Mowgli program doesn't interpret the ending numbers as reference for taxon multiples ##
###########################################################################
## Packages ##
library(ape)
library(phangorn)
library(geiger)
library(stringr)

## This function forces the write.tree command to write out branches with decimal places rather than scientific notation - the latter occasionally causes errors (e.g. in Ranger-DTL) ##
.write.tree2 <- function (phy, digits = 10, tree.prefix = "") {
  brl <- !is.null(phy$edge.length)
  nodelab <- !is.null(phy$node.label)
  phy$tip.label <- checkLabel(phy$tip.label)
  if (nodelab) 
    phy$node.label <- checkLabel(phy$node.label)
  f.d <- paste("%.", digits, "f", sep = "")
  cp <- function(x) {
    STRING[k] <<- x
    k <<- k + 1
  }
  add.internal <- function(i) {
    cp("(")
    desc <- kids[[i]]
    for (j in desc) {
      if (j > n) 
        add.internal(j)
      else add.terminal(ind[j])
      if (j != desc[length(desc)]) 
        cp(",")
    }
    cp(")")
    if (nodelab && i > n) 
      cp(phy$node.label[i - n])
    if (brl) {
      cp(":")
      cp(sprintf(f.d, phy$edge.length[ind[i]]))
    }
  }
  add.terminal <- function(i) {
    cp(phy$tip.label[phy$edge[i, 2]])
    if (brl) {
      cp(":")
      cp(sprintf(f.d, phy$edge.length[i]))
    }
  }
  n <- length(phy$tip.label)
  parent <- phy$edge[, 1]
  children <- phy$edge[, 2]
  kids <- vector("list", n + phy$Nnode)
  for (i in 1:length(parent)) kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])
  ind <- match(1:max(phy$edge), phy$edge[, 2])
  LS <- 4 * n + 5
  if (brl) 
    LS <- LS + 4 * n
  if (nodelab) 
    LS <- LS + n
  STRING <- character(LS)
  k <- 1
  cp(tree.prefix)
  cp("(")
  getRoot <- function(phy) phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]
  root <- getRoot(phy)
  desc <- kids[[root]]
  for (j in desc) {
    if (j > n) 
      add.internal(j)
    else add.terminal(ind[j])
    if (j != desc[length(desc)]) 
      cp(",")
  }
  if (is.null(phy$root.edge)) {
    cp(")")
    if (nodelab) 
      cp(phy$node.label[1])
    cp(";")
  }
  else {
    cp(")")
    if (nodelab) 
      cp(phy$node.label[1])
    cp(":")
    cp(sprintf(f.d, phy$root.edge))
    cp(";")
  }
  paste(STRING, collapse = "")
}
###########################################################################

setwd("/Users/aesin/Desktop/Mowgli/Species_tree")
species_tree <- read.tree("Ultrametric_species_tree.txt")

# tip_labels <- species_tree$tip.label
output_str <- ""
for (n in 1:length(species_tree$tip.label)) {
	label_name <- species_tree$tip.label[n]
	
	new_label <- str_replace(label_name, "(_[0-9]+$)", "\\1A")

	if (identical(label_name, new_label) == FALSE) {
		species_tree$tip.label[n] <- new_label
		output_str <- str_c(output_str, "\n", label_name, "\t", new_label)
	} else {
		output_str <- str_c(output_str, "\n", label_name, "\t", label_name)
	}
}

## The digits argument is given so that the edge lengths aren't written out as "1.0000000000" -- purely aesthetic ##
write.tree (species_tree, file = "Ultrametric_species_tree_tipfix.txt", digits = 1)

output_str <- str_trim(output_str, side = "both")
write(output_str, file = "Tip_fix_key.tsv")