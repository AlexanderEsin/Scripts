#!/usr/bin/Rscript
## Relabel gene tree tips to be used by mowgli. 

##The raxml run produces a tree with tip labels as follows:
## Binomial_name {Class} {Protein_ID.species_ID}
## R automatically collapses the white spaces in the tip labels, so "Binomial_name{Class}{Protein_ID.species_ID}"

## This script removes the class/protein_ID terms, returning just the binomial Binomial_name
## If this binomial name ends in "..._[number]" mowgli would misinterpret this as labelling duplicate taxa, so we replace it with:
## "Binomial_name..._[number]" --> "Binomial_name..._[number]A"

## THEN, for taxa that do appear multiple time, we want to label them as repeats - so we then ADD a "_[number]" onto the end of the tip name

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
##### SWITCHES #####

## Will this be input for treefix/ranger or mowgli? For mowgli the trees should be midpoint rooted ##
mowgli_mp_root_switch = 1; ## 1 = mp root / 0 = unrooted

## Should a map file be made in the Treefix-DTL format? ##
make_map_switch = 1; # 1 = ON / 0 = OFF


###########################################################################

directory = "/users/aesin/Desktop/Mowgli"
final_bs_dir_name = "Final_BS"
relabelled_dir_name = "Tips_relabelled_BS"
relabelled_nbs_dir_name = "Tips_relabelled_NBS"
key_dir_name = "Tip_keys"

final_bs_dir_path = paste0(directory, "/", final_bs_dir_name)
relabelled_dir_path = paste0(directory, "/", relabelled_dir_name)
relabelled_nbs_dir_path =paste0(directory, "/", relabelled_nbs_dir_name)
key_dir_path = paste0(directory, "/", key_dir_name)

dir.create(relabelled_dir_path, showWarnings = FALSE)
dir.create(key_dir_path, showWarnings = FALSE)
dir.create(relabelled_nbs_dir_path, showWarnings = FALSE)

if (make_map_switch == 1) {
	map_dir_name = "Species_maps"
	map_dir_path = paste0(directory, "/", map_dir_name)
	dir.create(map_dir_path, showWarnings = FALSE)
}

setwd(final_bs_dir_path)
input_trees <- Sys.glob(paste0(getwd(), "/*txt"))

counter = 1

for (tree in input_trees) {

	name = str_extract(tree, "[0-9]+")
	if (mowgli_mp_root_switch == 1) {
		output_tree_bs_name = paste0(relabelled_dir_path, "/Relabelled_mp_", name, ".txt")
		output_tree_nbs_name = paste0(relabelled_nbs_dir_path, "/Relabelled_mp_", name, ".txt")
	} else {
		output_tree_name = paste0(relabelled_dir_path, "/Relabelled_unr_", name, ".txt")
	}
	

	output_key_name = paste0(key_dir_path, "/TIP_KEY_", name, ".tsv")
	output_key = ""

	if (make_map_switch == 1) {
		output_map_name = paste0(map_dir_path, "/SMAP_", name, ".smap")
		output_map = ""
	}

	## Read and rewrite to correct labels containing double underscores ##
	tree_data <- read.tree(tree)
	write.tree(tree_data, file = tree)
	tree_data <- read.tree(tree)

	test_vector <- ""

	for (n in 1:length(tree_data$tip.label)) {

		## Remove the bracketed terms ##
		old_label <- tree_data$tip.label[n]
		new_label <- gsub("\\{.+\\}", "", old_label)

		## Attempt to fix any trailing "_[number]" labels ##
		num_fix_label <- str_replace(new_label, "(_[0-9]+$)", "\\1A")

		## Check whether the num_fix_label is the same as the original trimmed label - i.e. was there a trailing number that was modified? If so, make the fixed label the new output ##
		if (identical(new_label, num_fix_label) == FALSE) {
			new_label <- num_fix_label
		}
		test_vector <- c(test_vector, new_label)

		## Check whether the binomial label appears more than once, if so - number the binomial name accordingly. Otherwise, label the tip as _1 (otherwise Mowgli gets upset) ##
		multi_taxa <- length(grep(num_fix_label, test_vector))
		new_label <- paste0(num_fix_label, "_", multi_taxa)
			

		## Rename the tip in the tree object ##
		tree_data$tip.label[n] <- new_label

		## Append the output_key string with the relevant translation ##
		output_key <- str_c(output_key, "\n", old_label, "\t", new_label)

		if (make_map_switch == 1) {
			## Add the entry to the map file IF it is NOT already there ##
			if (str_detect(output_map, num_fix_label) == FALSE) {
				output_map <- str_c(output_map, "\n", paste0(num_fix_label, "_*", "\t", num_fix_label))
			}
		}
		
	}
	if (mowgli_mp_root_switch == 1) {
		tree_data <- midpoint(tree_data)
	}

	if (exists('node.label', where = tree_data) == TRUE) {
			write.tree(tree_data, file = output_tree_bs_name)
	} else {
			write.tree(tree_data, file = output_tree_nbs_name)
	}

	## Trim and write out the key information ##
	output_key <- str_trim(output_key, side = "both")
	write(output_key, file = output_key_name)

	if (make_map_switch == 1) {
		## Trim and write out the map information ##
		output_map <- str_trim(output_map, side = "both")
		write(output_map, file = output_map_name)
	}
	
	print(counter)
	counter = counter + 1
	
}