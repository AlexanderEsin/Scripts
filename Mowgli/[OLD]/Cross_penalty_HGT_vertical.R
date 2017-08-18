#!/usr/bin/Rscript
## Compare the results values at different transfer penalty costs and extract those groups which are positive for HGT signal into Anoxy/Geo at all the penalties ##
library(stringr)
library(plyr)
library(dplyr)
library(reshape2)
library(lazyeval)
library(gtools)

MyMerge <- function(x, y){
	df <- merge(x, y, by = "V1")
}

###########################################################################
## Set whether the mowgli data parser excluded all those trees that contained IG tips only ##
IG_inclusion <- TRUE
## Set whether we are testing all the scenarios or just 1 + 2 (those where the root is NOT in Anoxy/Geo) ##
all_scenarios <- FALSE; ## FALSE gives us ONLY scenarios 1 + 2

###########################################################################
## Get and order the results files by dictionary order (package = gtools) ##
if (IG_inclusion == TRUE) {
	in_path <- "/users/aesin/Desktop/Mowgli/Mowgli_outputs/Results_Full"
	dir.create("/users/aesin/Desktop/Mowgli/Mowgli_outputs/IG_intersect/All_scenarios", showWarnings = FALSE)
} else {
	in_path <- "/users/aesin/Desktop/Mowgli/Mowgli_outputs/Results_NO_IG"
	dir.create("/users/aesin/Desktop/Mowgli/Mowgli_outputs/NO_IG_intersect/All_scenarios", showWarnings = FALSE)
}

setwd(in_path)
results_files <- mixedsort(grep("results", dir(), value = TRUE))

## We want to do all the necessary permutations of cross-penalty HGT. E.g. for penalties 3, 4, 5 & 6 we want to do: 3 + 4, 3 + 4 + 5, 3 + 4 + 5 + 6 ##
results_counter = 2
while (results_counter <= length(results_files)) {

	## Subset the results files to do all the necessary permutations ##
	setwd(in_path)
	results_subset <- results_files[1:results_counter]

	## Prepare the penalty list to hold all the data and vectors for the penalty_group_num_df dataframe ##
	penalty_list <- vector("list", length(results_subset))
	penalty_values <- vector()
	group_nums <- vector()

	###########################################################################
	## Process each results file in turn ##
	# results_files = "T3_D2_L1_results.txt" "T4_D2_L1_results.txt" #

	i = 1
	for (file in results_subset) {
		transfer_penalty <- as.integer(str_sub(regmatches(file, regexpr("[[:alnum:]+]+", file)), 2))
		penalty_values <- c(penalty_values, transfer_penalty)

		## Assing data both to a new table (for individual overview) and into the penalty_list list for downstream merging ##
		name <- paste0("results_table_", transfer_penalty)
		read_table  <- read.table(file = file, header = FALSE, sep = "\t")
		assign(name, read_table)
		penalty_list[[i]] <- read_table
		colnames(penalty_list[[i]])[2] <- paste0("HGT_", transfer_penalty)
		colnames(penalty_list[[i]])[3] <- paste0("T", transfer_penalty, "_trans_type")
		i = i + 1

		## Get the total number of entries (i.e. groups processed by Mowgli) for each table ##
		group_nums <- c(group_nums, nrow(read_table))
	}

	###########################################################################
	## Print the values of the transfer penalties being compared ##
	message(cat("Transfer penalties being intersected:", penalty_values))

	###########################################################################
	## Create a data frame to show how many groups were reconciled to produce the HGT/Vertical splits at each penalty ##

	penalty_group_num_df <- as.data.frame(cbind(penalty_values, group_nums))
	colnames(penalty_group_num_df) <- c("HGT_Penalty", "Number_groups_queried")

	###########################################################################
	## Produce the table and sort it down to either those that all give positive HGT or positive vertical signal ##

	merged_results <- Reduce(MyMerge, penalty_list)

	# The set of conditions for the filter function are provided as a string #
	HGT_conditions <- paste(colnames(merged_results[,c(seq(2,ncol(merged_results), by = 2))]), "== 1", collapse = " & ")
	vertical_conditions <- paste(colnames(merged_results[,c(seq(2,ncol(merged_results), by = 2))]), "== 0", collapse = " & ")

	## If we want all the scenarios, then we do not need to filter further ##
	if (all_scenarios == TRUE) {
		# The condition string can be used with the standard evaluation version of the filter (filter_) function in dplyr #
		true_HGT_intersect <- filter_(merged_results, HGT_conditions)$V1
		true_vertical_intersect <- filter_(merged_results, vertical_conditions)$V1
	} else {
		## Here we need to furhter filter based on the scenario breakdown. Obviously the vertical doesn't change ##
		true_vertical_intersect <- filter_(merged_results, vertical_conditions)$V1

		true_HGT_intersect <- filter_(merged_results, HGT_conditions)
		## Second filtering step here ##
		scenario_filter <- paste(colnames(merged_results[,c(seq(3,ncol(merged_results), by = 2))]), "< 3", collapse = " & ")
		true_HGT_intersect <- filter_(true_HGT_intersect, scenario_filter)$V1
	}



	###########################################################################
	## Prepare string for the output name based on the penalties queried ##

	penalty_value_str <- paste(penalty_values, collapse = "_t")

	if (IG_inclusion == TRUE) {
		if (all_scenarios == TRUE) {
			setwd("/users/aesin/Desktop/Mowgli/Mowgli_outputs/IG_intersect/All_scenarios")
		} else {
			dir.create("/users/aesin/Desktop/Mowgli/Mowgli_outputs/IG_intersect/Scenarios_1_2", showWarnings = FALSE)
			setwd("/users/aesin/Desktop/Mowgli/Mowgli_outputs/IG_intersect/Scenarios_1_2")
		}	
	} else {
		if (all_scenarios == TRUE) {
			setwd("/users/aesin/Desktop/Mowgli/Mowgli_outputs/NO_IG_intersect/All_scenarios")
		} else {
			dir.create("/users/aesin/Desktop/Mowgli/Mowgli_outputs/NO_IG_intersect/Scenarios_1_2", showWarnings = FALSE)
			setwd("/users/aesin/Desktop/Mowgli/Mowgli_outputs/NO_IG_intersect/Scenarios_1_2")
		}
	}
	# write.table(true_HGT_intersect, file = paste0("true_HGT_intersect_t", penalty_value_str, ".tsv"), sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
	# write.table(true_vertical_intersect, file = paste0("true_vertical_intersect_t", penalty_value_str, ".tsv"), sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
	# write.table(penalty_group_num_df, file = paste0("Groups_reconciled_per_penalty_t", penalty_value_str, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

	results_counter = results_counter + 1
}