#!/usr/bin/Rscript
## Compare the results values at different transfer penalty costs and extract those groups which are positive for HGT signal into Anoxy/Geo at all the penalties ##
library(stringr)
library(plyr)
library(dplyr)
library(reshape2)
library(lazyeval)
library(gtools)

MergeByGroup <- function(x, y){
	df <- merge(x, y, by = "Group")
}

###########################################################################
## Set whether the mowgli data parser excluded all those trees that contained IG tips only ##
IG_inclusion	<- TRUE
## Set whether we are testing all the scenarios or just 1 + 2 (those where the root is NOT in Anoxy/Geo) ##
all_scenarios	<- FALSE; ## FALSE gives us ONLY scenarios 1 + 2

###########################################################################
## Get and order the results files by dictionary order (package = gtools) ##
if (IG_inclusion == TRUE) {
	in_path	<- "/users/aesin/Desktop/Mowgli/Mowgli_outputs/Results_Full"
	dir.create("/users/aesin/Desktop/Mowgli/Consistent_HGT_Vertical/Consistent_Full/All_scenarios/No_more_HGT", showWarnings = FALSE)
} else {
	in_path	<- "/users/aesin/Desktop/Mowgli/Mowgli_outputs/Results_NO_IG"
	dir.create("/users/aesin/Desktop/Mowgli/Consistent_HGT_Vertical/Consistent_NO_IG/All_scenarios/No_more_HGT", showWarnings = FALSE)
}

setwd(in_path)
results_files	<- mixedsort(grep("results", dir(), value = TRUE))

while (length(results_files) > 1) {
	## We want to do all the necessary permutations of cross-penalty HGT. E.g. for horizontal we want to do: 3 + 4, 3 + 4 + 5, 3 + 4 + 5 + 6. For vertical, we want to do it from the other end (most permissive vertical is at the highest transfer penalty) ##
	results_counter = 2
	while (results_counter <= length(results_files)) {

		###########################################################################
		## Subset the results files to do all the necessary permutations ##
		setwd(in_path)
		results_horizontal_subset	<- results_files[1:results_counter]
		results_vertical_subset		<- results_files[((length(results_files) + 1) - results_counter):length(results_files)]

		## Prepare the penalty list to hold all the data and vectors for the penalty_group_num_df dataframe ##
		hor_penalty_list	<- vector("list", length(results_horizontal_subset))
		hor_penalty_values	<- vector()
		hor_group_nums		<- vector()

		ver_penalty_list	<- vector("list", length(results_vertical_subset))
		ver_penalty_values	<- vector()
		ver_group_nums		<- vector()

		###########################################################################
		## Process each results file in turn ##
		# results_files = "T3_D2_L1_results.txt" "T4_D2_L1_results.txt" #

		for (index in 1:length(results_horizontal_subset)) {
			## Pull out the individual transfer penalty value for each file in the vertical and horizontal subsets ##
			hor_trans_penalty	<- as.integer(str_sub(regmatches(results_horizontal_subset[index], regexpr("[[:alnum:]+]+", results_horizontal_subset[index])), 2))
			ver_trans_penalty	<- as.integer(str_sub(regmatches(results_vertical_subset[index], regexpr("[[:alnum:]+]+", results_vertical_subset[index])), 2))

			hor_penalty_values	<- c(hor_penalty_values, hor_trans_penalty)
			ver_penalty_values	<- c(ver_penalty_values, ver_trans_penalty)

			## Assign data both to a new table (for individual overview) and then into the penalty_list lists for downstream merging ##
			hor_name	<- paste0("results_table_", hor_trans_penalty)
			hor_table	<- read.table(file = results_horizontal_subset[index], header = FALSE, sep = "\t")

			ver_name	<- paste0("results_table_", ver_trans_penalty)
			ver_table	<- read.table(file = results_vertical_subset[index], header = FALSE, sep = "\t")

			assign(hor_name, hor_table)
			assign(ver_name, ver_table)

			hor_penalty_list[[index]]	<- hor_table
			ver_penalty_list[[index]]	<- ver_table

			colnames(hor_penalty_list[[index]])[1:3]	<- c("Group", paste0("HGT_", hor_trans_penalty), paste0("T", hor_trans_penalty, "_trans_type"))
			colnames(ver_penalty_list[[index]])[1:3]	<- c("Group", paste0("HGT_", ver_trans_penalty), paste0("T", ver_trans_penalty, "_trans_type"))

			## Get the total number of entries (i.e. groups processed by Mowgli) for each table ##
			hor_group_nums	<- c(hor_group_nums, nrow(hor_table))
			ver_group_nums	<- c(ver_group_nums, nrow(ver_table))

		}

		###########################################################################
		## Print the values of the transfer penalties being compared ##

		message("\nTransfer penalties being intersected for HGT: ", paste(hor_penalty_values, collapse = " + "))
		message("Transfer penalties being intersected for Vertical: ", paste(ver_penalty_values, collapse = " + "))

		###########################################################################
		## Create a data frame to show how many groups were reconciled to produce the HGT/Vertical splits at each penalty ##

		hor_penalty_num_group_df	<- as.data.frame(cbind(hor_penalty_values, hor_group_nums))
		ver_penalty_num_group_df	<- as.data.frame(cbind(ver_penalty_values, ver_group_nums))

		colnames(hor_penalty_num_group_df)	<- c("HGT_Penalty", "Number_groups_queried")
		colnames(ver_penalty_num_group_df)	<- c("HGT_Penalty", "Number_groups_queried")

		###########################################################################
		## Produce the table and sort it down to either those that all give positive HGT or positive vertical signal ##

		hor_merged_results	<- Reduce(MergeByGroup, hor_penalty_list)
		ver_merged_results	<- Reduce(MergeByGroup, ver_penalty_list)

		## The set of conditions for the filter function are provided as a string ##
		hor_conditions	<- paste(colnames(hor_merged_results[,c(seq(2,ncol(hor_merged_results), by = 2))]), "== 1", collapse = " & ")
		ver_conditions	<- paste(colnames(ver_merged_results[,c(seq(2,ncol(ver_merged_results), by = 2))]), "== 0", collapse = " & ")

		#################################
		### We also want to extract those groups at each penalty for which HGT is no longer predicted at the next highest penalty

		## Get the last penalty (highest value), and change its condition to be absent ##
		all_penalties_split	<- strsplit(hor_conditions, "&")[[1]]
		last_penalty		<- str_trim(all_penalties_split[(length(all_penalties_split))])
		other_penalties		<- unlist(all_penalties_split)[1:(length(all_penalties_split) - 1)]
		str_sub(last_penalty, -1, -1) <- "0"

		## Prepare the filter string ##
		no_more_HGT_conditions <- paste(c(other_penalties, last_penalty), collapse = "& ")

		#################################

		## If we want all the scenarios, then we do not need to filter further ##
		if (all_scenarios == TRUE) {

			## The condition string can be used with the standard evaluation version of the filter (filter_) function in dplyr ##
			true_hor_intersect	<- filter_(hor_merged_results, hor_conditions)$Group
			true_ver_intersect	<- filter_(ver_merged_results, ver_conditions)$Group

			## Last penalty is now not HGT ##
			no_more_hor_intersect <- filter_(hor_merged_results, no_more_HGT_conditions)$Group

		} else {

			## Here we need to further filter based on the scenario breakdown. Obviously the vertical doesn't change ##
			true_ver_intersect	<- filter_(ver_merged_results, ver_conditions)$Group

			## Secondary filtering step ##
			scenario_filter		<- paste(colnames(hor_merged_results)[c(seq(3,ncol(hor_merged_results), by = 2))], "< 3", collapse = " & ")

			true_hor_intersect	<- filter_(hor_merged_results, hor_conditions)
			true_hor_intersect	<- filter_(true_hor_intersect, scenario_filter)$Group

			## Last penalty is now not HGT ##
			no_more_hor_scenario_filter		<- paste(colnames(hor_merged_results)[c(seq(3,(ncol(hor_merged_results) - 2), by = 2))], "< 3", collapse = " & ")
			message(no_more_hor_scenario_filter)

			no_more_hor_intersect	<- filter_(hor_merged_results, no_more_HGT_conditions)
			no_more_hor_intersect	<- filter_(no_more_hor_intersect, no_more_hor_scenario_filter)$Group
		}



		###########################################################################
		## Prepare string for the output name based on the penalties queried ##

		hor_penalty_value_str	<- paste(hor_penalty_values, collapse = "_t")
		ver_penalty_value_str	<- paste(ver_penalty_values, collapse = "_t")
		no_more_hor_value_str	<- paste(hor_penalty_values[1:(length(hor_penalty_values) - 1)], collapse = "_t")

		if (IG_inclusion == TRUE) {

			if (all_scenarios == TRUE) {
				setwd("/users/aesin/Desktop/Mowgli/Consistent_HGT_Vertical/Consistent_Full/All_scenarios")

			} else {
				dir.create("/users/aesin/Desktop/Mowgli/Consistent_HGT_Vertical/Consistent_Full/Scenarios_1_2", showWarnings = FALSE)
				dir.create("/users/aesin/Desktop/Mowgli/Consistent_HGT_Vertical/Consistent_Full/Scenarios_1_2/No_more_HGT", showWarnings = FALSE)
				setwd("/users/aesin/Desktop/Mowgli/Consistent_HGT_Vertical/Consistent_Full/Scenarios_1_2")
			}

		} else {

			if (all_scenarios == TRUE) {
				setwd("/users/aesin/Desktop/Mowgli/Consistent_HGT_Vertical/Consistent_NO_IG/All_scenarios")

			} else {
				dir.create("/users/aesin/Desktop/Mowgli/Consistent_HGT_Vertical/Consistent_NO_IG/Scenarios_1_2", showWarnings = FALSE)
				dir.create("/users/aesin/Desktop/Mowgli/Consistent_HGT_Vertical/Consistent_NO_IG/Scenarios_1_2/No_more_HGT", showWarnings = FALSE)
				setwd("/users/aesin/Desktop/Mowgli/Consistent_HGT_Vertical/Consistent_NO_IG/Scenarios_1_2")
			}
		}

		write.table(true_hor_intersect, file = paste0("true_HGT_intersect_t", hor_penalty_value_str, ".tsv"), sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)
		write.table(true_ver_intersect, file = paste0("true_vertical_intersect_t", ver_penalty_value_str, ".tsv"), sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)


		write.table(hor_penalty_num_group_df, file = paste0("Groups_reconciled_per_penalty_t", hor_penalty_value_str, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
		write.table(ver_penalty_num_group_df, file = paste0("Groups_reconciled_per_penalty_t", ver_penalty_value_str, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

		write.table(no_more_hor_intersect, file = paste0("No_more_HGT/no_more_HGT_intersect_t", no_more_hor_value_str, ".tsv"), sep = "\n", quote = FALSE, col.names = FALSE, row.names = FALSE)

		results_counter = results_counter + 1
	}

	results_files <- results_files[-length(results_files)]
}
