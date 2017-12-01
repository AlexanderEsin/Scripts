---
title: Comparing FastTree and RAxML based reconciliations
date: 24 Nov 2017
output:
  html_document:
    toc: true
    toc_depth: 3
    toc_float: true
    number_sections: true
    code_folding: hide
    theme: cosmo
---

```{r global_options, include = FALSE}
knitr::opts_chunk$set(	fig.width	= 10, 
						fig.height	= 7, 
						fig.path	= "/Users/aesin/Desktop/FastTree/Analysis/Figures", 
						fig.align	= 'center', 
						dpi			= 300, 
						cache.path	= "/Users/aesin/Desktop/FastTree/Analysis/Cache", 
						warning		= TRUE, 
						message		= TRUE,
						tidy		= TRUE)

```

# Packages, functions {.tabset .tabset-fade}
```{r packages, warning = FALSE, message = FALSE}
library(dplyr)
library(stringr)
library(ggplot2)


reprocessEventTips	<- function(event_list) {
	tips_ordered	<- apply(event_list, 1, function(row) {
		vector_tips	<- unlist(str_split(row[4], " "))
		order_tips	<- vector_tips[order(vector_tips)]
		num_tips	<- length(order_tips)
		order_tips	<- paste(order_tips, collapse = " ")
		return(data.frame(Group = row[1], Donor_nodes = row[2], Receptor_nodes = row[3], Tips = order_tips, Num.Tips = num_tips, stringsAsFactors = FALSE))
	})
	tips_ordered_df	<- bind_rows(tips_ordered)
	return(tips_ordered_df)
}

findSplitEvents		<- function(query_events, subject_events) {


	split_events <- lapply(1:nrow(query_events), function(row_index) {
		row				<- query_events[row_index,]
		## Define the necessary columns
		query_group		<- row$Group
		query_num_tips	<- row$Num.Tips
		in_subject		<- subject_events[which(subject_events$Group == query_group),]

		## If the group in the subject set (e.g. RAxML) is not found
		## in the test set (e.g. FT) - return NULL
		if (nrow(in_subject) == 0) {
			return(NULL)
		}

		## Here test whether every tip in the query set reappears in
		## multiple events predicted in the subject set. A success is when
		## every query tip is found in 2 or more events
		query_tips		<- unlist(str_split(row$Tips, " "))
		found_tips_l	<- lapply(query_tips, function(tip) {
			tip_found_index	<- grep(tip, in_subject$Tips)
			if (length(tip_found_index) > 0) {
				return(data.frame(Found = TRUE, Subject_index = tip_found_index, stringsAsFactors = FALSE))
			} else {
				return(data.frame(Found = FALSE, Subject_index = NA, stringsAsFactors = FALSE))
			}
		})
		found_tips_df	<- bind_rows(found_tips_l)

		## Check that all the query tips have been found, if not return
		if (!all(found_tips_df$Found)) {
			return(NULL)
		}

		## Get the events in the subject set that correspond
		## to the tips from the query set
		unique_indices	<- unique(found_tips_df$Subject_index)
		in_subject_trim	<- in_subject[unique_indices,]

		## Finally we want to make sure that the tip sets correspond exactly
		## E.g. tips (A B C) in query might be (A B) (C) in subject,
		## which is acceptable. However (A B) (C D) in subject would not
		subject_trim_tips	<- unlist(str_split(in_subject_trim$Tips, " "))
		if (isTRUE(SameElements(query_tips, subject_trim_tips))) {
			row$type				<- "Query"
			in_subject_trim$type	<- "Subject"
			return(bind_rows(row, in_subject_trim))
		} else {
			return(NULL)
		}
	})

	## Bind together a dataframe - splut_events_df contains all the
	## query and subject events that match
	split_events	<- split_events[lapply(split_events, length) > 0]
	if (length(split_events) == 0) {
		return(NULL)
	}

	split_events_df	<- bind_rows(split_events)

	## Take the initial lists, and remove any split events successfully
	## found
	query_reduct	<- setdiff(query_events, split_events_df[which(split_events_df$type == "Query"),-ncol(split_events_df)])
	subject_reduct	<- setdiff(subject_events, split_events_df[which(split_events_df$type == "Subject"),-ncol(split_events_df)])

	## Number of tips (proteins) that we've reconciled in this split
	## analysis
	split_tips		<- sum(split_events_df$Num.Tips[which(split_events_df$type == "Query")])

	return(list(Query_reduced = query_reduct, Subject_reduced = subject_reduct, Tips.reconciled = split_tips, all_reconciled_df = split_events_df))
}

findSameTips	<- function(events_a, events_b) {

	identicalTips	<- lapply(1:nrow(events_a), function(row_index) {
		row			<- events_a[row_index,]
		group		<- row$Group
		in_subject	<- events_b[which(events_b$Group == group),]

		## If the group in the subject set (e.g. RAxML) is not found
		## in the test set (e.g. FT) - return 0
		if (nrow(in_subject) == 0) {
			return(0)
		}

		row_tips	<- unlist(str_split(row$Tips, " "))
		found_tips_l	<- lapply(row_tips, function(tip) {

			tip_found_index	<- grep(tip, in_subject$Tips)
			if (length(tip_found_index) > 0) {
				return(1)
			} else {
				return(0)
			}
		})
		num_tips_found	<- Reduce("+", found_tips_l)
		return(num_tips_found)
	})

	total_tips_found	<- Reduce("+", identicalTips)
	return(total_tips_found)
}

SameElements <- function(a, b) return(identical(sort(a), sort(b)))

```

# 1. Read in data
```{r penalty_read_in_dta, warning = FALSE, message = FALSE, cache = TRUE}
penalty <- 5

## Read in events
RAX_lHGT_events <- read.table(file = file.path("/Users/aesin/Desktop/FastTree/RAxML_outputs/Refined_events/Events", paste0("T", penalty, "_full_lHGT_events.tsv")), sep = "\t", header = T, stringsAsFactors = F)
RAX_sHGT_events <- read.table(file = file.path("/Users/aesin/Desktop/FastTree/RAxML_outputs/Refined_events/Events", paste0("T", penalty, "_full_sHGT_events.tsv")), sep = "\t", header = T, stringsAsFactors = F)

FT_lHGT_events <- read.table(file = file.path("/Users/aesin/Desktop/FastTree/FastTree_outputs/Refined_events/Events", paste0("T", penalty, "_full_lHGT_events.tsv")), sep = "\t", header = T, stringsAsFactors = F)
FT_sHGT_events <- read.table(file = file.path("/Users/aesin/Desktop/FastTree/FastTree_outputs/Refined_events/Events", paste0("T", penalty, "_full_sHGT_events.tsv")), sep = "\t", header = T, stringsAsFactors = F)

RAX_vert_groups	<- read.table(file = "/Users/aesin/Desktop/FastTree/RAxML_outputs/Constant_events/Ver_const_t4_t5_t6.tsv", sep = "\n", header = FALSE, stringsAsFactors = FALSE)$V1
FT_vert_groups	<- read.table(file = "/Users/aesin/Desktop/FastTree/FastTree_outputs/Constant_events/Ver_const_t4_t5_t6.tsv", sep = "\n", header = FALSE, stringsAsFactors = FALSE)$V1

RAX_const_HGT	<- read.table(file = "/Users/aesin/Desktop/FastTree/RAxML_outputs/Constant_events/HGT_full_tbl_t3_t4_t5.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
FT_const_HGT	<- read.table(file = "/Users/aesin/Desktop/FastTree/FastTree_outputs/Constant_events/HGT_full_tbl_t3_t4_t5.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

RAX_const_Vert	<- read.table(file = "/Users/aesin/Desktop/FastTree/RAxML_outputs/Constant_events/Ver_full_tbl_t3_t4_t5_t6.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
FT_const_Vert	<- read.table(file = "/Users/aesin/Desktop/FastTree/FastTree_outputs/Constant_events/Ver_full_tbl_t3_t4_t5_t6.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
```


# 2. Process lHGT data
```{r lHGT_primary_data, warning = FALSE, message = FALSE, cache = TRUE}

## Process the lHGT events. Reorder the tips in the same
## way for each dataset so the columns can be directly compared
RAX_lHGT_order	<- reprocessEventTips(RAX_lHGT_events)
FT_lHGT_order	<- reprocessEventTips(FT_lHGT_events)

rax_lHGT_event_num	<- nrow(RAX_lHGT_order)
rax_lHGT_group_num	<- length(unique(RAX_lHGT_order$Group))
total_lHGT_rax_tips	<- sum(RAX_lHGT_order$Num.Tips)

ft_lHGT_event_num	<- nrow(FT_lHGT_order)
ft_lHGT_group_num	<- length(unique(FT_lHGT_order$Group))
total_lHGT_ft_tips	<- sum(FT_lHGT_order$Num.Tips)

lhgt_event_num_df		<- data.frame(Method = c("RAxML", "FastTree"), Number = c(rax_lHGT_event_num, ft_lHGT_event_num), Statistic = rep("Number of Events", 2), stringsAsFactors = FALSE)
lhgt_group_num_df		<- data.frame(Method = c("RAxML", "FastTree"), Number = c(rax_lHGT_group_num, ft_lHGT_group_num), Statistic = rep("Number of Groups", 2), stringsAsFactors = FALSE)
lhgt_tip_num_df			<- data.frame(Method = c("RAxML", "FastTree"), Number = c(total_lHGT_rax_tips, total_lHGT_ft_tips), Statistic = rep("Number of Proteins", 2), stringsAsFactors = FALSE)
lhgt_stat_data_tbl		<- bind_rows(list(lhgt_event_num_df, lhgt_group_num_df, lhgt_tip_num_df))
```

```{r lHGT_primary_stats_plot, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}
lHGT_stats.p	<- 	ggplot(data = lhgt_stat_data_tbl, aes(x = Method, y = Number, fill = Statistic)) +
					geom_bar(stat = "identity") +
					facet_wrap(~Statistic, scales = "free") +
					ggtitle("Number of events, groups, and unique proteins between RAxML and FT") +
					theme(plot.title = element_text(hjust = 0.5))
print(lHGT_stats.p)
```


```{r lHGT_common_events, warning = FALSE, message = FALSE, cache = TRUE}
## Completely identical events
lHGT_identical	<- inner_join(RAX_lHGT_order, FT_lHGT_order)
num_lHGT_ident	<- nrow(lHGT_identical)
lhgt_prop_rax	<- num_lHGT_ident / rax_lHGT_event_num
lhgt_prop_ft	<- num_lHGT_ident / ft_lHGT_event_num

## Ignoring donor node identity
lHGT_exdonor		<- inner_join(RAX_lHGT_order[,-2], FT_lHGT_order[,-2])
num_lHGT_exdon		<- nrow(lHGT_exdonor)
lhgt_exdon_prop_rax	<- (num_lHGT_exdon - num_lHGT_ident) / rax_lHGT_event_num
lhgt_exdon_prop_ft	<- (num_lHGT_exdon - num_lHGT_ident) / ft_lHGT_event_num

## Ignoring exact donor and receptor node predictions
## I.e. only group and tips identical
lHGT_exdon_rec			<- inner_join(RAX_lHGT_order[,-2:-3], FT_lHGT_order[,-2:-3])
num_lHGT_exdon_rec		<- nrow(lHGT_exdon_rec)
lhgt_exdon_rec_rax		<- (num_lHGT_exdon_rec - num_lHGT_exdon) / rax_lHGT_event_num
lhgt_exdon_rec_ft		<- (num_lHGT_exdon_rec - num_lHGT_exdon) / ft_lHGT_event_num


## Place data into data frame. Use factors for 'type' so
## bars can be arranged in custom order (according to factor)
factors					<- c("Split", "Ignoring Receptor", "Ignoring Donor", "Completely Identical")
types					<- factor(c(factors[-1], levels = factors)
lHGT_rax_event_ident_df	<- data.frame(Method = rep("RAxML", 3), Fraction = c(lhgt_prop_rax, lhgt_exdon_prop_rax, lhgt_exdon_rec_rax), Type = types, stringsAsFactors = TRUE)
lHGT_ft_event_ident_df	<- data.frame(Method = rep("FastTree", 3), Fraction = c(lhgt_prop_ft, lhgt_exdon_prop_ft, lhgt_exdon_rec_ft), Type = types, stringsAsFactors = TRUE)
lHGT_event_ident_table	<- bind_rows(list(lHGT_rax_event_ident_df, lHGT_ft_event_ident_df))
```


```{r lHGT_common_plot, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}
## Plot
lHGT_events_reconciled.p	<-	ggplot(lHGT_event_ident_table, aes(x = Method, y = Fraction, fill = Type)) +
								geom_bar(stat = "identity") +
								scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.05)) +
								ggtitle("Proportion of all RAxML / FT events found in the other") +
								theme(plot.title = element_text(hjust = 0.5))
print(lHGT_events_reconciled.p)
```



```{r lHGT_split_events, warning = FALSE, message = FALSE, cache = TRUE}
## For both sets, remove the events that we consider identical
## 'exdon_rec'
RAX_lHGT_missing		<- setdiff(RAX_lHGT_order, semi_join(RAX_lHGT_order, lHGT_exdon_rec))
FT_lHGT_missing			<- setdiff(FT_lHGT_order, semi_join(FT_lHGT_order, lHGT_exdon_rec))

## Find events that are split from RAxML to FT and vice versa, i.e.:
## situations where 1 event with 10 tips (in RAX) might be split into 2 events
## with 6 and 4 tips (in FT)
RAX_split_in_FT_lHGT	<- findSplitEvents(RAX_lHGT_missing, FT_lHGT_missing)
RAX_lHGT_missing_2		<- RAX_split_in_FT_lHGT$Query_reduced
FT_lHGT_missing_2		<- RAX_split_in_FT_lHGT$Subject_reduced

FT_split_in_RAX_lHGT	<- findSplitEvents(FT_lHGT_missing_2, RAX_lHGT_missing_2)
RAX_lHGT_missing_fin	<- FT_split_in_RAX_lHGT$Subject_reduced
FT_lHGT_missing_fin		<- FT_split_in_RAX_lHGT$Query_reduced

## Count the number of extra/missing lHGT events even after
## accounting for split events
num_RAXlHGT_missing		<- nrow(RAX_lHGT_missing_fin)
num_FTlHGT_missing		<- nrow(FT_lHGT_missing_fin)

lhgt_prop_split_rax		<- (nrow(RAX_lHGT_missing) - num_RAXlHGT_missing) / rax_lHGT_event_num
lhgt_prop_split_ft		<- (nrow(FT_lHGT_missing) - num_FTlHGT_missing) / ft_lHGT_event_num

types					<- factor("Split", levels = factors)
lhgt_prop_split_df		<- data.frame(Method = c("RAxML", "FastTree"), Fraction = c(lhgt_prop_split_rax, lhgt_prop_split_ft), Type = types, stringsAsFactors = TRUE)
lHGT_event_ident_table	<- rbind(lHGT_event_ident_table, lhgt_prop_split_df)
```


```{r lHGT_split_events_plot, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}
lHGT_events_split.p		<-	ggplot(lHGT_event_ident_table, aes(x = Method, y = Fraction, fill = Type)) +
							geom_bar(stat = "identity") +
							scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.05)) +
							ggtitle("Proportion of all RAxML / FT events found in the other") +
							theme(plot.title = element_text(hjust = 0.5))
print(lHGT_events_split.p)
```



```{r lHGT_tips_accounted, warning = FALSE, message = FALSE, cache = TRUE}
## How many tips are missing
RAX_lHGT_tip_missing	<- sum(RAX_lHGT_missing_fin$Num.Tips)
FT_lHGT_tip_missing		<- sum(FT_lHGT_missing_fin$Num.Tips)


## Lastly, let's see how many proteins are actually missing. E.g.
## if a gene family is predicted in both the RAX and FT sets, but they are
## not identical or split predictions, how many lHGT proteins from lHGT
## are also predicted in the FT set?
tips_shared_lHGT_fin	<- findSameTips(RAX_lHGT_missing_fin, FT_lHGT_missing_fin)

## RAxML predicted lHGT tips account for in FT
RAX_lHGT_tips_accounted	<- (total_lHGT_rax_tips - (RAX_lHGT_tip_missing - tips_shared_lHGT_fin)) / total_lHGT_rax_tips
## FT predicted lHGT tips account for in RAxML
FT_lHGT_tips_accounted	<- (total_lHGT_ft_tips - (FT_lHGT_tip_missing - tips_shared_lHGT_fin)) / total_lHGT_ft_tips

## Put into dataframe for plotting
lHGT_tips_stats_df		<- data.frame(Method = c("RAxML", "FastTree"), Number.Of.Tips = c(RAX_lHGT_tips_accounted, FT_lHGT_tips_accounted))
```

```{r lHGT_tips_plot, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}
lHGT_tip_number.p	<-	ggplot(lHGT_tips_stats_df, aes(x = Method, y = Number.Of.Tips, fill = Method)) +
						geom_bar(stat = "identity") +
						scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.05)) +
						ggtitle("Proportion of proteins shared between RAxML / FT predictions") +
						theme(plot.title = element_text(hjust = 0.5))
```



# 3. Process sHGT data
```{r sHGT_primary_data, warning = FALSE, message = FALSE, cache = TRUE}

## Process the sHGT events. Reorder the tips in the same
## way for each dataset so the columns can be directly compared
RAX_sHGT_order	<- reprocessEventTips(RAX_sHGT_events)
FT_sHGT_order	<- reprocessEventTips(FT_sHGT_events)

rax_sHGT_event_num	<- nrow(RAX_sHGT_order)
rax_sHGT_group_num	<- length(unique(RAX_sHGT_order$Group))
total_sHGT_rax_tips	<- sum(RAX_sHGT_order$Num.Tips)

ft_sHGT_event_num	<- nrow(FT_sHGT_order)
ft_sHGT_group_num	<- length(unique(FT_sHGT_order$Group))
total_sHGT_ft_tips	<- sum(FT_sHGT_order$Num.Tips)

shgt_event_num_df		<- data.frame(Method = c("RAxML", "FastTree"), Number = c(rax_sHGT_event_num, ft_sHGT_event_num), Statistic = rep("Number of Events", 2), stringsAsFactors = FALSE)
shgt_group_num_df		<- data.frame(Method = c("RAxML", "FastTree"), Number = c(rax_sHGT_group_num, ft_sHGT_group_num), Statistic = rep("Number of Groups", 2), stringsAsFactors = FALSE)
shgt_tip_num_df			<- data.frame(Method = c("RAxML", "FastTree"), Number = c(total_sHGT_rax_tips, total_sHGT_ft_tips), Statistic = rep("Number of Proteins", 2), stringsAsFactors = FALSE)
shgt_stat_data_tbl		<- bind_rows(list(shgt_event_num_df, shgt_group_num_df, shgt_tip_num_df))
```


```{r sHGT_primary_stats_plot, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}
sHGT_stats.p	<- 	ggplot(data = shgt_stat_data_tbl, aes(x = Method, y = Number, fill = Statistic)) +
					geom_bar(stat = "identity") +
					facet_wrap(~Statistic, scales = "free") +
					ggtitle("Number of events, groups, and unique proteins between RAxML and FT") +
					theme(plot.title = element_text(hjust = 0.5))
pritn(sHGT_stats.p)					
```

```{r sHGT_common_events, warning = FALSE, message = FALSE, cache = TRUE}
## Completely identical events
sHGT_identical	<- inner_join(RAX_sHGT_order, FT_sHGT_order)
num_sHGT_ident	<- nrow(sHGT_identical)
shgt_prop_rax	<- num_sHGT_ident / rax_sHGT_event_num
shgt_prop_ft	<- num_sHGT_ident / ft_sHGT_event_num

## Ignoring donor node identity
sHGT_exdonor		<- inner_join(RAX_sHGT_order[,-2], FT_sHGT_order[,-2])
num_sHGT_exdon		<- nrow(sHGT_exdonor)
shgt_exdon_prop_rax	<- (num_sHGT_exdon - num_sHGT_ident) / rax_sHGT_event_num
shgt_exdon_prop_ft	<- (num_sHGT_exdon - num_sHGT_ident) / ft_sHGT_event_num

## Ignoring exact donor and receptor node predictions
## I.e. group and tips identical
sHGT_exdon_rec		<- inner_join(RAX_sHGT_order[,-2:-3], FT_sHGT_order[,-2:-3])
num_sHGT_exdon_rec	<- nrow(sHGT_exdon_rec)
shgt_exdon_rec_rax		<- (num_sHGT_exdon_rec - num_sHGT_exdon) / rax_sHGT_event_num
shgt_exdon_rec_ft		<- (num_sHGT_exdon_rec - num_sHGT_exdon) / ft_sHGT_event_num

## Into a dataframe...
factors					<- c("Split", "Ignoring Receptor", "Ignoring Donor", "Completely Identical")
types					<- factor(factors[-1], levels = factors)
sHGT_rax_event_ident_df	<- data.frame(Method = rep("RAxML", 3), Fraction = c(shgt_prop_rax, shgt_exdon_prop_rax, shgt_exdon_rec_rax), Type = types, stringsAsFactors = TRUE)
sHGT_ft_event_ident_df	<- data.frame(Method = rep("FastTree", 3), Fraction = c(shgt_prop_ft, shgt_exdon_prop_ft, shgt_exdon_rec_ft), Type = types, stringsAsFactors = TRUE)
sHGT_event_ident_table	<- bind_rows(list(sHGT_rax_event_ident_df, sHGT_ft_event_ident_df))
```


```{r sHGT_split_events, warning = FALSE, message = FALSE, cache = TRUE}
## For both sets, remove the events that we consider identical
## 'exdon_rec'
RAX_sHGT_missing		<- setdiff(RAX_sHGT_order, semi_join(RAX_sHGT_order, sHGT_exdon_rec))
FT_sHGT_missing			<- setdiff(FT_sHGT_order, semi_join(FT_sHGT_order, sHGT_exdon_rec))


## Find events that are split from RAxML to FT and vice versa, i.e.:
## situations where 1 event with 10 tips (in RAX) might be split into 2 events
## with 6 and 4 tips (in FT)
RAX_split_in_FT_sHGT	<- findSplitEvents(RAX_sHGT_missing, FT_sHGT_missing)
RAX_sHGT_missing_2		<- RAX_split_in_FT_sHGT$Query_reduced
FT_sHGT_missing_2		<- RAX_split_in_FT_sHGT$Subject_reduced

FT_split_in_RAX_sHGT	<- findSplitEvents(FT_sHGT_missing_2, RAX_sHGT_missing_2)
RAX_sHGT_missing_fin	<- FT_split_in_RAX_sHGT$Subject_reduced
FT_sHGT_missing_fin		<- FT_split_in_RAX_sHGT$Query_reduced



shgt_prop_split_rax		<- (nrow(RAX_sHGT_missing) - nrow(RAX_sHGT_missing_fin)) / rax_sHGT_event_num
shgt_prop_split_ft		<- (nrow(FT_sHGT_missing) - nrow(FT_sHGT_missing_fin)) / ft_sHGT_event_num

types					<- factor("Split", levels = factors)
shgt_prop_split_df		<- data.frame(Method = c("RAxML", "FastTree"), Fraction = c(shgt_prop_split_rax, shgt_prop_split_ft), Type = types, stringsAsFactors = TRUE)
sHGT_event_ident_table	<- rbind(sHGT_event_ident_table, shgt_prop_split_df)
```


```{r sHGT_split_events_plot, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}
sHGT_events_split.p	<-	ggplot(sHGT_event_ident_table, aes(x = Method, y = Fraction, fill = Type)) +
						geom_bar(stat = "identity") +
						scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.05)) +
						ggtitle("Proportion of all RAxML / FT events found in the other") +
						theme(plot.title = element_text(hjust = 0.5))
print(sHGT_events_split.p)
```



```{r lHGT_tips_accounted, warning = FALSE, message = FALSE, cache = TRUE}
## How many tips are missing
RAX_sHGT_tip_missing	<- sum(RAX_sHGT_missing_fin$Num.Tips)
FT_sHGT_tip_missing		<- sum(FT_sHGT_missing_fin$Num.Tips)

## Lastly, let's see how many proteins are actually missing. E.g.
## if a gene family is predicted in both the RAX and FT sets, but they are
## not identical or split predictions, how many lHGT proteins from lHGT
## are also predicted in the FT set?
tips_shared_sHGT_fin		<- findSameTips(RAX_sHGT_missing_fin, FT_sHGT_missing_fin)

## RAxML predicted lHGT tips account for in FT
RAX_sHGT_tips_accounted	<- (total_sHGT_rax_tips - (RAX_sHGT_tip_missing - tips_shared_sHGT_fin)) / total_sHGT_rax_tips
## FT predicted lHGT tips account for in RAxML
FT_sHGT_tips_accounted	<- (total_sHGT_ft_tips - (FT_sHGT_tip_missing - tips_shared_sHGT_fin)) / total_sHGT_ft_tips

## Put into dataframe
sHGT_tips_stats_df		<- data.frame(Method = c("RAxML", "FastTree"), Number.Of.Tips = c(RAX_sHGT_tips_accounted, FT_sHGT_tips_accounted))
```


```{r sHGT_tips_plot, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}
sHGT_tip_number.p	<-	ggplot(sHGT_tips_stats_df, aes(x = Method, y = Number.Of.Tips)) +
						geom_bar(stat = "identity") +
						scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.05)) +
						ggtitle("Proportion of proteins shared between RAxML / FT predictions") +
						theme(plot.title = element_text(hjust = 0.5))
print(sHGT_tip_number.p)
```						





# 4. Fate of remaining lHGT events predicted by RAxML
Where are the missing RAxML lHGT events? We expect 241 lHGT events based on the RAxML, but even accounting for split tips (one RAxML events into two FT events, or two+ RAxML events into one FT event), we still have 53 outstanding events.
## RAxML lHGTs -> FastTree sHGTs

```{r RAxML_lHGTs_to_sHGTs, warning = FALSE, message = FALSE, cache = TRUE}
## Some RAxML lHGT events could be predicted as sHGT in FT
RAX_lHGT_to_FT_sHGT	<- inner_join(RAX_lHGT_missing_fin[,-2:-3], FT_sHGT_missing_fin[,-2:-3])
RAX_lHGT_remaining	<- setdiff(RAX_lHGT_missing_fin, semi_join(RAX_lHGT_missing_fin, RAX_lHGT_to_FT_sHGT))

# 4 RAxML lHGT events are found to be sHGT events in FT prediction
num_lHGT_to_sHGT	<- nrow(RAX_lHGT_missing_fin) - nrow(RAX_lHGT_remaining)
fract_lHGT_to_sHGT	<- num_lHGT_to_sHGT / num_RAXlHGT_missing
```

## RAxML lHGTs -> FastTree Vertical
```{r RAxML_lHGTs_to_Vertical, warning = FALSE, message = FALSE, cache = TRUE}
## Some RAxML lHGT events could be predicted as Vertical in FT
RAX_lHGT_to_FT_vert_full	<- RAX_lHGT_remaining[which(RAX_lHGT_remaining$Group %in% as.numeric(FT_vert_groups)),]

# 7 events, corresponding to 7 groups 
num_RAX_lHGT_to_FT_vert		<- nrow(RAX_lHGT_to_FT_vert_full)
fract_lHGT_to_FT_vert		<- num_RAX_lHGT_to_FT_vert / num_RAXlHGT_missing
```

## RAxML lHGTs -> Inconstant & Inconsistent
```{r RAxML_lHGTs_to_inconst_inconsist, warning = FALSE, message = FALSE, cache = TRUE}
## The remaining events, corresponding to groups are either not-constant or do not have
## consistently predicted lHGT or sHGT
RAX_lHGT_remaining			<- setdiff(RAX_lHGT_remaining, RAX_lHGT_to_FT_vert_full)
RAX_lhgt_remain_group		<- unique(RAX_lHGT_remaining$Group)

## We have 42 events remaining, from 33 groups
## Find whether HGT was constantly predicted for these remaining groups across penalties
RAX_lHGT_const_remain_tbl		<- FT_const_HGT[which(FT_const_HGT$Gene.Family %in% RAX_lhgt_remain_group),c(1,2,4,6)]
RAX_lHGT_const_remain_tbl$rSum	<- rowSums(RAX_lHGT_const_remain_tbl[,-1])

# There are 8 groups which are not constant with HGT prediction (or predict root)
RAX_lHGT_inconst_remain		<- RAX_lHGT_const_remain_tbl[which(RAX_lHGT_const_remain_tbl$rSum < 3),]
RAX_lHGT_inconst_groups		<- RAX_lHGT_inconst_remain$Gene.Family
RAX_lHGT_inconst_events		<- RAX_lHGT_remaining[which(RAX_lHGT_remaining$Group %in% RAX_lHGT_inconst_groups),]
num_RAX_inconst_events		<- nrow(RAX_lHGT_inconst_events)
fract_RAX_inconst_events	<- num_RAX_inconst_events / num_RAXlHGT_missing

## We now have 34 events remaining, consisting of 119 tips that do not show CONSISTENT lHGT
RAX_lHGT_remaining			<- setdiff(RAX_lHGT_remaining, RAX_lHGT_inconst_events)
num_RAX_inconsist_events	<- nrow(RAX_lHGT_remaining)
fract_RAX_inconsist_events	<- num_RAX_inconsist_events / num_RAXlHGT_missing
```

```{r RAxML_lHGT_fate_plot, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}
factors					<- c("Long To Short", "Long to Vertical", "Inconstant", "Inconsistent")
types					<- factor(factors, levels = factors)
RAX_lHGT_fates_df		<- data.frame(Method = rep("RAxML", 4), Fraction = c(fract_lHGT_to_sHGT, fract_lHGT_to_FT_vert, fract_RAX_inconst_events, fract_RAX_inconsist_events), Type = types)

RAX_lHGT_fates.p		<-	ggplot(RAX_lHGT_fates_df, aes(x = Method, y = Fraction, fill = Type)) +
							geom_bar(stat = "identity") +
							scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.05)) +
							ggtitle("Fate of RAxML lHGT events") +
							theme(plot.title = element_text(hjust = 0.5))
print(RAX_lHGT_fates.p)
```



# 4. Fate of extra lHGT events predicted by FastTree
221 lHGT events are predicted by FastTree. After trimming all events shared with the RAxML set (even with splitting), we still have 29 extra lHGT events. Where did these events come from?

## FastTree lHGTS <- RAxML shGTS
```{r FastTree_lHGTS_from_sHGTS, warning = FALSE, message = FALSE, cache = TRUE}
## Some FT lHGT events could come from RAxML sHGTS
FT_lHGT_from_RAX_sHGT	<- inner_join(FT_lHGT_missing_fin[,-2:-3], RAX_sHGT_missing_fin[,-2:-3])
FT_lHGT_remaining		<- setdiff(FT_lHGT_missing_fin, semi_join(FT_lHGT_missing_fin, FT_lHGT_from_RAX_sHGT))

# 6 FT lHGT events are coming from RAX sHGT. 23 remain unexplained
num_FT_lHGT_from_sHGT		<- nrow(FT_lHGT_missing_fin) - nrow(FT_lHGT_remaining)
fract_FT_lHGT_from_sHGT		<- num_lHGT_from_sHGT / num_FTlHGT_missing
```

## FastTree lHGTS <- RAxML Vertical
```{r FastTree_lHGTS_from_Vertical, warning = FALSE, message = FALSE, cache = TRUE}
## Some could have from the RAxML vertical set
FT_lHGT_from_RAX_vert		<- FT_lHGT_remaining[which(FT_lHGT_remaining$Group %in% as.numeric(RAX_vert_groups)),]
FT_lHGT_from_RAX_v_groups	<- unique(FT_lHGT_from_RAX_vert$Group)

# 5 events, corresponding to 5 groups 
num_FT_lHGT_from_RAX_vert	<- nrow(FT_lHGT_from_RAX_vert)
num_FT_lHGT_grp_from_RAX_v	<- length(FT_lHGT_from_RAX_v_groups)
fract_FT_lHGT_from_RAX_vert	<- num_FT_lHGT_from_RAX_vert / num_FTlHGT_missing
```

Every event picked from RAxML must have been consistent. However, the encompassing gene families could have had inconsistent events that might become consistent in the FastTree reconciliations
```{r FastTree_lHGTs_from_inconst_inconsist, warning = FALSE, message = FALSE, cache = TRUE}
## The remaining events, corresponding to groups are either not-constant or do not have
## consistently predicted lHGT or sHGT
FT_lHGT_remaining			<- setdiff(FT_lHGT_remaining, FT_lHGT_from_RAX_vert)
FT_lHGT_remain_group		<- unique(FT_lHGT_remaining$Group)
num_FT_lHGT_remain			<- nrow(FT_lHGT_remaining)
fract_FT_lHGT_inconsistent	<- num_FT_lHGT_remain / num_FTlHGT_missing

## There can be no incostant RAxML events by our definition group choice
fract_FT_lHGT_inconstant	<- 0
````

```{r FastTree_lHGT_fate_plot, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}
factors					<- c("Long To Short", "Long to Vertical", "Inconstant", "Inconsistent")
types					<- factor(factors, levels = factors)
FT_lHGT_fates_df		<- data.frame(Method = rep("FastTree", 4), Fraction = c(fract_FT_lHGT_from_RAX_vert, fract_FT_lHGT_from_sHGT, fract_FT_lHGT_inconstant, fract_FT_lHGT_inconsistent), Type = types)

both_lHGT_fates_df		<- bind_rows(RAX_lHGT_fates_df, FT_lHGT_fates_df)

both_lHGT_fates.p		<-	ggplot(both_lHGT_fates_df, aes(x = Method, y = Fraction, fill = Type)) +
						geom_bar(stat = "identity") +
						scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.05)) +
						ggtitle("Fate of RAxML and FastTree lHGT events") +
						theme(plot.title = element_text(hjust = 0.5))
print(both_lHGT_fates.p)
```


# 5. Comparing Vertical groups
We picked 300 gene families with a vertical Geobacillus history - this means that no transfer was constantly detected in the family across transfer penalties: 4, 5, 6
```{r confirm_vert_assignment, warning = FALSE, message = FALSE, cache = TRUE}
## 300 RAxML groups (as selected) & 264 FT vert groups
num_RAX_vert_groups	<- length(RAX_vert_groups)
num_FT_vert_groups	<- length(FT_vert_groups)
## How many are shared?
shared_vert_groups	<- intersect(RAX_vert_groups, FT_vert_groups)
num_vert_shared		<- length(shared_vert_groups)

method_types		<- factor(c("RAxML", "FastTree", "Shared"), levels = c("RAxML", "FastTree", "Shared"))
vert_stats_df		<- data.frame(Method = method_types, Number.Of.Groups = c(num_RAX_vert_groups, num_FT_vert_groups, num_vert_shared))
```

```{r vert_assignment_plot, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}

vert_stats.p	<-	ggplot(vert_stats_df, aes(x = Method, y = Number.Of.Groups, fill = Method)) +
					geom_bar(stat = "identity") +
					ggtitle("RAxML and FastTree vertical events") +
					theme(plot.title = element_text(hjust = 0.5))
print(vert_stats.p)
```

In the FastTree set we constantly find `r num_FT_vert_groups` vertical groups, of these we find `r num_vert_shared` are shared with the RAxML set.
Thus, we need to account for `r num_RAX_vert_groups - num_vert_shared` RAxML groups that are no longer vertical and `r num_FT_vert_groups - num_vert_shared` extra FastTree predicted vertical groups.
We already know from above that `r num_FT_lHGT_grp_from_RAX_v` FastTree lHGT groups came from vertical RAxML groups.

```{r missing_RAxML_vert_fates, warning = FALSE, message = FALSE, cache = TRUE}
## RAxML groups that are not predicted as vertical in FT set
missing_RAX_vert	<- RAX_vert_groups[!RAX_vert_groups %in% shared_vert_groups]
num_all_RAX_missing	<- length(missing_RAX_vert)

## Some groups may be predicted to have a root in AG - in which case they are discounted (inconstant)
## Check across penalties 4, 5, 6
FT_vert_const_root_tbl	<- FT_const_Vert[which(FT_const_Vert$Gene.Family %in% missing_RAX_vert),c(1,5,7,9)]
FT_AG_root_groups		<- FT_vert_const_root_tbl$Gene.Family[apply(FT_vert_const_root_tbl[,-1], 1, any)]
# 2 former RAxML vertical groups are predicted as having tree root in AG
num_FT_AG_root_groups	<- length(FT_AG_root_groups)
# Update missing_RAX list
missing_RAX_vert		<- missing_RAX_vert[-which(missing_RAX_vert %in% FT_AG_root_groups)]


## All the entries in this table must be inconstant for verticality (i.e. not all 0s). We want to find those that are
## inconstant also for HGT (must be across penalties 4, 5 only). Such a set would be neither predicted as HGT or vertical
FT_vert_const_tbl		<- FT_const_Vert[which(FT_const_Vert$Gene.Family %in% missing_RAX_vert),c(1,2,4,6,8)]
FT_vert_const_tbl$rSum	<- rowSums(FT_vert_const_tbl[,2:4])

## Anything with a score of less than 3 in these columns would be inconsistent vertically and for HGT
FT_vert_inconst			<- FT_vert_const_tbl[which(FT_vert_const_tbl$rSum < 3),]
num_FT_vert_inconst		<- nrow(FT_vert_inconst)
FT_vert_inconst_groups	<- FT_vert_inconst$Gene.Family

## Which missing RAxML vert groups are neither FT lHGTS, are AG root, or inconstant for either HGT or Vert (so sHGTs or inconsistent)
missing_RAX_remain		<- missing_RAX_vert[-which(missing_RAX_vert %in% as.numeric(c(FT_lHGT_from_RAX_v_groups, FT_AG_root_groups, FT_vert_inconst_groups)))]
## 11 groups remain to be accounted for
num_missing_RAX_remain	<- length(missing_RAX_remain)


## Are these 11 groups sHGTs?
RAX_vert_to_sHGT		<- missing_RAX_remain[which(missing_RAX_remain %in% as.numeric(unique(FT_sHGT_order$Group)))]
## 10 are, so 1 must be inconsistent
num_RAX_vert_to_sHGT	<- length(RAX_vert_to_sHGT)
num_RAX_vert_inconsist	<- num_missing_RAX_remain - num_RAX_vert_to_sHGT

## Factors
fates				<- c("RAxML Vert Missing", "Inconstant", "AG_root", "To sHGT", "To lHGT", "Inconsistent")
fates.factor		<- factor(fates, levels = fates)
types				<- c("All", rep("Not HGT/Vert", 2), rep("as HGT", 3))
type.factor			<- factor(types, levels = unique(types))

RAxML_vert_fate_df	<-	data.frame(
						Fate	= fates.factor, 
						Number	= c(num_all_RAX_missing,
									num_FT_vert_inconst,
									num_FT_AG_root_groups,
									num_RAX_vert_to_sHGT,
									length(FT_lHGT_from_RAX_v_groups),
									num_RAX_vert_inconsist),
						Type 	= type.factor
						)

```

```{r RAxML_vert_missing_fates.p, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}
RAX_missing_fates.p	<-	ggplot(RAxML_vert_fate_df, aes(x = Fate, y = Number, fill = Type)) +
						geom_bar(stat = "identity") +
						ggtitle("Fate of missing RAxML vertical events") +
						theme(plot.title = element_text(hjust = 0.5))
print(RAX_missing_fates.p)
```


```{r FT_extra_vert_groups, warning = FALSE, message = FALSE, cache = TRUE, dev = 'png'}
## Where do the extra vertical groups come from?
## 7 come from RAxML lHGTs / 7 come from RAxML sHGTs
extra_FT_vert		<- FT_vert_groups[!FT_vert_groups %in% shared_vert_groups]
num_extra_FT_vert	<- length(extra_FT_vert)

extra_FT_from_lHGT		<- intersect(extra_FT_vert, as.numeric(RAX_lHGT_order$Group))
num_ext_FT_from_lHGT	<- length(extra_FT_from_lHGT)

extra_FT_from_sHGT		<- intersect(extra_FT_vert, as.numeric(RAX_sHGT_order$Group))
num_ext_FT_from_sHGT	<- length(extra_FT_from_sHGT)

extra_FT_vert_df		<- data.frame(Fate = c("All extra", "From sHGT", "From lHGT"), Number = c(num_extra_FT_vert, num_ext_FT_from_sHGT, num_ext_FT_from_lHGT))

extra_FT_fates.p	<-	ggplot(extra_FT_vert_df, aes(x = Fate, y = Number, fill = Fate)) +
						geom_bar(stat = "identity") +
						ggtitle("Fate of extra FastTree vertical events") +
						theme(plot.title = element_text(hjust = 0.5))
print(extra_FT_fates.p)
```
















