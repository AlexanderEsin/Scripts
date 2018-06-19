#!/usr/bin/env Rscript

require(pacman)
p_load("wesanderson")

# Data input paths
master_path			<- "/Users/aesin/Desktop/Geo_again"
genome_path			<- file.path(master_path, "Genomes")
position_path		<- file.path(master_path, "HGT_position")
sporulation_path	<- file.path(master_path, "Sporulation")
taxdmp_path			<- file.path(genome_path, "taxdmp")
giProcess_path		<- file.path(position_path, "GI_data")

# Paths for scripts and functions
scripts_path		<- "/Users/aesin/Documents/Scripts/Geo_again/HGT_position"
run_path			<- file.path(scripts_path, "Run")
functions_path		<- file.path(scripts_path, "Functions")

# Master variables
dataTypes			<- c("All", "lHGT", "sHGT", "Ver")
dataTypes_withAge	<- c(dataTypes, "Old", "Recent")
penalty_list		<- as.character(c(3, 4, 5, 6))

# Penalties (HGT + Ver) that are used for most cases (e.g. COG analysis)
verPenalty			<- "3"
hgtPenalty			<- "4"

# All prot db
allProtDB_path		<- file.path(master_path, "All_prot_db_new")

# Location of reused dataObjects
positionData_path	<- file.path(position_path, "DataObjects")
supergroupData_path	<- file.path(positionData_path, "SuperGroups")
if (!dir.exists(supergroupData_path)) dir.create(supergroupData_path, recursive = TRUE)

# Outlier genomes (genomic rearrangements)
outlierTaxid		<- c(294699, 544556)

# Location of figure outputs
figureOutput_path	<- file.path(position_path, "Figures", "23Genomes")

# Standard colors for various things
dataTypeCols		<- list(
	All = wes_palette("IsleofDogs2")[5],
	HGT = wes_palette("Darjeeling1")[2],
	Ver = wes_palette("Darjeeling1")[3],
	Old = wes_palette("GrandBudapest1")[2],
	Recent = wes_palette("Darjeeling1")[5]
)
