#!/usr/bin/env Rscript

require(pacman)
p_load(tidyverse, fs, wesanderson, viridis)

# Data input paths
master_path			<- "/Users/aesin/Desktop/Staph"

mowgli_path			<- file.path(master_path, "Mowgli")
genome_path			<- file.path(master_path, "Core_genomes")
position_path		<- file.path(master_path, "HGT_position")
# time_path			<- file.path(master_path, "HGT_time")
# sporulation_path	<- file.path(master_path, "Sporulation")

taxdmp_path			<- file.path(master_path, "Genomes", "taxdmp")
# giProcess_path		<- file.path(position_path, "GI_data")

# Paths for scripts and functions
scripts_path		<- "/Users/aesin/Documents/Scripts/Staph/HGT_position"
run_path			<- file.path(scripts_path, "Run")
functions_path		<- file.path(scripts_path, "Functions")

# Master variables
dataTypes			<- c("All", "lHGT", "sHGT", "Ver")
# dataTypes_withAge	<- c(dataTypes, "Old", "Recent")
penalty_list		<- as.character(c(3, 4))

# Penalties (HGT + Ver) that are used for most cases (e.g. COG analysis)
verPenalty			<- "3"
hgtPenalty			<- "4"

# All prot db
allProtDB_path		<- file.path(master_path, "ForStaph_prot_db")

# Location of reused dataObjects
positionData_path	<- file.path(position_path, "DataObjects")
dir_create(positionData_path)
# timeData_path		<- file.path(time_path, "Data")
# supergroupData_path	<- file.path(positionData_path, "SuperGroups")
# if (!dir.exists(supergroupData_path)) dir.create(supergroupData_path, recursive = TRUE)
# Outlier genomes (genomic rearrangements)
# outlierTaxid		<- c(294699, 544556)
outlierTaxid		<- NULL

# Location of figure outputs
numGenomes			<- 25 - length(outlierTaxid)
figureOutput_path	<- file.path(position_path, "Figures", paste0(numGenomes, "Genomes"))
# timeOutput_path		<- file.path(time_path, "Figures")
if (!dir.exists(figureOutput_path)) dir.create(figureOutput_path, recursive = TRUE)
# if (!dir.exists(timeOutput_path)) dir.create(timeOutput_path, recursive = TRUE)

# Standard colors for various things
greenScale2	<- viridis(3, begin = 0.4, end = 0.87, direction = -1)
boundCol	<- wes_palette("Darjeeling1")[1]
dataTypeCols		<- list(
	All = wes_palette("IsleofDogs2")[5],
	HGT = greenScale2[2],
	Ver = wes_palette("Darjeeling1")[3],
	Old = greenScale2[3]
)

# Recent = greenScale2[1],
# Other = wes_palette("IsleofDogs1")[6]