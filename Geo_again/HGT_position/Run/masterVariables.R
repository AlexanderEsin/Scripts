#!/usr/bin/env Rscript
# Data input paths
master_path			<- "/Users/aesin/Desktop/Geo_again"
genome_path			<- file.path(master_path, "Genomes")
position_path		<- file.path(master_path, "HGT_position")
taxdmp_path			<- file.path(genome_path, "taxdmp")

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

# Location of figure outputs
figureOutput_path	<- file.path(position_path, "Figures")

# Genomic subcompartments list
subDivision_list	<- list(
	Ori = c(0, 0.05),
	nearOri = c(0.05, 0.15),
	farOri = c(0.15, 0.25),
	farTer = c(0.25, 0.35),
	nearTer = c(0.35, 0.45),
	Ter = c(0.45, 0.5)
)