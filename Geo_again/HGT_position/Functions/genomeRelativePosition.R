#!/usr/bin/env Rscript

genomeRelativePosition	<- function(genePosition, oriStart, oriEnd, oriStrand, genomeLength) {

	if (identical(oriStrand, "Forward")) {
		# Get the position relative to the origin
		relGenePosition	<- genePosition - oriStart

		# If the relative location is negative then it lies upstream of the origin
		if (relGenePosition < 0) {
			relGenePosition	<- relGenePosition + genomeLength
		}
	} else if (identical(oriStrand, "Reverse")) {
		# Get the position relative to the origin
		relGenePosition	<- genePosition - oriEnd

		# If it's negative it lies downstream of the origin. Take the absolute value and leave as is #
		# If it's positive it lies upstream. We subtract it's position from the genome size. This will result in chromStart being > chromEnd #
		relGenePosition	<- ifelse(relGenePosition < 0, abs(relGenePosition), genomeLength - relGenePosition)
	} else {
		stop("Origin strand information is missing")
	}

	## Get the fractional position
	fractGenePosition	<- round(relGenePosition / genomeLength, digits = 8)
	return(fractGenePosition)
}

genomeRelativePosition_format	<- function(row) {

	oriStrand	<- row[["oriStrand"]]
	geneStart	<- as.numeric(row[["gene_start"]])
	genome_l	<- as.numeric(row[["genome_l"]])

	if (identical(oriStrand, "Forward")) {

		oriStart	<- as.numeric(row[["oriStart"]])

		# Get the position relative to the origin
		relGenePosition	<- geneStart - oriStart

		# If the relative location is negative then it lies upstream of the origin
		if (relGenePosition < 0) {
			relGenePosition	<- relGenePosition + genome_l
		}
	} else if (identical(oriStrand, "Reverse")) {

		oriEnd	<- as.numeric(row[["oriEnd"]])

		# Get the position relative to the origin
		relGenePosition	<- geneStart - oriEnd

		# If it's negative it lies downstream of the origin. Take the absolute value and leave as is #
		# If it's positive it lies upstream. We subtract it's position from the genome size. This will result in chromStart being > chromEnd #
		relGenePosition	<- ifelse(relGenePosition < 0, abs(relGenePosition), genome_l - relGenePosition)
	} else {
		stop("Origin strand information is missing")
	}

	## Get the fractional position
	fractGenePosition	<- round(relGenePosition / genome_l, digits = 8)
	return(fractGenePosition)
}