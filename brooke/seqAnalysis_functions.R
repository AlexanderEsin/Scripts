## // Split into n chunks // ##
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 


## // Functions to create random DNA sequences 
## based on the same pattern as the library // ##

# Pattern for the nonet repeats
OctNonet_pat	<- function(times, bases = c("A", "C", "G", "T")) {

	OctNonet_multi	<- lapply(1:times, function(x) {
		# Make the 8 x 9 base sets
		total_seq <- lapply(1:8, function(x) {
			c(	sample(bases[2:4], 1),
				sample(bases[1:4], 1),
				bases[4],
				sample(bases[1:4], 1),
				sample(bases[2:4], 2, replace = TRUE),
				sample(bases[1:3], 1),
				sample(bases[1:4], 1),
				sample(bases[3:4], 1)
			)
		})
		total_seq_str	<- str_c(unlist(total_seq), collapse = "")
	})
	return(OctNonet_multi)
}

# Pattern for the trip repeats
SoloTrip_pat	<- function(times, bases = c("A", "C", "G", "T")) {

	SoloTrip_multi	<- lapply(1:times, function(x) {
		total_seq	<-c (
			sample(bases[1:4], 1),
			sample(bases[1:4], 1),
			sample(bases[c(2,4)], 1)
		)
		total_seq_str	<- str_c(total_seq, collapse = "")
	})	
	return(SoloTrip_multi)
}

# Assemble each random sequence according to library rules
fullAssem_pat	<- function(times, midConsensus = "GGAGATAGTACCTCGTCTACACCT") {
	full_list	<- Map(list,
		OctNonet_pat(times),
		SoloTrip_pat(times),
		OctNonet_pat(times),
		SoloTrip_pat(times),
		rep(midConsensus, times),
		OctNonet_pat(times),
		SoloTrip_pat(times),
		OctNonet_pat(times),
		SoloTrip_pat(times))
	
	full_str	<- lapply(full_list, function(sub_list) {
		full_str	<- str_c(unlist(sub_list), collapse = "")
	})
	return(unlist(full_str))
}

# Wrapper function that produces a DNA stringset of desired number of random sequences
makeRandSeqs	<- function(seqNumber = 100000) {
	subsets		<- seqNumber / 10000

	seqSubsets	<- lapply(1:subsets, function(subset) {

		seqSet	<- DNAStringSet(fullAssem_pat(10000))
		message(paste0("\tMade ", as.integer(subset*10000), " // ", as.integer(seqNumber), " random sequences"))
		return(seqSet)
	})
	combSets	<- Reduce(c, seqSubsets)
	names(combSets)	<- seq(1:length(combSets))
	return(combSets)
}



## // Update functions for the dagLogo package - mostly for speed // ##
formatSequence_AE		<- function (seq, proteome) {

	width <- unique(unlist(lapply(seq, nchar)))
	if (length(width) > 1) {
		stop("seq must be characters with same length", call. = FALSE)
	}
	upstreamOffset		<- floor(width/2)
	downstreamOffset	<- width - upstreamOffset - 1
	
	center	<- upstreamOffset + 1

	## Define the anchor AA
	anchorAA	<- str_sub(seq, center, center)
	anchorPos	<- rep(1, length(seq))

	dat <- data.frame(IDs = NA, anchorAA = anchorAA, anchorPos = anchorPos, 
		peptide = seq, anchor = anchorAA)

	seqchar	<- do.call(rbind, str_split(seq, ""))

	new("dagPeptides", data = dat, peptides = seqchar, upstreamOffset = upstreamOffset, downstreamOffset = downstreamOffset, type = "")
}

buildBackgroundModel_AE	<- function (dagPeptides, permutationSize = 30L, rand.seed = 1, replacement = FALSE, proteome) {
	if (missing(dagPeptides) || class(dagPeptides) != "dagPeptides") {
		stop("dagPeptides should be an object of dagPeptides.",
			 call. = FALSE)
	}

	permutationSize <- as.integer(permutationSize)
	if (permutationSize < 2) stop("permutationSize should be greater than 1")

	length  <- dagPeptides@upstreamOffset + dagPeptides@downstreamOffset + 1   
	matches <- proteome@proteome$SEQUENCE  


	set.seed(rand.seed)
	n <- nrow(dagPeptides@data)

	if (length(matches) < n) stop("too few matches in background. Please try different parameters.", call. = FALSE)


	background <- lapply(seq_len(permutationSize), function(p) {
		s <- sample(matches, n, replace = replacement, prob = NULL)
		return(str_split_fixed(s, "", n = Inf))
	})

	new("dagBackground", background = background, permutationSize = permutationSize)
}

testDAU_AE				<- function (dagPeptides, dagBackground, group = c("null", "classic", "charge", "chemistry", "hydrophobicity")) {

	if (missing(dagPeptides) || class(dagPeptides) != "dagPeptides") {
		stop("dagPeptides should be an object of dagPeptides.\n\n\tPlease try ?fetchSequence to get help.", call. = FALSE)
	}
	if (missing(dagBackground) || class(dagBackground) != "dagBackground") {
		stop("dagBackground should be an object of dagBackground.\n\n\tPlease try ?buildBackgroundModel to get help.", call. = FALSE)
	}
	
	group   <- match.arg(group)
	exp     <- dagPeptides@peptides
	bg      <- dagBackground@background
	
	OneLet_AA   <- c("A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
	

	## Dataframes encoding what different AAs means in different contexts
	triLet_df   <- data.frame(
		Type = c("Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"),
		AA = OneLet_AA,
		stringsAsFactors = FALSE
	)
		
	classic_df	<- data.frame(
		Type = c(
			rep("nonpolar_aliphatic", 6),
			rep("polar_uncharged", 6),
			rep("aromatic", 3),
			rep("positively_charged", 3),
			rep("negatively_charged", 2)
			),
		AA = c("A", "G", "L", "M", "I", "V", "C", "P", "Q", "N", "S", "T", "F", "W", "Y", "H", "K", "R", "D", "E"),
		stringsAsFactors = FALSE
	)

	charge_df	<- data.frame(
		Type = c(
			rep("positive", 3),
			rep("neutral", 15),
			rep("negative", 2)
			),
		AA = c("H", "K", "R", "A", "C", "F", "G", "I", "L", "M", "N", "P", "Q", "S", "T", "V", "W", "Y", "D", "E"),
		stringsAsFactors = FALSE
	)

	chemist_df	<- data.frame(
		Type = c(
			rep("hydrophobic", 8),
			rep("polar", 5),
			rep("basic", 3),
			rep("neutral", 2),
			rep("acidic", 2)
			),
		AA = c("A", "F", "I", "L", "M", "P", "V", "W", "C", "G", "S", "T", "Y", "H", "K", "R", "N", "Q", "D", "E"),
		stringsAsFactors = FALSE
	)

	hydroph_df	<- data.frame(
		Type = c(
			rep("hydrophilic", 6),
			rep("neutral", 6),
			rep("hydrophobic", 8)
			),
		AA = c("D", "E", "K", "N", "Q", "R", "A", "G", "H", "P", "S", "T", "C", "F", "I", "L", "M", "V", "W", "Y"),
		stringsAsFactors = FALSE
	)
	
	typeList	<- list(null = triLet_df, classic = classic_df, charge = charge_df, chemistry = chemist_df, hydrophobicity = hydroph_df)

	
	exp_processed	<- apply(exp, 2, function(site) {
		tbl <- table(site)
		names(tbl)	<- typeList[[group]]$Type[match(names(tbl), typeList[[group]]$AA)]
		red_tbl		<- sapply(unique(typeList[[group]]$Type), function(type) sum(tbl[which(names(tbl) == type)]) / sum(tbl))
		return(red_tbl)
	})

	bg_processed	<- lapply(bg, function(replicate) {
		apply(replicate, 2, function(site) {
			tbl <- table(site)
			names(tbl)	<- typeList[[group]]$Type[match(names(tbl), typeList[[group]]$AA)]
			red_tbl		<- sapply(unique(typeList[[group]]$Type), function(type) sum(tbl[which(names(tbl) == type)]) / sum(tbl))
			return(red_tbl)
		})
	})
	
	if (ncol(exp_processed) != ncol(bg_processed[[1]])) {
		stop("the length of background is different from inputs", call. = FALSE)
	}
	

	bgPerSite		<- lapply(1:(ncol(exp_processed)), function(site) {
		siteCols	<- do.call(cbind, lapply(bg_processed, function(replicate) {
			replicate[, site]
		}))
	})

	bgPerSite_Std	<- do.call(cbind, lapply(bgPerSite, function(site) {
		apply(site, 1, sd, na.rm = TRUE)
	}))

	bgPerSite_mu	<- do.call(cbind, lapply(bgPerSite, function(site) {
		apply(site, 1, mean, na.rm = TRUE)
	}))



	exp_processed[is.na(exp_processed)]	<- 0
	perSite_diff	<- exp_processed - bgPerSite_mu
	perSite_diff[is.na(perSite_diff)]	<- 0

	perSite_zscore	<- perSite_diff / bgPerSite_Std

	siteNames		<- as.character(1:(ncol(exp_processed)))
	colnames(perSite_diff)	<- colnames(perSite_zscore)	<- siteNames

	perSite_pval	<- 2 * pnorm(-abs(perSite_zscore))

	new("testDAUresults", group = group, difference = perSite_diff, zscore = perSite_zscore, 
		pvalue = perSite_pval, background = bgPerSite_mu, motif = exp_processed, upstream = dagPeptides@upstreamOffset, 
		downstream = dagPeptides@downstreamOffset)
}

dagHeatmap_AE			<- function (testDAUresults, type = c("diff", "zscore"), ...) {

	if (missing(testDAUresults) || class(testDAUresults) != "testDAUresults") {
		stop("testDAUresults should be an object of testDAUresults\n\n\tPlease try ?testDAU to get help.", call. = FALSE)
	}

	type    <- match.arg(type)
	if (type == "diff") {
		dat <- testDAUresults@difference
	} else {
		dat <- testDAUresults@zscore
	}
	rownames(dat) <- str_to_title(str_replace(rownames(dat), "_", " "))
	dat     <- dat[order(rowSums(dat), decreasing = TRUE),]

	pheatmap::pheatmap(dat, ..., cluster_rows = FALSE, cluster_cols = FALSE, scale = "column")
}



## // Custom hydrophobicity scale to plug into featureHydro function // ##
featureHydro_AE	<- function (seq) {
	L = unique(sapply(seq, function(x) {
		length(unlist(strsplit(x, split = "")))
	}))
	
	H = c(0.31, -0.07, 0.56, 1.13, 0.24, 0.23, -0.17, -0.01, -0.14, -0.13, 1.85, 0.94, -0.45, -0.96, -2.02, -0.58, -1.23, -0.42, -0.99, -0.81)
	names(H) = c("I", "V", "L", "F", "C", "M", "A", "G", "T", "S", "W", "Y", "P", "H", "E", "Q", "D", "N", "K", "R")

	hydro = sapply(seq, function(x) {
		H[unlist(strsplit(x, split = ""))]
	})
	rownames(hydro) = paste("H:", 1:L, sep = "")

	return(t(hydro))
}


## // Codon usage analysis - the CAI function // ##
# The default function for calculating CAI is slow. We can make it faster by
# calculating some constant variables outside the function (i.e. we don't need
# to recalculate the same value for each sequence processed)
cai.faster  <- function (seq, w, exclude) {
	nncod <- uco(seq)
	nncod <- nncod[-exclude]
	sigma <- nncod %*% log(w)
	exp(sigma/sum(nncod))
}
# The vectorize function allows us to provide a vector or list as an argument,
# rather than looping over each sequence
cai.faster_vf <- Vectorize(cai.faster, vectorize.args = "seq", SIMPLIFY = TRUE)


## // Shine-Dalgarno // ##

# Given a P-set of interest and an approximate position of interest, e.g. 115 for "P1"
# plot a sequence logo spanning +/-30 bp of that position, e.g. to look for ATGs
# Example call: 
# P1_115_plot <-plotSDSeqLogo("P1", 115)
# P1_115_plot

plotSDSeqLogo	<- function(pset, position) {
	nucleot_seqs	<- SD_AllPosits_list[[pset]]$SD_contain_seqs
	subselect		<- subseq(nucleot_seqs, start = (position - 30), end = (position + 30))

	nucleo_chars	<- as.character(subselect)
	splitVector		<- split(1:60, ceiling(seq_along(1:60)/12))

	splitNucleoCassette <- lapply(splitVector, function(cassette) {
		cassette_range	<- c(min(cassette), max(cassette))
		sequence_split 	<- lapply(nucleo_chars, function(char_seq) {
			str_sub(char_seq, cassette_range[1], cassette_range[2])
		})
		return(unlist(sequence_split))
	})
	seqLogoPlot	<- ggseqlogo(data = splitNucleoCassette, method = "probability", ncol = 1, nrow = 5)
	return(seqLogoPlot)
}


## // RNAfold in-R processing // ##
readRNAfold	<- function (file) {

    if (!file.exists(file)) {
        stop(paste("The input file:", file, "does not exist."))
    }

    con <- file(file)
    lines <- readLines(con)
    close(con)


    if (class(lines) == "try-error" || length(lines) == 0) {
        return(data.frame(Structure = character(0), deltaG = numeric(0), stringsAsFactors = FALSE))
    }

    seq.idx <- seq(2, length(lines), 3)
    seqs <- lines[seq.idx]

    str.idx <- seq.idx + 1
    str <- lines[str.idx]

    str.s <- stringr::str_split(str, pattern = "  \\(  ")
    expr <- "([\\.\\(\\)|x><]+)[[:space:]]\\([[[:space:]]]*(.[^\\)^[[:space:]]]*)"

    m <- stringr::str_match(str, expr)
    delta.G <- as.numeric(m[, 3])
    structure <- m[, 2]
    result <- data.frame(Structure = structure, deltaG = delta.G, 
        stringsAsFactors = FALSE)
    return(result)
}





































