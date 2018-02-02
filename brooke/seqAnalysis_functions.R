
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
		message(paste0("Made ", as.integer(subset*10000), " // ", as.integer(seqNumber), " random sequences"))
		return(seqSet)
	})
	combSets	<- Reduce(c, seqSubsets)
	names(combSets)	<- seq(1:length(combSets))
	return(combSets)
}

## // Update functions for the dagLogo package - mostly for speed // ##
formatSequence_AE	<- function (seq, proteome, upstreamOffset, downstreamOffset) {

    if (missing(proteome) || class(proteome) != "Proteome") {
        stop("proteome should be an object of Proteome. \n\n             Try ?prepareProteome to get help", 
            call. = FALSE)
    }
    if (missing(seq)) {
        stop("seq is required parameter.", call. = FALSE)
    }
    width <- unique(unlist(lapply(seq, nchar)))
    if (length(width) > 1) {
        stop("seq must be characters with same length", call. = FALSE)
    }
    if (missing(upstreamOffset) && missing(downstreamOffset)) {
        upstreamOffset <- floor(width/2)
        downstreamOffset <- width - upstreamOffset - 1
    }
    else {
        if (missing(downstreamOffset)) {
            downstreamOffset <- width - upstreamOffset - 1
        }
        else {
            upstreamOffset <- width - downstreamOffset - 1
        }
    }

    center <- upstreamOffset + 1

    ## Define the anchor AA
    anchorAA <- unlist(lapply(seq, function(.ele) substr(.ele, center, center)))
    anchorPos <- rep(1, length(seq))
    anchor  <- anchorAA

    dat <- data.frame(IDs = IDs, anchorAA = anchorAA, anchorPos = anchorPos, 
        peptide = seq, anchor = anchor)


    dat$upstream    <- str_sub(seq, 1, upstreamOffset)
    dat$downstream  <- str_sub(seq, upstreamOffset + 2, upstreamOffset + downstreamOffset + 1)


    seqchar.upstream <- do.call(rbind, lapply(dat$upstream, function(.seq) {
        .seq <- c(rep("NA", upstreamOffset), unlist(lapply(1:nchar(.seq), 
            function(i) substr(.seq, i, i))))
        .seq <- .seq[(length(.seq) - upstreamOffset + 1):length(.seq)]
        .seq
    }))


    seqchar.downstream <- do.call(rbind, lapply(dat$downstream, 
        function(.seq) {
            .seq <- c(unlist(lapply(1:nchar(.seq), function(i) substr(.seq, 
                i, i))), rep("NA", downstreamOffset))
            .seq <- .seq[1:downstreamOffset]
            .seq
    }))

    seqchar <- cbind(seqchar.upstream, as.character(dat$anchorAA), 
        seqchar.downstream)

    new("dagPeptides", data = dat, peptides = seqchar, upstreamOffset = upstreamOffset, 
        downstreamOffset = downstreamOffset, type = type)
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