# load in data and extract landmarks
setwd('~/Desktop/MacArthur/de-novo-splice/')
decoys = read.table('decoy/data/5raw_feats.txt', header=T, sep='\t')
upstream = apply(decoys[,c("RefLoc", "AltLoc")], 1, min)
downstream = apply(decoys[,c("RefLoc", "AltLoc")], 1, max)
exon_start = decoys$ExonStart
intron_end = decoys$IntronEnd
exonic = decoys$AltLoc < decoys$RefLoc

get_from_ref <- function(elem) {
	ragged = make_ragged(as.character(decoys[,elem]))
	from_ref = mapply(subtract, ragged, as.list(decoys$RefLoc))
	return(from_ref)
}

count <- function(v, lb, ub) { sum(v > lb & v <= ub)}

find_bp <- function(r, from_ref, ref_loc, hard_lb=exon_start+4, hard_ub=intron_end-50) {
	best = -Inf
	for (i in r) {
		lb_use = hard_lb - ref_loc
		ub_use = hard_ub - ref_loc
		n_exon = mapply(count, from_ref, lb_use, i) 
		n_intron = mapply(count, from_ref, i, ub_use)
		d_exon = i - lb_use - 1
		time_exon = sum(d_exon[d_exon > 0])
		d_intron = ub_use - i 
		time_intron = sum(d_intron[d_intron > 0])
		result = rateratio.test(c(n_exon, n_intron), c(time_exon, time_intron))
		score = result$p.value
		if (score > best) {
			best = score
			bp = i
		}
	}
	return(bp)
}

