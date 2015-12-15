# load in data and extract landmarks
setwd('~/Desktop/MacArthur/de-novo-splice/')
decoys = read.table('decoy/data/5raw_feats.txt', header=T, sep='\t')
upstream = apply(decoys[,c("RefLoc", "AltLoc")], 1, min)
downstream = apply(decoys[,c("RefLoc", "AltLoc")], 1, max)
exon_start = decoys$ExonStart
intron_end = decoys$IntronEnd
exonic = decoys$AltLoc < decoys$RefLoc

get_matrix <- function(elem) {
  ragged = make_ragged(as.character(decoys[,elem]))
  from_up = mapply(subtract, ragged, as.list(upstream))
  from_down = mapply(subtract, ragged, as.list(downstream))
  f(elem, 10:-40, from_up) # upstream
  g1(elem, -10:40, from_up) # middle up
  g2(elem, 10:-150, from_down) # middle down
  h(elem, -10:40, from_down) # downstream
}

# compute distance of ESE/ESS elements to landmarks
get_matrix('RescueESE')
get_matrix('FAS.ESS')
get_matrix('Gtrip')
get_matrix('PESE')
get_matrix('PESS')

# ESEfinder motifs
get_matrix('SRp40')
get_matrix('SF2_ASF')
get_matrix('SC35')
get_matrix('SRp55')

# compute zeros of the first derivative. if more than one, take the better one
