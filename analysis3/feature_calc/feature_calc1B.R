library(magrittr)
library(plyr)
library(dplyr)

# call on data frame as a whole
derive_feats_df <- function(df) {
  eref = apply(cbind(df$ExonStart + 3, df$RefLoc - 100), 1, max)
  ealt = apply(cbind(df$ExonStart + 3, df$AltLoc - 100), 1, max)
  iref = apply(cbind(df$IntronEnd - 20, df$RefLoc + 100), 1, min)
  ialt = apply(cbind(df$IntronEnd - 20, df$AltLoc + 100), 1, min)
  df %>%
    mutate(AltDist = AltLoc - RefLoc) %>%
    mutate(intronic = AltDist > 0) %>%
    mutate(ProximalDist = (AltDist*intronic + -AltDist*(!intronic)) - 9) %>%
    mutate(DistalDist = (ialt - AltLoc - 6)*intronic + (AltLoc - ealt - 3)*(!intronic)) %>%
    mutate(ExonicDist = RefLoc - eref - 3) %>%
    mutate(IntronicDist = iref - RefLoc - 6) %>%
    mutate(ExonicLength = RefLoc - ExonStart) %>%
    mutate(IntronicLength = IntronEnd - RefLoc) -> new_df
  return(new_df)
}

# apply to each row of data frame
derive_feats_row <- function(slice) {
  ese_loc = strsplit(as.character(slice[,"ESELoc"]), ",")[[1]] %>% as.numeric
  ess_loc = strsplit(as.character(slice[,"ESSLoc"]), ",")[[1]] %>% as.numeric
  gtrip_loc = strsplit(as.character(slice[,"GtripLoc"]), ",")[[1]] %>% as.numeric
  # gtrip_len = strsplit(as.character(slice[,"GtripLen"]), ",")[[1]] %>% as.numeric
  ref_loc = slice[,'RefLoc']
  alt_loc = slice[,'AltLoc']
  exon_start = slice[,'ExonStart']
  intron_end = slice[,'IntronEnd']
  # define bounds for where to look for SREs
  eref = max(exon_start + 3, ref_loc - 100)
  ealt = max(exon_start + 3, alt_loc - 100)
  iref = min(intron_end - 20, ref_loc + 100)
  ialt = min(intron_end - 20, alt_loc + 100)
  
  if (as.logical(slice[,"intronic"])) {
    slice <- transform(slice,
                       ProximalESE = sum((ese_loc > ref_loc + 6) & (ese_loc + 6 < alt_loc - 3)),
                       ProximalESS = sum((ess_loc > ref_loc + 6) & (ess_loc + 6 < alt_loc - 3)),
                       ProximalGtrip = sum((gtrip_loc > ref_loc + 6) & (gtrip_loc + 3 < alt_loc - 3)),
                       DistalESE = sum((ese_loc > alt_loc + 6) & (ese_loc + 6 < ialt)),
                       DistalESS = sum((ess_loc > alt_loc + 6) & (ess_loc + 6 < ialt)),
                       DistalGtrip = sum((gtrip_loc > alt_loc + 6) & (gtrip_loc + 3 < ialt)))
  } else {
    slice <- transform(slice,
                       ProximalESE = sum((ese_loc > alt_loc + 6) & (ese_loc + 6 < ref_loc - 3)),
                       ProximalESS = sum((ess_loc > alt_loc + 6) & (ess_loc + 6 < ref_loc - 3)),
                       ProximalGtrip = sum((gtrip_loc > alt_loc + 6) & (gtrip_loc + 3 < ref_loc - 3)),
                       DistalESE = sum((ese_loc > ealt) & (ese_loc + 6 < alt_loc - 3)),
                       DistalESS = sum((ess_loc > ealt) & (ess_loc + 6 < alt_loc - 3)),
                       DistalGtrip = sum((gtrip_loc > ealt) & (gtrip_loc + 3 < alt_loc - 3)))
  }
  transform(slice,
            ExonicESE = sum((ese_loc > eref) & (ese_loc + 6 < ref_loc - 3)),
            ExonicESS = sum((ess_loc > eref) & (ess_loc + 6 < ref_loc - 3)),
            ExonicGtrip = sum((gtrip_loc > eref) & (gtrip_loc + 3 < ref_loc - 3)),
            IntronicESE = sum((ese_loc > ref_loc + 6) & (ese_loc + 6 < iref)),
            IntronicESS = sum((ess_loc > ref_loc + 6) & (ess_loc + 6 < iref)),
            IntronicGtrip = sum((gtrip_loc > ref_loc + 6) & (gtrip_loc + 3 < iref)))
}

main <- function() {
  setwd('~/Desktop/MacArthur/de-novo-splice/')
  # just de novo events
  novo_feats = derive_feats_df(read.table('novo/data/de_novo5_feat.overlap.txt', header=T, sep='\t'))
  novo_feats = adply(novo_feats, 1, derive_feats_row)
  # all DBASS events
  dbass_feats = derive_feats_df(read.table('dbass5/data/dbass5_feat.overlap.txt', header=T, sep='\t'))
  dbass_feats = adply(dbass_feats, 1, derive_feats_row)
  rbind(dbass_feats, novo_feats) %>% duplicated(fromLast=T) -> dup_idx
  cryptic_feats = dbass_feats[!dup_idx[1:nrow(dbass_feats)],]
  # decoys
  decoy_feats = derive_feats_df(read.table('decoy/data/5pr_decoy_feat.overlap.txt', header=T, sep='\t'))
  use = subset(decoy_feats, AltDist > -310 & AltDist < 197)
  tmp = decoy_feats
  use_idx = tmp$AltDist > -310 & tmp$AltDist < 197
  use = adply(use, 1, derive_feats_row)
  iuse = subset(use, intronic)
  euse = subset(use, !intronic) 
}

if (getOption('run.main', default=TRUE)) {
  main()
}