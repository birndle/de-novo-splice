library(magrittr)
library(plyr)
library(dplyr)

# call on data frame as a whole
derive_feats_df <- function(df) {
  df %>%
    mutate(AltDist = AltLoc - RefLoc) %>%
    mutate(intronic = AltDist > 0) %>%
    mutate(ProximalDist = AltDist*intronic + -AltDist*(!intronic)) %>%
    mutate(DistalDist = 100) %>%
    mutate(ExonicDist = 100) %>%
    mutate(IntronicDist = 100) %>%
    mutate(ExonicLength = RefLoc - ExonStart) %>%
    mutate(IntronicLength = IntronEnd - RefLoc) -> new_df
  return(new_df)
}

# apply to each row of data frame
derive_feats_row <- function(slice) {
  ese_loc = strsplit(as.character(slice[,"ESELoc"]), ",")[[1]] %>% as.numeric
  ess_loc = strsplit(as.character(slice[,"ESSLoc"]), ",")[[1]] %>% as.numeric
  gtrip_loc = strsplit(as.character(slice[,"GtripLoc"]), ",")[[1]] %>% as.numeric
  ref_loc = slice[,'RefLoc']
  alt_loc = slice[,'AltLoc']
  exon_start = slice[,'ExonStart']
  intron_end = slice[,'IntronEnd']

  if (as.logical(slice[,"intronic"])) {
    slice <- transform(slice, 
    ProximalESE = sum((ese_loc > ref_loc) & (ese_loc < alt_loc)), 
    ProximalESS = sum((ess_loc > ref_loc) & (ess_loc < alt_loc)),
    ProximalGtrip = sum((gtrip_loc > ref_loc) & (gtrip_loc < alt_loc)),
    DistalESE = sum((ese_loc > alt_loc) & (ese_loc < (alt_loc + 100))),
    DistalESS = sum((ess_loc > alt_loc) & (ess_loc < (alt_loc + 100))),
    DistalGtrip = sum((gtrip_loc > alt_loc) & (gtrip_loc < (alt_loc + 100))))
  } else {
    slice <- transform(slice, 
    ProximalESE = sum((ese_loc > alt_loc) & (ese_loc < ref_loc)),
    ProximalESS = sum((ess_loc > alt_loc) & (ess_loc < ref_loc)),
    ProximalGtrip = sum((gtrip_loc > alt_loc) & (gtrip_loc < ref_loc)),
    DistalESE = sum((ese_loc > (alt_loc - 100)) & (ese_loc < alt_loc)),
    DistalESS = sum((ess_loc > (alt_loc - 100)) & (ess_loc < alt_loc)),
    DistalGtrip = sum((gtrip_loc > (alt_loc - 100)) & (gtrip_loc < alt_loc)))
  }
  transform(slice, 
  ExonicESE = sum((ese_loc > (ref_loc - 100)) & (ese_loc < ref_loc)),
  ExonicESS = sum((ess_loc > (ref_loc - 100)) & (ess_loc < ref_loc)),
  ExonicGtrip = sum((gtrip_loc > (ref_loc - 100)) & (gtrip_loc < ref_loc)),
  IntronicESE = sum((ese_loc > ref_loc) & (ese_loc < (ref_loc + 100))),
  IntronicESS = sum((ess_loc > ref_loc) & (ess_loc < (ref_loc + 100))),
  IntronicGtrip = sum((gtrip_loc > ref_loc) & (gtrip_loc < (ref_loc + 100))))
}

main <- function() {
  setwd('~/Desktop/MacArthur/de-novo-splice/')
  # just de novo events
  novo_feats = derive_feats_df(read.table('novo/data/de_novo5_feat_overlap.txt', header=T, sep='\t'))
  novo_feats = adply(novo_feats, 1, derive_feats_row)
  # all DBASS events
  dbass_feats = derive_feats_df(read.table('dbass5/data/dbass5_feat_overlap.txt', header=T, sep='\t'))
  dbass_feats = adply(dbass_feats, 1, derive_feats_row)
  rbind(dbass_feats, novo_feats) %>% duplicated(fromLast=T) -> dup_idx
  cryptic_feats = dbass_feats[!dup_idx[1:nrow(dbass_feats)],]
  # decoys
  decoy_feats = derive_feats_df(read.table('decoy/data/both/5pr_decoy_feat_overlap.txt', header=T, sep='\t'))
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