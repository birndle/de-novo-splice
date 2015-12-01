library(magrittr)
library(dplyr)

get_proximal_feats <- function(slice) {
  count_inner <- function(cell, intronic) {
    if (is.na(cell)) {
      return(0)
    } else {
      loc = strsplit(as.character(cell), ",")[[1]] %>% as.numeric
      if (intronic) {
        return(sum((loc > slice$RefLoc) & (loc < slice$AltLoc)))
      } else {
        return(sum((loc > slice$AltLoc) & (loc < slice$RefLoc)))
      }
    }
  }
  feats = slice[,!(colnames(slice) %in% other)]
  sapply(feats, count_inner, slice$intronic)
}

get_distal_feats <- function(slice) {
  count_outer <- function(cell, intronic) {
    if (is.na(cell)) {
      return(0)
    } else {
      ialt = min(slice$IntronEnd, slice$AltLoc + 100)
      ealt = max(slice$ExonStart, slice$AltLoc - 100)
      loc = strsplit(as.character(cell), ",")[[1]] %>% as.numeric
      if (intronic) {
        return(sum((loc > slice$AltLoc) & (loc < ialt)))
      } else {
        return(sum((loc > ealt) & (loc < slice$AltLoc)))
      }
    }
  }
  feats = slice[,!(colnames(slice) %in% other)]
  sapply(feats, count_outer, slice$intronic)
}

get_exonic_feats <- function(slice) {
  count_exonic <- function(cell, intronic) {
    if (is.na(cell)) {
      return(0)
    } else {
      eref = max(slice$ExonStart, slice$RefLoc - 100)
      loc = strsplit(as.character(cell), ",")[[1]] %>% as.numeric
      return(sum((loc > eref) & (loc < slice$RefLoc)))
    }
  }
  feats = slice[,!(colnames(slice) %in% other)]
  sapply(feats, count_exonic)
}

get_intronic_feats <- function(slice, k) {
  count_intronic <- function(cell) {
    if (is.na(cell)) {
      return(0)
    } else {
      iref = min(slice$IntronEnd, slice$RefLoc + 100)
      loc = strsplit(as.character(cell), ",")[[1]] %>% as.numeric
      return(sum((loc > slice$RefLoc) & (loc < iref)))
    }
  }
  feats = slice[,!(colnames(slice) %in% other)]
  sapply(feats, count_intronic)
}

other = c('AltLoc', 'RefLoc', 'ExonStart', 'IntronEnd', 
          'AltMES', 'RefMES', 'AltDist', 'intronic',
          'ProximalDist', 'DistalDist','ExonicDist','IntronicDist', 
          'ExonicLength', 'IntronicLength')

# --- ESE ---
# de novo
ese_raw_novo = derive_feats_df(read.table('../novo/data/raw/dbass5_ese_feats.txt', header=T))
ese_novo = list()
ese_novo$proximal = adply(ese_raw_novo, 1, get_proximal_feats)
ese_novo$distal = adply(ese_raw_novo, 1, get_distal_feats)
ese_novo$exonic = adply(ese_raw_novo, 1, get_exonic_feats)
ese_novo$intronic = adply(ese_raw_novo, 1, get_intronic_feats)
# cryptic
ese_raw_cryp = derive_feats_df(read.table('../cryptic/data/raw/cryptic5_ese_feats.txt', header=T))
ese_cryp = list()
ese_cryp$proximal = adply(ese_raw_cryp, 1, get_proximal_feats)
ese_cryp$distal = adply(ese_raw_cryp, 1, get_distal_feats)
ese_cryp$exonic = adply(ese_raw_cryp, 1, get_exonic_feats)
ese_cryp$intronic = adply(ese_raw_cryp, 1, get_intronic_feats)
# decoys
# ese_raw_decoy = derive_feats_df(read.table('decoy/data/both/raw/5pr_ese_feats.txt', header=T))
# # get index object for use subset of decoys, general to all feature data frames
# tmp = ese_raw_decoy
# use_idx = tmp$AltDist > -310 & tmp$AltDist < 197
# ese_raw_decoy = ese_raw_decoy[use_idx,]
# 
# ese_decoy = list()
# ese_decoy$proximal = adply(ese_raw_decoy, 1, get_proximal_feats)
# ese_decoy$distal = adply(ese_raw_decoy, 1, get_distal_feats)
# ese_decoy$exonic = adply(ese_raw_decoy, 1, get_exonic_feats)
# ese_decoy$intronic = adply(ese_raw_decoy, 1, get_intronic_feats)

# --- ESS ---
# de novo
ess_raw_novo = derive_feats_df(read.table('../novo/data/raw/dbass5_ess_feats.txt', header=T))
ess_novo = list()
ess_novo$proximal = adply(ess_raw_novo, 1, get_proximal_feats)
ess_novo$distal = adply(ess_raw_novo, 1, get_distal_feats)
ess_novo$exonic = adply(ess_raw_novo, 1, get_exonic_feats)
ess_novo$intronic = adply(ess_raw_novo, 1, get_intronic_feats)
# cryptic
ess_raw_cryp = derive_feats_df(read.table('../cryptic/data/raw/cryptic5_ess_feats.txt', header=T))
ess_cryp = list()
ess_cryp$proximal = adply(ess_raw_cryp, 1, get_proximal_feats)
ess_cryp$distal = adply(ess_raw_cryp, 1, get_distal_feats)
ess_cryp$exonic = adply(ess_raw_cryp, 1, get_exonic_feats)
ess_cryp$intronic = adply(ess_raw_cryp, 1, get_intronic_feats)
# decoy
# ess_raw_decoy = derive_feats_df(read.table('decoy/data/both/raw/5pr_ess_feats.txt', header=T))
# ess_raw_decoy = ess_raw_decoy[use_idx,]
# ess_decoy = list()
# ess_decoy$proximal = adply(ess_raw_decoy, 1, get_proximal_feats)
# ess_decoy$distal = adply(ess_raw_decoy, 1, get_distal_feats)
# ess_decoy$exonic = adply(ess_raw_decoy, 1, get_exonic_feats)
# ess_decoy$intronic = adply(ess_raw_decoy, 1, get_intronic_feats)

# --- G-triplets ---
# de novo
gtrip_raw_novo = derive_feats_df(read.table('../novo/data/raw/dbass5_gtrip_feats.txt', header=T, sep="\t"))
gtrip_novo = list()
gtrip_novo$proximal = adply(gtrip_raw_novo, 1, get_proximal_feats)
gtrip_novo$distal = adply(gtrip_raw_novo, 1, get_distal_feats)
gtrip_novo$exonic = adply(gtrip_raw_novo, 1, get_exonic_feats)
gtrip_novo$intronic = adply(gtrip_raw_novo, 1, get_intronic_feats)
# cryptic
gtrip_raw_cryp = derive_feats_df(read.table('../cryptic/data/raw/cryptic5_gtrip_feats.txt', header=T, sep="\t"))
gtrip_cryp = list()
gtrip_cryp$proximal = adply(gtrip_raw_cryp, 1, get_proximal_feats)
gtrip_cryp$distal = adply(gtrip_raw_cryp, 1, get_distal_feats)
gtrip_cryp$exonic = adply(gtrip_raw_cryp, 1, get_exonic_feats)
gtrip_cryp$intronic = adply(gtrip_raw_cryp, 1, get_intronic_feats)
# decoy
# gtrip_raw_decoy = derive_feats_df(read.table('decoy/data/both/raw/5pr_gtrip_feats.txt', header=T, sep="\t"))
# gtrip_raw_decoy = gtrip_raw_decoy[use_idx,]
# gtrip_decoy = list()
# gtrip_decoy$proximal = adply(gtrip_raw_decoy, 1, get_proximal_feats)
# gtrip_decoy$distal = adply(gtrip_raw_decoy, 1, get_distal_feats)
# gtrip_decoy$exonic = adply(gtrip_raw_decoy, 1, get_exonic_feats)
# gtrip_decoy$intronic = adply(gtrip_raw_decoy, 1, get_intronic_feats)
