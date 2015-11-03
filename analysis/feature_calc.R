library(magrittr)
library(dplyr)

setwd('~/Desktop/MacArthur/DBASS/analysis')
feats = read.table('../novo/data/dbass5_de_novo_feats.txt', header=T, sep='\t')
feats$intronic = feats$DeNovoDist > 0
decoy_feats = read.table('../decoy/5pr_decoy_feats.txt', header=T, sep='\t')
decoy_feats$intronic = decoy_feats$DeNovoDist > 0

s = subset(decoy_feats, (DeNovoDist > -100 & DeNovoDist < -15) | (DeNovoDist < 100 & DeNovoDist > 15))
use = subset(s, DeNovoMES > 0)
iuse = subset(use, intronic)
euse = subset(use, !intronic)

sre_cols = c("InnerESE", "InnerESS", "InnerGtrip", "OuterESE", "OuterESS", 
             "OuterGtrip", "ExonicESE", "ExonicESS", "ExonicGtrip", 
             "IntronicESE", "IntronicESS", "IntronicGtrip")

ci <- function(k, n) {
  p = k/n
  margin = 1.96*p*(1-p)/sqrt(n)
  return(c(p, margin))
}

pois_ci <- function(k, n) {
  p = k/n
  se = sqrt(p/n)
  margin = 1.96*se
  return(c(p, margin))
}
