library(magrittr)
library(plyr)
library(dplyr)

options(run.main=F)
source(sprintf('feature_calc/feature_calc1%s.R', v_calc))
source(sprintf('feature_calc/feature_calc2%s.R', v_calc))

# just de novo events
novo_feats = derive_feats_df(read.table('../novo/data/de_novo5_feat.overlap.txt', header=T, sep='\t'))
novo_feats = adply(novo_feats, 1, derive_feats_row)

# all DBASS events
dbass_feats = derive_feats_df(read.table('../dbass5/data/dbass5_feat.overlap.txt', header=T, sep='\t'))
dbass_feats = adply(dbass_feats, 1, derive_feats_row)
rbind(dbass_feats, novo_feats) %>% duplicated(fromLast=T) -> dup_idx
cryptic_feats = dbass_feats[!dup_idx[1:nrow(dbass_feats)],]

# decoys
use = read.table(sprintf('../decoy/data/%s/5pr_decoy_feat.use%s.txt', vers, v_calc), header=T, sep='\t')
euse = subset(use, !intronic)
iuse = subset(use, intronic)

novo_feats = novo_feats[,colnames(use)]
cryptic_feats = cryptic_feats[,colnames(use)]

other = c('AltLoc', 'RefLoc', 'ExonStart', 'IntronEnd', 
          'AltMES', 'RefMES', 'AltDist', 'intronic',
          'ProximalDist', 'DistalDist','ExonicDist','IntronicDist', 
          'ExonicLength', 'IntronicLength')

# load in raw features
# decoys
ese_decoy = list()
ese_decoy$proximal = read.table(sprintf('../decoy/data/%s/cooked%s/ese_prox.txt', vers, v_calc), sep='\t', header=T)
ese_decoy$distal = read.table(sprintf('../decoy/data/%s/cooked%s/ese_dist.txt', vers, v_calc), sep='\t', header=T)
ese_decoy$exonic = read.table(sprintf('../decoy/data/%s/cooked%s/ese_exonic.txt', vers, v_calc), sep='\t', header=T)
ese_decoy$intronic = read.table(sprintf('../decoy/data/%s/cooked%s/ese_intronic.txt', vers, v_calc), sep='\t', header=T)

ess_decoy = list()
ess_decoy$proximal = read.table(sprintf('../decoy/data/%s/cooked%s/ess_prox.txt', vers, v_calc), sep='\t', header=T)
ess_decoy$distal = read.table(sprintf('../decoy/data/%s/cooked%s/ess_dist.txt', vers, v_calc), sep='\t', header=T)
ess_decoy$exonic = read.table(sprintf('../decoy/data/%s/cooked%s/ess_exonic.txt', vers, v_calc), sep='\t', header=T)
ess_decoy$intronic = read.table(sprintf('../decoy/data/%s/cooked%s/ess_intronic.txt', vers, v_calc), sep='\t', header=T)

gtrip_decoy = list()
gtrip_decoy$proximal = read.table(sprintf('../decoy/data/%s/cooked%s/gtrip_prox.txt', vers, v_calc), sep='\t', header=T)
gtrip_decoy$distal = read.table(sprintf('../decoy/data/%s/cooked%s/gtrip_dist.txt', vers, v_calc), sep='\t', header=T)
gtrip_decoy$intronic = read.table(sprintf('../decoy/data/%s/cooked%s/gtrip_intronic.txt', vers, v_calc), sep='\t', header=T)
gtrip_decoy$exonic = read.table(sprintf('../decoy/data/%s/cooked%s/gtrip_exonic.txt', vers, v_calc), sep='\t', header=T)

# de novo
# ese_novo = list()
# ese_novo$proximal = read.table('novo/data/cooked/ese_prox.txt', sep='\t', header=T)
# ese_novo$distal = read.table('novo/data/cooked/ese_dist.txt', sep='\t', header=T)
# ese_novo$exonic = read.table('novo/data/cooked/ese_exonic.txt', sep='\t', header=T)
# ese_novo$intronic = read.table('novo/data/cooked/ese_intronic.txt', sep='\t', header=T)
# 
# ess_novo = list()
# ess_novo$proximal = read.table('novo/data/cooked/ess_prox.txt', sep='\t', header=T)
# ess_novo$distal = read.table('novo/data/cooked/ess_dist.txt', sep='\t', header=T)
# ess_novo$exonic = read.table('novo/data/cooked/ess_exonic.txt', sep='\t', header=T)
# ess_novo$intronic = read.table('novo/data/cooked/ess_intronic.txt', sep='\t', header=T)
# 
# gtrip_novo = list()
# gtrip_novo$proximal = read.table('novo/data/cooked/gtrip_prox.txt', sep='\t', header=T)
# gtrip_novo$distal = read.table('novo/data/cooked/gtrip_dist.txt', sep='\t', header=T)
# gtrip_novo$intronic = read.table('novo/data/cooked/gtrip_intronic.txt', sep='\t', header=T)
# gtrip_novo$exonic = read.table('novo/data/cooked/gtrip_exonic.txt', sep='\t', header=T)

# cryptic
# ese_cryp = list()
# ese_cryp$proximal = read.table('cryptic/data/cooked/ese_prox.txt', sep='\t', header=T)
# ese_cryp$distal = read.table('cryptic/data/cooked/ese_dist.txt', sep='\t', header=T)
# ese_cryp$exonic = read.table('cryptic/data/cooked/ese_exonic.txt', sep='\t', header=T)
# ese_cryp$intronic = read.table('cryptic/data/cooked/ese_intronic.txt', sep='\t', header=T)
# 
# ess_cryp = list()
# ess_cryp$proximal = read.table('cryptic/data/cooked/ess_prox.txt', sep='\t', header=T)
# ess_cryp$distal = read.table('cryptic/data/cooked/ess_dist.txt', sep='\t', header=T)
# ess_cryp$exonic = read.table('cryptic/data/cooked/ess_exonic.txt', sep='\t', header=T)
# ess_cryp$intronic = read.table('cryptic/data/cooked/ess_intronic.txt', sep='\t', header=T)
# 
# gtrip_cryp = list()
# gtrip_cryp$proximal = read.table('cryptic/data/cooked/gtrip_prox.txt', sep='\t', header=T)
# gtrip_cryp$distal = read.table('cryptic/data/cooked/gtrip_dist.txt', sep='\t', header=T)
# gtrip_cryp$intronic = read.table('cryptic/data/cooked/gtrip_intronic.txt', sep='\t', header=T)
# gtrip_cryp$exonic = read.table('cryptic/data/cooked/gtrip_exonic.txt', sep='\t', header=T)
