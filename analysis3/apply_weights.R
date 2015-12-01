# add columns with weighted sum of ESE/ESS elements
source('weight_hexamers.R')

# -- de novo --
# ESE
novo_feats[,"ProximalESEWeighted"] = as.matrix(ese_novo$proximal[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
novo_feats[,"DistalESEWeighted"] = as.matrix(ese_novo$distal[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
novo_feats[,"ExonicESEWeighted"] = as.matrix(ese_novo$exonic[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
novo_feats[,"IntronicESEWeighted"] = as.matrix(ese_novo$intronic[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
# ESS
novo_feats[,"ProximalESSWeighted"] = as.matrix(ess_novo$proximal[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector
novo_feats[,"DistalESSWeighted"] = as.matrix(ess_novo$distal[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector
novo_feats[,"ExonicESSWeighted"] = as.matrix(ess_novo$exonic[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector
novo_feats[,"IntronicESSWeighted"] = as.matrix(ess_novo$intronic[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector

# -- cryptic --
# ESE
cryptic_feats[,"ProximalESEWeighted"] = as.matrix(ese_cryp$proximal[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
cryptic_feats[,"DistalESEWeighted"] = as.matrix(ese_cryp$distal[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
cryptic_feats[,"ExonicESEWeighted"] = as.matrix(ese_cryp$exonic[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
cryptic_feats[,"IntronicESEWeighted"] = as.matrix(ese_cryp$intronic[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
# ESS
cryptic_feats[,"ProximalESSWeighted"] = as.matrix(ess_cryp$proximal[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector
cryptic_feats[,"DistalESSWeighted"] = as.matrix(ess_cryp$distal[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector
cryptic_feats[,"ExonicESSWeighted"] = as.matrix(ess_cryp$exonic[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector
cryptic_feats[,"IntronicESSWeighted"] = as.matrix(ess_cryp$intronic[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector
# sort cryptic columns to match de novo columns
cryptic_feats = cryptic_feats[,colnames(novo_feats)]

# -- decoy --
# ESE
use[,"ProximalESEWeighted"] = as.matrix(ese_decoy$proximal[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
use[,"DistalESEWeighted"] = as.matrix(ese_decoy$distal[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
use[,"ExonicESEWeighted"] = as.matrix(ese_decoy$exonic[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
use[,"IntronicESEWeighted"] = as.matrix(ese_decoy$intronic[,ese_cols]) %*% as.matrix(ese_w) %>% as.vector
# ESS
use[,"ProximalESSWeighted"] = as.matrix(ess_decoy$proximal[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector
use[,"DistalESSWeighted"] = as.matrix(ess_decoy$distal[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector
use[,"ExonicESSWeighted"] = as.matrix(ess_decoy$exonic[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector
use[,"IntronicESSWeighted"] = as.matrix(ess_decoy$intronic[,ess_cols]) %*% as.matrix(ess_w) %>% as.vector

