# write raw features
# --- de novo ---
# ESE
write.table(ese_novo$proximal, file='novo/data/cooked/ese_prox.txt', quote=F, sep="\t")
write.table(ese_novo$distal, file='novo/data/cooked/ese_dist.txt', quote=F, sep="\t")
write.table(ese_novo$exonic, file='novo/data/cooked/ese_exonic.txt', quote=F, sep="\t")
write.table(ese_novo$intronic, file='novo/data/cooked/ese_intronic.txt', quote=F, sep="\t")
# ESS
write.table(ess_novo$intronic, file='novo/data/cooked/ess_intronic.txt', quote=F, sep="\t")
write.table(ess_novo$exonic, file='novo/data/cooked/ess_exonic.txt', quote=F, sep="\t")
write.table(ess_novo$distal, file='novo/data/cooked/ess_dist.txt', quote=F, sep="\t")
write.table(ess_novo$proximal, file='novo/data/cooked/ess_prox.txt', quote=F, sep="\t")
# G-triplets
write.table(gtrip_novo$proximal, file='novo/data/cooked/gtrip_prox.txt', quote=F, sep="\t")
write.table(gtrip_novo$distal, file='novo/data/cooked/gtrip_dist.txt', quote=F, sep="\t")
write.table(gtrip_novo$intronic, file='novo/data/cooked/gtrip_intronic.txt', quote=F, sep="\t")
write.table(gtrip_novo$exonic, file='novo/data/cooked/gtrip_exonic.txt', quote=F, sep="\t")
# --- cryptic ---
# ESE
write.table(ese_cryp$proximal, file='cryptic/data/cooked/ese_prox.txt', quote=F, sep="\t")
write.table(ese_cryp$distal, file='cryptic/data/cooked/ese_dist.txt', quote=F, sep="\t")
write.table(ese_cryp$exonic, file='cryptic/data/cooked/ese_exonic.txt', quote=F, sep="\t")
write.table(ese_cryp$intronic, file='cryptic/data/cooked/ese_intronic.txt', quote=F, sep="\t")
# ESS
write.table(ess_cryp$intronic, file='cryptic/data/cooked/ess_intronic.txt', quote=F, sep="\t")
write.table(ess_cryp$exonic, file='cryptic/data/cooked/ess_exonic.txt', quote=F, sep="\t")
write.table(ess_cryp$distal, file='cryptic/data/cooked/ess_dist.txt', quote=F, sep="\t")
write.table(ess_cryp$proximal, file='cryptic/data/cooked/ess_prox.txt', quote=F, sep="\t")
# G-triplets
write.table(gtrip_cryp$proximal, file='cryptic/data/cooked/gtrip_prox.txt', quote=F, sep="\t")
write.table(gtrip_cryp$distal, file='cryptic/data/cooked/gtrip_dist.txt', quote=F, sep="\t")
write.table(gtrip_cryp$intronic, file='cryptic/data/cooked/gtrip_intronic.txt', quote=F, sep="\t")
write.table(gtrip_cryp$exonic, file='cryptic/data/cooked/gtrip_exonic.txt', quote=F, sep="\t")

# --- decoy ---
# write aggregated decoy features
feat_file = 'both/5pr_decoy_feat.useC.txt'
write.table(use, file=sprintf('decoy/data/%s', feat_file), quote=F, sep='\t')
# ESE
vers = 'both/cookedC'
write.table(ese_decoy$proximal, file=sprintf('decoy/data/%s/ese_prox.txt', vers), quote=F, sep="\t")
write.table(ese_decoy$distal, file=sprintf('decoy/data/%s/ese_dist.txt', vers), quote=F, sep="\t")
write.table(ese_decoy$intronic, file=sprintf('decoy/data/%s/ese_intronic.txt', vers), quote=F, sep="\t")
write.table(ese_decoy$exonic, file=sprintf('decoy/data/%s/ese_exonic.txt', vers), quote=F, sep="\t")
# ESS
write.table(ess_decoy$distal, file=sprintf('decoy/data/%s/ess_dist.txt', vers), quote=F, sep="\t")
write.table(ess_decoy$exonic, file=sprintf('decoy/data/%s/ess_exonic.txt', vers), quote=F, sep="\t")
write.table(ess_decoy$proximal, file=sprintf('decoy/data/%s/ess_prox.txt', vers), quote=F, sep="\t")
write.table(ess_decoy$intronic, file=sprintf('decoy/data/%s/ess_intronic.txt', vers), quote=F, sep="\t")
# G-triplets
write.table(gtrip_decoy$proximal, file=sprintf('decoy/data/%s/gtrip_prox.txt', vers), quote=F, sep="\t")
write.table(gtrip_decoy$distal, file=sprintf('decoy/data/%s/gtrip_dist.txt', vers), quote=F, sep="\t")
write.table(gtrip_decoy$intronic, file=sprintf('decoy/data/%s/gtrip_intronic.txt', vers), quote=F, sep="\t")
write.table(gtrip_decoy$exonic, file=sprintf('decoy/data/%s/gtrip_exonic.txt', vers), quote=F, sep="\t")
