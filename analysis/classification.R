library("party")
feats = read.table('dbass5/dbass5_de_novo_feats.txt', header=T, sep='\t')
feats$intronic = feats$DeNovoDist > 0
feats$denovo = rep(TRUE,nrow(feats))
decoy_feats = read.table('decoys/5pr_decoy_feats.txt', header=T, sep='\t')
clipped = subset(decoy_feats, DeNovoDist > -400 & DeNovoDist < 400)
clipped$intronic = clipped$DeNovoDist > 0
clipped$denovo = rep(FALSE,nrow(clipped))
m = rbind(clipped, feats)

# Inner G triplets in intron
intronic = subset(m, intronic)
intronic$MES_diff = intronic$DeNovoMES-intronic$AuthenticMES
# intronic$Gtripfreq = intronic$InnerGtrip/(intronic$DeNovoLocation - intronic$AuthenticLocation)
intronic$Gtripfreq = intronic$InnerGtrip/(intronic$AuthenticLocation - intronic$ExonStart)
tr1 <- ctree(denovo ~ MES_diff + DeNovoDist + Gtripfreq, data=intronic)

# Inner ESS in exon
exonic = subset(m, !intronic)
exonic$MES_diff = exonic$DeNovoMES-exonic$AuthenticMES
exonic$ESSfreq = exonic$InnerESS/(exonic$AuthenticLocation - exonic$DeNovoLocation)
tr2 = ctree(denovo ~ MES_diff + DeNovoDist + ESSfreq, data=exonic)
