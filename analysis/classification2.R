library(tree)
decoy_feats = read.table('decoys/5pr_decoy_feats.txt', header=T, sep='\t')
decoy_feats = subset(decoy_feats, DeNovoDist > -400 & DeNovoDist < 400)
X = decoy_feats[,c("AuthenticMES", "DeNovoMES", "DeNovoDist", "InnerESE")]
X$y = X$DeNovoDist < 0
X$MES_diff = (X$DeNovoDist > 0) * (X$DeNovoMES-X$AuthenticMES) + (X$DeNovoDist < 0) * (X$AuthenticMES-X$DeNovoMES)
X$dist = abs(X$DeNovoDist)
X$ese_freq = X$InnerESE/X$dist
tr = tree(y ~ MES_diff + ese_freq, data=X)
plot(tr); text(tr, cex=0.8)

feats = read.table('dbass5/dbass5_de_novo_feats.txt', header=T, sep='\t')
X2 = feats[,c("AuthenticMES", "DeNovoMES", "DeNovoDist", "InnerESE")]
X2$y = X2$DeNovoDist > 0
X2$MES_diff = (X2$DeNovoDist > 0) * (X2$DeNovoMES-X2$AuthenticMES) + (X2$DeNovoDist < 0) * (X2$AuthenticMES-X2$DeNovoMES)
X2$dist = abs(X2$DeNovoDist)
X2$ese_freq = X2$InnerESE/X2$dist
y = as.data.frame(predict(tr, X2))
colnames(y) = 'conf'
y$yhat = yhat$conf > 0.5
y$y == X2$y

roc = as.data.frame(matrix(nrow=100,ncol=2))
colnames(roc) = c('sens', '1-spec')
for (i in seq(0.01,1,0.01)) {
  yhat = y$conf > i
  sens = sum(yhat*X2$y)/sum(X2$y)
  spec = sum((!yhat)*(!X2$y))/sum(!X2$y)
  roc[i*100,"sens"] = sens
  roc[i*100,"1-spec"] = 1-spec
}
plot(roc[,'1-spec'], roc[,'sens'], pch=19,xlab='1 - sensitivity', ylab='Specificity', bty='n')
abline(0,1,lty=2)

