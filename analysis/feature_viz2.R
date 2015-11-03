# Distribution of decoys around authentic sites + enrichment of inner ESS w.r.t DecoyDist

source('~/Desktop/MacArthur/DBASS/analysis/feature_calc.R')

# distribution of decoy sites w.r.t authentic sites
hist(decoy_feats$DeNovoDist, main='Distribution of decoys around authentic')
hist(use$DeNovoDist, breaks=100, main='Distribution of decoys around authentic')

# distributions of intron and exon lengths
intron_lengths = s$IntronEnd - s$AuthenticLocation
s$IntronLength = intron_lengths
hist(intron_lengths[intron_lengths < 12000], breaks=100, main='Distribution of intron lengths')
exon_lengths = s$AuthenticLocation - s$ExonStart
hist(exon_lengths[exon_lengths < 1000], breaks=100, main='Distribution of exon lengths')
s$ExonLength = exon_lengths

# Enrichment of inner ESS for intronic decoys w.r.t DecoyDist
bg = sum(use$OuterESS) / sum(use$IntronEnd - use$DeNovoLocation)
b = seq(0, 200, 10)
bins = cut(use$DeNovoDist, breaks=b, labels=F)
x = 1:max(bins, na.rm=T)
y = vector(length=length(x), mode='numeric') # ess density in intron / ess density in exon
for (i in x) {
  d = use[bins == i & !is.na(bins), c("DeNovoDist", "InnerESS")] %>% colSums
  d = d["InnerESS"]/d["DeNovoDist"]
  y[i] = d/bg
}
plot(x,y, xaxt='n', bty='n', main='Inner vs. Outer ESS enrichment for intronic decoys',
     xlab='Distance between authentic and decoy',
     ylab='Fold difference in ESS density')
abline(h=1, lty=2)
axis(1, at=x, labels=b[1:(length(b)-1)])
