setwd('~/Desktop/MacArthur/de-novo-splice/')

# get DBASS data
novo_feats = read.table('novo/data/5raw_feats.txt', header=T, sep='\t')
novo_feats$mes_diff = novo_feats$AltMES - novo_feats$RefMES
mean(subset(novo_feats, mes_diff > -15)$mes_diff) # 0.294
var(subset(novo_feats, mes_diff > -15)$mes_diff) # 10.790
cryptic_feats = read.table('cryptic/data/5raw_feats.txt', header=T, sep='\t')
cryptic_feats$mes_diff = cryptic_feats$AltMES - cryptic_feats$RefMES

# combine de novo and cryptic samples
alt = rbind(cryptic_feats, novo_feats)
mes_alt = unique(alt$mes_diff)
mean(mes_alt) # 3.532
var(mes_alt) # 49.58

# get decoys
decoy_feats = read.table('decoy/data/5raw_feats.txt', header=T, sep='\t')
decoy_feats$mes_diff = decoy_feats$AltMES - decoy_feats$RefMES

# downsample decoys so that training set is balanced
exonic = decoy_feats$AltLoc < decoy_feats$RefLoc
intronic = decoy_feats$AltLoc > decoy_feats$RefLoc
target = sum(exonic)
keep = sample(which(intronic), target)
balanced = rbind(decoy_feats[exonic,], decoy_feats[keep,])

# fit a gaussian
x = novo_feats$mes_diff[novo_feats$mes_diff > -8 & novo_feats$mes_diff < 10]
hist(x, breaks=20, prob=T, xlim=c(-9, 10.5))
# merge density curve bins
num_bins = 9
lines(density(x, n=num_bins), col='blue')
fit = density(x, n=num_bins)
x_lb = fit$x[-length(fit$x)]
x_ub = fit$x[-1]
y_lb = fit$y[-length(fit$y)]
y_ub = fit$y[-1]
# compute trapezoidal area of each bin, use as pdf function
bin_size = x_ub[1] - x_lb[1]
pdf = (y_lb + y_ub)/2 # compute relative area of each bin
breaks = fit$x
points(breaks[-1] - bin_size/2, pdf)

# get a feasible distribution to target in downsampling
y = balanced$mes_diff
binned = cut(y, breaks=breaks, labels=F)
counts = sapply(1:length(pdf), function(x) {sum(binned == x, na.rm=T)}) # get number of items in each bin
n_pot = counts/pdf # calculate how many total samples, assuming we max out one of the bins
pot = sapply(n_pot, function(n) { n*pdf }) # get corresponding distributions of samples
m = apply(pot, 2, function(col) { col <= counts })
j = which(apply(m, 2, all))
apply(as.matrix(pot[,j]), 2, function(col) {all(col > 1)}) # TRUE
colSums(as.matrix(pot[,j])) # > 1000
desired = ceiling(pot[,j]) # desired counts for each bin (i.e. target for downsampling)

# DOWNSAMPLE
score = 0
while (score < 0.47) {
  biased = data.frame()
  for (i in (1:length(pdf))) {
    idx = which(binned == i)
    keep = sample(idx, desired[i])
    biased = rbind(biased, balanced[keep,])
  }
  score = mean(biased$RefLoc > biased$AltLoc) 
  print(score)
}
mean(biased$RefLoc > biased$AltLoc)
hist(biased$mes_diff)
write.table(biased, file='decoy/data/5raw_feats.downsampled.txt', row.names=F, quote=F, sep='\t')

