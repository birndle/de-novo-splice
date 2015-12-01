# source('~/Desktop/MacArthur/DBASS/analysis/feature_calc.R')

smooth <- function(df, x, xlab, b=seq(0, 0.3, 0.05)) {
  bins = cut(x, breaks=b, labels=F)
  y = vector(mode='numeric', length=max(bins, na.rm=T))
  for (i in 1:max(bins, na.rm=T)) {
    bin = df[bins==i & !is.na(bins),]
    y[i] = median(bin$AltMES - bin$RefMES)
  }
  x = 1:max(bins, na.rm=T)
  plot(x, y, ylab="Authentic - Decoy MES", xlab=xlab)
  cor.test(x, y)
}

# --- ESE compensation for exonic decoys ---
# inner ESE density vs. relative ss strength (i.e. authentic - decoy MES)
tmp = subset(use, !intronic & (AltMES - RefMES) > -10 & ProximalDist > 0)
ese_dens = tmp$ProximalESE
y = tmp$AltMES - tmp$RefMES
# plot relationship
plot(ese_dens, y, ylab='Decoy - Authentic MES',
     main='ESE compensating for exonic decoy') # f(ese_dens) = mes_diff
cor.test(ese_dens, y)
# fit regression line
l = lm(y ~ ese_dens)
abline(l, col='red', lty=2)
stats = cor.test(ese_dens, y)
margin = (stats$conf.int[2] - stats$conf.int[1])/2
legend("topright", cex=0.9, bty='n',
       sprintf("r = %.3f +/- %.3f", stats$estimate, margin))

# learn individual feature weights by running multiple linear regression
m = subset(ese_decoy$proximal, !intronic & (AltMES - RefMES) > -10 & ProximalDist > 0)
y_true = m$AltMES - m$RefMES
ese_cols = colnames(m)[!(colnames(m) %in% other)]
X = m[,ese_cols]
fit = lm(y_true ~ ., data=X)
ese_w = fit$coefficients[2:length(fit$coefficients)]
ese_w[is.na(ese_w)] = 0

# plot relationship with weighted sum
x_ese = subset(ese_decoy$proximal, !intronic & (AltMES - RefMES) > -10 & ProximalDist > 0, select=ese_cols) %>% 
  as.matrix() %*% as.matrix(ese_w)
plot(x_ese, y)
cor.test(x_ese,y)
l = lm(y ~ x_ese)
abline(l, col='red', lty=2)

# --- ESS compensation for intronic decoys ---
# inner ESS density vs. relative ss strength (i.e. authentic - decoy MES)
tmp = subset(use, intronic & (AltMES - RefMES) > -10 & ProximalDist > 0)
ess_dens = tmp$ProximalESS
y = tmp$AltMES - tmp$RefMES
plot(ess_dens, y, ylab='Decoy - Authentic MES',
     main='ESS compensating for intronic decoys') 
cor.test(ess_dens, y)

# fit linear regression line
l = lm(y ~ ess_dens)
abline(l, col='red', lty=2)
stats = cor.test(ess_dens, y)
margin = (stats$conf.int[2] - stats$conf.int[1])/2
legend("topright", cex=0.9, bty='n',
       sprintf("r = %.3f +/- %.3f", stats$estimate, margin))

# learn individual feature weights by running multiple linear regression
m = subset(ess_decoy$proximal, intronic & (AltMES - RefMES) > -10 & ProximalDist > 0)
y_true = m$AltMES - m$RefMES
ess_cols = colnames(m)[!(colnames(m) %in% other)]
X = m[,ess_cols]
fit = lm(y_true ~ ., data=X)
ess_w = fit$coefficients[2:length(fit$coefficients)]

# plot relationship with weighted sum
x_ess = subset(ess_decoy$proximal, intronic & (AltMES - RefMES) > -10 & ProximalDist > 0, select=ess_cols) %>% 
  as.matrix() %*% as.matrix(ess_w)
plot(x_ess, y)
cor.test(x_ess,y)
l = lm(y ~ x_ess)
abline(l, col='red', lty=2)
