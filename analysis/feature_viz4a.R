source('~/Desktop/MacArthur/DBASS/analysis/feature_calc.R')

smooth1 <- function(df, x, sre, b=seq(0, 10, 0.1)) {
  bins = cut(x, breaks=b, labels=F)
  y = vector(mode='numeric', length=max(bins, na.rm=T))
  for (i in 1:max(bins, na.rm=T)) {
    bin = df[bins==i & !is.na(bins),]
    y[i] = sum(bin[,sre])/sum(abs(bin$AuthenticLocation - bin$DeNovoLocation))
  }
  x = 1:max(bins, na.rm=T)
  plot(x, y, ylab=paste(sre, " Density"), xlab="Authentic MES - Decoy MES bins")
  cor.test(x, y)
}

smooth2 <- function(df, y, xlab, b=seq(0, 0.2, 0.01)) {
  bins = cut(x, breaks=b, labels=F)
  y = vector(mode='numeric', length=max(bins, na.rm=T))
  for (i in 1:max(bins, na.rm=T)) {
    bin = df[bins==i & !is.na(bins),]
    y[i] = mean(bin$AuthenticMES - bin$DeNovoMES)
  }
  x = 1:max(bins, na.rm=T)
  plot(x, y, ylab="Authentic - Decoy MES", xlab=xlab)
  cor.test(x, y)
}

# exonic decoys
# inner ESE density vs. relative ss strength (i.e. authentic - decoy MES)
ese_dens = euse$InnerESE/(euse$AuthenticLocation - euse$DeNovoLocation)
mes_diff = euse$AuthenticMES - euse$DeNovoMES
plot(mes_diff, ese_dens) # f(mes_diff) = ese_dens
cor.test(mes_diff, ese_dens)
plot(ese_dens, mes_diff) # f(ese_dens) = mes_diff
# fit regerssion line
cor.test(ese_dens, mes_diff)
l = lm(mes_diff ~ ese_dens)
abline(l, col='red', lty=2)
stats = cor.test(ese_dens, mes_diff)
margin = stats$conf.int[2] - stats$conf.int[1]
legend("bottomright", cex=0.8, bty='n',
       sprintf("r = %.3f +/- %.3f", stats$estimate, margin))

# smooth by binning
smooth1(euse, mes_diff, "InnerESE")
smooth2(euse, ese_dens,"Inner ESE density bins")

# intronic decoys
# inner ESS density vs. relative ss strength (i.e. authentic - decoy MES)
ess_dens = iuse$InnerESS/(iuse$DeNovoLocation - iuse$AuthenticLocation)
mes_diff = iuse$AuthenticMES - iuse$DeNovoMES
plot(mes_diff, ess_dens) # f(mes_diff) = ess_dens
cor.test(mes_diff, ess_dens)
plot(ess_dens, mes_diff) # f(ess_dens) = mes_diff
# fit linear regression
l = lm(mes_diff ~ ess_dens)
abline(l, col='red', lty=2)
stats = cor.test(ess_dens, mes_diff)
margin = stats$conf.int[2] - stats$conf.int[1]
legend("bottomright", cex=0.8, bty='n',
       sprintf("r = %.3f +/- %.3f", stats$estimate, margin))

# smooth by binning
smooth1(iuse, mes_diff, "InnerESS", b=seq(-2, 10, 0.1))
smooth2(iuse, ess_dens, "Inner ESS density bins", b=seq(0, 0.08, 0.005))
