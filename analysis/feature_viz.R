# SRE density comparison for decoy vs. de novo

source('~/Desktop/MacArthur/DBASS/analysis/feature_calc.R')

# compare De Novo to Authentic MES scores
hist(feats$DeNovoMES-feats$AuthenticMES, breaks=20, prob=T,
     main='De Novo ss strength vs. Authentic ss strength',
     xlab='De Novo MES - Authentic MES',
     col=rgb(1,0,0,0.5), ylim=c(0,0.2))
hist(use$DeNovoMES-use$AuthenticMES, breaks=10, prob=T, add=T, col=rgb(0,0,1,0.5))
legend('topleft', bty='n', c('De Novo', 'Decoy'), fill=c('red', 'blue'))

hist(feats$AuthenticMES, breaks=20, prob=T, main='Authentic ss strength', 
     xlab='Authentic MES',col=rgb(1,0,0,0.5), ylim=c(0,0.3))
hist(use$AuthenticMES, breaks=50, prob=T, add=T, col=rgb(0,0,1,0.5))
legend('topleft', bty='n', c('De Novo', 'Decoy'), fill=c('red', 'blue'))

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper)) {
    stop("vectors must be same length")
  }
  arrows(x,y+upper,x,y-lower, angle=90, code=3, length=length)
}

plot_bars <- function(novo, decoy, title, xax) {
  bars = rbind(novo[,1], decoy[,1])
  margins = rbind(novo[,2], decoy[,2])
  ymax = max(bars) + max(margins)
  bp = barplot(bars, beside=T, main=title, col=c('#9D1309','#AAAAAA'), 
               ylim=c(0, ymax), space=c(0, 2),
               names.arg=xax)
  error.bar(bp, bars, margins)
}

get_plot_data <- function(novo, decoy, sections, sre) {
  no = matrix(nrow=0, ncol=2)
  de = matrix(nrow=0, ncol=2)
  for (loc in sections) {
    col1 = paste(loc, sre, sep="")
    col2 = paste(loc, "Dist", sep="")
    k = sum(novo[col1])
    n = sum(novo[col2])
    no = rbind(no, pois_ci(k, n))
    k = sum(decoy[col1])
    n = sum(decoy[col2])
    de = rbind(de, pois_ci(k,n))
  }
  return(cbind(no, de))
}

# ESEs for intronic alternative
subset(feats, intronic) %>% 
  mutate(ExonicDist=AuthenticLocation-ExonStart, InnerDist=DeNovoLocation-AuthenticLocation, OuterDist=IntronEnd-DeNovoLocation) %>% 
  select(ExonicESE, InnerESE, OuterESE, ExonicDist, InnerDist, OuterDist) %>%
  colSums -> novo
subset(use, intronic) %>% 
  mutate(ExonicDist=AuthenticLocation-ExonStart, InnerDist=DeNovoLocation-AuthenticLocation, OuterDist=IntronEnd-DeNovoLocation) %>% 
  select(ExonicESE, InnerESE, OuterESE, ExonicDist, InnerDist, OuterDist) %>%
  colSums -> decoy
sections = c("Exonic", "Inner", "Outer")
m = get_plot_data(novo, decoy, sections, "ESE")
plot_bars(m[,1:2], m[,3:4], 'ESE for intronic alternative', xax=sections)

# ESE for exonic alternative
subset(feats, !intronic) %>% 
  mutate(OuterDist=DeNovoLocation-ExonStart, InnerDist=AuthenticLocation-DeNovoLocation, IntronicDist=IntronEnd-AuthenticLocation) %>% 
  select(OuterESE, InnerESE, IntronicESE, OuterDist, InnerDist, IntronicDist) %>%
  colSums -> novo
subset(use, !intronic) %>% 
  mutate(OuterDist=DeNovoLocation-ExonStart, InnerDist=AuthenticLocation-DeNovoLocation, IntronicDist=IntronEnd-AuthenticLocation) %>% 
  select(OuterESE, InnerESE, IntronicESE, OuterDist, InnerDist, IntronicDist) %>%
  colSums -> decoy
sections = c("Outer", "Inner", "Intronic")
m = get_plot_data(novo, decoy, sections, "ESE")
plot_bars(m[,1:2], m[,3:4], 'ESE for exonic alternative', xax=sections)

# ESS for intronic alternative
subset(feats, intronic) %>%
  mutate(OuterDist=IntronEnd-DeNovoLocation, InnerDist=DeNovoLocation-AuthenticLocation, ExonicDist=AuthenticLocation-ExonStart) %>%
  select(OuterESS, InnerESS, ExonicESS, OuterDist, InnerDist, ExonicDist) %>%
  colSums -> novo
subset(use, intronic) %>%
  mutate(OuterDist=IntronEnd-DeNovoLocation, InnerDist=DeNovoLocation-AuthenticLocation, ExonicDist=AuthenticLocation-ExonStart) %>%
  select(OuterESS, InnerESS, ExonicESS, OuterDist, InnerDist, ExonicDist) %>%
  colSums -> decoy
sections = c("Exonic", "Inner", "Outer")
m = get_plot_data(novo, decoy, sections, "ESS")
plot_bars(m[,1:2], m[,3:4], 'ESS for intronic alternative', xax=sections)

# ESS for exonic alternative
subset(feats, !intronic) %>% 
  mutate(OuterDist=DeNovoLocation-ExonStart, InnerDist=AuthenticLocation-DeNovoLocation, IntronicDist=IntronEnd-AuthenticLocation) %>% 
  select(OuterESS, InnerESS, IntronicESS, OuterDist, InnerDist, IntronicDist) %>%
  colSums -> novo
subset(use, !intronic) %>% 
  mutate(OuterDist=DeNovoLocation-ExonStart, InnerDist=AuthenticLocation-DeNovoLocation, IntronicDist=IntronEnd-AuthenticLocation) %>% 
  select(OuterESS, InnerESS, IntronicESS, OuterDist, InnerDist, IntronicDist) %>%
  colSums -> decoy
sections = c("Outer", "Inner", "Intronic")
m = get_plot_data(novo, decoy, sections, "ESS")
plot_bars(m[,1:2], m[,3:4], 'ESS for exonic alternative', xax=sections)

# G-triplet for intronic alternative
subset(feats, intronic) %>%
  mutate(OuterDist=IntronEnd-DeNovoLocation, InnerDist=DeNovoLocation-AuthenticLocation, ExonicDist=AuthenticLocation-ExonStart) %>%
  select(OuterGtrip, InnerGtrip, ExonicGtrip, OuterDist, InnerDist, ExonicDist) %>%
  colSums -> novo
subset(use, intronic) %>%
  mutate(OuterDist=IntronEnd-DeNovoLocation, InnerDist=DeNovoLocation-AuthenticLocation, ExonicDist=AuthenticLocation-ExonStart) %>%
  select(OuterGtrip, InnerGtrip, ExonicGtrip, OuterDist, InnerDist, ExonicDist) %>%
  colSums -> decoy
sections = c("Exonic", "Inner", "Outer")
m = get_plot_data(novo, decoy, sections, "Gtrip")
plot_bars(m[,1:2], m[,3:4], 'G triplet for intronic alternative', xax=sections)

# G-triplet of exonic alternative
subset(feats, !intronic) %>% 
  mutate(OuterDist=DeNovoLocation-ExonStart, InnerDist=AuthenticLocation-DeNovoLocation, IntronicDist=IntronEnd-AuthenticLocation) %>% 
  select(OuterGtrip, InnerGtrip, IntronicGtrip, OuterDist, InnerDist, IntronicDist) %>%
  colSums -> novo
subset(use, !intronic) %>% 
  mutate(OuterDist=DeNovoLocation-ExonStart, InnerDist=AuthenticLocation-DeNovoLocation, IntronicDist=IntronEnd-AuthenticLocation) %>% 
  select(OuterGtrip, InnerGtrip, IntronicGtrip, OuterDist, InnerDist, IntronicDist) %>%
  colSums -> decoy
sections = c("Outer", "Inner", "Intronic")
m = get_plot_data(novo, decoy, sections, "Gtrip")
plot_bars(m[,1:2], m[,3:4], 'G triplet for exonic alternative', xax=sections)

