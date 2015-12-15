library(magrittr)
library(rateratio.test)

# test = c("2,6,2,52", "436,3636,2")
# ragged = make_ragged(test)
# lapply(ragged, function(v, x) {v - x}, list(6, 8))
# up = mapply(subtract, ragged, as.list(c(5,7)))
# down = mapply(subtract, ragged, as.list(c(1,4)))

make_ragged <- function(col) {
  l = strsplit(col, ",")
  l = lapply(l, as.numeric)
  return(l)
} 

# WE KNOW THERE ARE HUGE SPIKES IN ESE/ESS DENSITY AT THE SPLICE SITES (GET THESE PLOTS)
# HOW DO WE USE THAT INFORMATION TO DRAW GOOD BOUNDARIES?
# https://cran.r-project.org/web/packages/rateratio.test/vignettes/rateratio.test.pdf

# first optimize size of buffer between counting region and splice site
count = function(v, a, b) { sum(v > a & v < b) }

# upstream elements
# x = 0:-40
f <- function(elem, x, from_up, k=-50, hard_lb=exon_start+3, make_plot=T) {
  m = y1 = y2 = vector(mode='numeric', length=length(x)); i=1
  for (ub in x) {
    start = apply(cbind(upstream + ub + k, hard_lb), 1, max)
    dist = upstream + ub - start - 1
    dist[dist < 0] = 0
    n = mapply(count, from_up, start - upstream, ub)
    hits = c(sum(n[exonic], na.rm=T), sum(n[!exonic], na.rm=T))
    time = c(sum(dist[exonic]), sum(dist[!exonic]))
    t = rateratio.test(hits, time)
    m[i] = t$estimate[1]
    y1[i] = t$conf.int[1]
    y2[i] = t$conf.int[2]
    i = i + 1
  }
  if (make_plot) {
    title = sprintf('ratio of %s density, max_window_size = %s', elem, abs(k))
    yax = 'exonic / intronic'
    xax = 'nt from upstream splice site'
    plot(x, m, type='l', lty=2, main=title, xlab=xax, ylab=yax, ylim=c(min(y1),max(y2)))
    polygon(c(x,rev(x)),c(y2,rev(y1)),col=rgb(1,0,0,0.5))
  }
  return(m)
}

# downstream elements
# x = 0:40
h <- function(elem, x, from_down, k=50, hard_ub=intron_end-30, make_plot=T) {
  m = y1 = y2 = vector(mode='numeric', length=length(x)); i=1
  for (lb in x) {
    end = apply(cbind(downstream + lb + k, hard_ub), 1, min)
    dist = end - (downstream + lb) - 1
    dist[dist < 0] = 0
    n = mapply(count, from_down, lb, end - downstream)
    hits = c(sum(n[exonic], na.rm=T), sum(n[!exonic], na.rm=T))
    time = c(sum(dist[exonic]), sum(dist[!exonic]))
    t = rateratio.test(hits, time)
    m[i] = t$estimate[1]
    y1[i] = t$conf.int[1]
    y2[i] = t$conf.int[2]
    i = i + 1
  }
  if (make_plot) {
    title = sprintf('ratio of %s density, max_window_size = %s', elem, abs(k))
    yax = 'exonic / intronic'
    xax = 'nt from downstream splice site'
    plot(x, m, type='l', lty=2, main=title, xlab=xax, ylab=yax, ylim=c(min(y1),max(y2)))
    polygon(c(x,rev(x)),c(y2,rev(y1)),col=rgb(1,0,0,0.5))
  }
  return(m)
}

# x = 0:40
g1 <- function(elem, x, from_up, k=50, hard_ub=downstream-10, make_plot=T) {
  m = y1 = y2 = vector(mode='numeric', length=length(x)); i=1
  for (lb in x) {
    dist = hard_ub - (upstream + lb) - 1
    dist[dist < 0] = 0
    n = mapply(count, from_up, lb, hard_ub - upstream)
    hits = c(sum(n[exonic], na.rm=T), sum(n[!exonic], na.rm=T))
    time = c(sum(dist[exonic]), sum(dist[!exonic]))
    t = rateratio.test(hits, time)
    m[i] = t$estimate[1]
    y1[i] = t$conf.int[1]
    y2[i] = t$conf.int[2]
    i = i + 1
  }
  if (make_plot) {
    title = sprintf('ratio of %s density, proximal', elem, abs(k))
    yax = 'exonic / intronic'
    xax = 'nt from upstream splice site'
    plot(x, m, type='l', lty=2, main=title, xlab=xax, ylab=yax, ylim=c(min(y1),max(y2)))
    polygon(c(x,rev(x)),c(y2,rev(y1)),col=rgb(1,0,0,0.5))
  }
  return(m)
}

# x = 0:-40
g2 <- function(elem, x, from_down, k=50, hard_lb=upstream+3, make_plot=T) {
  m = y1 = y2 = vector(mode='numeric', length=length(x)); i=1
  for (ub in x) {
    dist = (downstream + ub) - hard_lb - 1
    dist[dist < 0] = 0
    n = mapply(count, from_down, hard_lb - downstream, ub)
    hits = c(sum(n[exonic], na.rm=T), sum(n[!exonic], na.rm=T))
    time = c(sum(dist[exonic]), sum(dist[!exonic]))
    t = rateratio.test(hits, time)
    m[i] = t$estimate[1]
    y1[i] = t$conf.int[1]
    y2[i] = t$conf.int[2]
    i = i + 1
  }
  if (make_plot) {
    title = sprintf('ratio of %s density, proximal', elem, abs(k))
    yax = 'exonic / intronic'
    xax = 'nt from downstream splice site'
    plot(x, m, type='l', lty=2, main=title, xlab=xax, ylab=yax, ylim=c(min(y1),max(y2)))
    polygon(c(x,rev(x)),c(y2,rev(y1)),col=rgb(1,0,0,0.5))
  }
  return(m)
}


