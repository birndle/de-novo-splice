error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper)) {
    stop("vectors must be same length")
  }
  arrows(x,y+upper,x,y-lower, angle=90, code=3, length=length)
}

plot_bars <- function(novo, decoy, title, xax, colors=c('#9D1309','#AAAAAA')) {
  bars = rbind(novo[,1], decoy[,1])
  margins = rbind(novo[,2], decoy[,2])
  ymax = max(bars) + max(margins)
  bp = barplot(bars, beside=T, main=title, col=colors, 
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

get_plot_data2 <- function(decoy, sections, sre, pois=T) {
  no = matrix(nrow=0, ncol=2)
  for (loc in sections) {
    kcol = paste(loc, sre, sep="")
    distcol = paste(loc, "Dist", sep="")
    if (pois) {
      k = sum(decoy[,kcol])
      n = sum(decoy[,distcol])
      no = rbind(no, pois_ci(k, n))
    } else {
      k = decoy[,kcol]
      no = rbind(no, mean_ci(k))
    }
  }
  return(no)
}

mean_ci <- function(k) {
  mu = mean(k)
  n = length(k)
  s = sqrt(var(k))
  se = s/sqrt(n)
  margin = 1.96*se
  return(c(mu, margin))
}

ci <- function(k, n) {
  p = k/n
  margin = 1.96*p*(1-p)/sqrt(n)
  return(c(p, margin))
}

pois_ci <- function(k, n) {
  p = k/n
  if (p < 0) {
    se = 0
  } else {
    se = sqrt(p/n) 
  }
  margin = 1.96*se
  return(c(p, margin))
}
