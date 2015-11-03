# Inner ESS distribution for intronic decoys

source('~/Desktop/MacArthur/DBASS/analysis/feature_cacl.R')

plot_norm <- function(h, binned, before=T, x='Dist from decoy') {
  if (before) {
    cdf = ecdf(binned)
    norm = cdf(1:(length(h$breaks)-1)) # for both
  } else {
    cdf = ecdf(-binned)
    norm = cdf(-1:-(length(h$breaks)-1))
  }
  norm = norm*length(binned)
  y = h$counts/norm
  barplot(y, names.arg=h$mids, xlab=x, ylab='Normalized frequency of ESSs',
          main='Normalized distribution of Inner ESSs for intronic decoys')
}

get_inner_ess_dist <- function(row, rel) {
  ess_loc = strsplit(row["ESSLoc"], ",")[[1]] %>% as.numeric
  inner = ess_loc > as.numeric(row["AuthenticLocation"]) & ess_loc < as.numeric(row["DeNovoLocation"])
  ess_loc = ess_loc[inner]
  return(ess_loc - as.numeric(row[rel]))
}

get_inner_ess_perc <- function(row, start=0, end=0) {
  ess_loc = strsplit(row["ESSLoc"], ",")[[1]] %>% as.numeric
  inner = ess_loc > (as.numeric(row["AuthenticLocation"]) + start) & ess_loc < (as.numeric(row["DeNovoLocation"]) - end)
  ess_loc = ess_loc[inner]
  p = (ess_loc - (as.numeric(row["AuthenticLocation"]) + start)) / (as.numeric(row["DeNovoDist"]) - (start+end))
  return(p)
}

# before decoy
ess_dist = apply(iuse, 1, get_inner_ess_dist, "DeNovoLocation") %>% do.call(c, .)
high_res_bins=seq(0.5, 50.5, by=1)
h = hist(ess_dist, plot=T, breaks=-high_res_bins, main='Distribution of Inner ESS for intronic decoys')
binned = cut(-iuse$DeNovoDist, breaks=h$breaks, labels=F) 
plot_norm(h, binned, before=T)

# after authentic
ess_dist = apply(iuse, 1, get_inner_ess_dist, "AuthenticLocation") %>% do.call(c, .)
high_res_bins=seq(0.5, 50.5, by=1)
h = hist(ess_dist, plot=T, breaks=high_res_bins, main='Distribution of Inner ESS for intronic decoys')
binned = cut(iuse$DeNovoDist, breaks=h$breaks, labels=F) 
plot_norm(h, binned, before=F, x='Dist from Authentic')

# percent dist between authentic and decoy
ess_perc1 = apply(iuse, 1, get_inner_ess_perc, start=6, end=9) %>% do.call(c, .)
hist(ess_perc1, breaks=100,
     xlab='Percentile distance between authentic and decoy splice site')
