get_dist <- function(row, rel, sre="ESSLoc") {
  ess_loc = strsplit(row[sre], ",")[[1]] %>% as.numeric
  return(ess_loc - as.numeric(row[rel]))
}

ess_dist = subset(use, intronic & AltDist > 50 & AltDist < 100) %>% apply(., 1, get_dist, "RefLoc") %>% do.call(c, .)
hist(ess_dist[(ess_dist < 50 & ess_dist > 6) | (ess_dist > -50 & ess_dist < -4)], 
     breaks=seq(-50.5, 50.5, 1), prob=T, col=rgb(0,0,1,0.5))

ese_dist = subset(use, !intronic) %>% apply(., 1, get_dist, "RefLoc", sre="ESELoc") %>% do.call(c, .)
hist(ese_dist[(ese_dist < 50 & ese_dist > 6) | (ese_dist > -50 & ese_dist < -6)], 
     breaks=seq(-50.5, 50.5, 1), col=rgb(1,0,0,0.5), prob=T)
