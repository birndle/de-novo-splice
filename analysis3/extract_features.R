feat_cols = c("ProximalESEWeighted", "ProximalESSWeighted", "ProximalESE", "ProximalESS", "ProximalGtrip",
              "DistalESEWeighted", "DistalESSWeighted", "DistalESE", "DistalESS", "DistalGtrip",
              "ExonicESEWeighted", "ExonicESSWeighted", "ExonicESE", "ExonicESS", "ExonicGtrip",
              "IntronicESEWeighted", "IntronicESSWeighted", "IntronicESE", "IntronicESS", "IntronicGtrip",
              "ProximalDist", "DistalDist", "ExonicDist", "IntronicDist")
new_feat_names = c("ese_prox_w", "ess_prox_w", "ese_prox", "ess_prox", "gtrip_prox",
               "ese_dist_w", "ess_dist_w", "ese_dist", "ess_dist", "gtrip_dist",
               "ese_exon_w", "ess_exon_w", "ese_exon", "ess_exon", "gtrip_exon",
               "ese_intron_w", "ess_intron_w", "ese_intron", "ess_intron", "gtrip_intron",
               "prox_d", "dist_d", "exon_d", "intron_d")
               
get_feature_matrix <- function(m, y, ntestdecoys=3000, viz=F, deterministic=F) {
  X = m[,feat_cols]
  # X = as.data.frame(apply(X, 2, function(col) { col[col < 0] = 0; return(col)}))
  if (FALSE) {
    e = 0.001
    X = (X + e)/(use[,"ProximalDist"] + e)
  }
  mes_diff = (m$AltMES - m$RefMES)*m$intronic + (m$RefMES - m$AltMES)*(!m$intronic)
  X = cbind(X, mes_diff, y, m$intronic)
  colnames(X) = c(new_feat_names, 'mes_diff', 'y', 'intronic')
  
  # visualize features
  if (viz) {
    n = length(y)
    col = rep('red', n)
    col[y] = 'blue'
    # blue if downstream gets used, red if upstream
    plot(X$mes_diff, X$ese_prox_w, col=col) 
    plot(X$mes_diff, X$ess_prox_w, col=col) 
    # plot3d(X$mes_diff, X$ese, X$ess, col=col) # blue if upstream gets used, red if downstream
  }
  
  # divide decoys into testing and training sets
  if (deterministic) {
    set.seed(1)
  }
  test_decoy_idx = sample(1:nrow(X), ntestdecoys)
  test_decoys = X[test_decoy_idx,]
  train_decoy_idx = (1:nrow(X))[!(1:nrow(X) %in% test_decoy_idx)]
  train = X[train_decoy_idx,]
  return(list(train=as.data.frame(train), test=as.data.frame(test_decoys)))
}

get_density_feats <- function(tmp) {
  # compute densities for SRE elements in order to make valid comparisons
  # unweighted
  prox = tmp[,c("ese_prox", "ess_prox", "gtrip_prox")]
  prox_d = tmp[,"prox_d"]
  prox_d[prox_d < 1] = 1
  tmp[,c("ese_prox_dens", "ess_prox_dens", "gtrip_prox_dens")] = prox/prox_d
  up = tmp$intronic*tmp[,c("ese_exon", "ess_exon", "gtrip_exon")] + (!tmp$intronic)*tmp[,c("ese_dist", "ess_dist", "gtrip_dist")]
  up_d = tmp$intronic*tmp$exon_d + (!tmp$intronic)*tmp$dist_d
  up_d[up_d < 1] = 1
  tmp[,c("ese_up", "ess_up", "gtrip_up")] = up # absolute 
  tmp[,c("ese_up_dens", "ess_up_dens", "gtrip_up_dens")] = up/up_d
  down = tmp$intronic*tmp[,c("ese_dist", "ess_dist", "gtrip_dist")] + (!tmp$intronic)*tmp[,c("ese_intron", "ess_intron", "gtrip_intron")]
  down_d = tmp$intronic*tmp$dist_d + (!tmp$intronic)*tmp$intron_d
  down_d[down_d < 1] = 1
  tmp[,c("ese_down", "ess_down", "gtrip_down")] = down # absolute
  tmp[,c("ese_down_dens", "ess_down_dens", "gtrip_down_dens")] = down/down_d
  # weighted
  prox = tmp[,c("ese_prox_w", "ess_prox_w")]
  prox_d = tmp[,"prox_d"]
  prox_d[prox_d < 1] = 1
  tmp[,c("ese_prox_dens_w", "ess_prox_dens_w")] = prox/prox_d
  up = tmp$intronic*tmp[,c("ese_exon_w", "ess_exon_w")] + (!tmp$intronic)*tmp[,c("ese_dist", "ess_dist")]
  up_d = tmp$intronic*tmp$exon_d + (!tmp$intronic)*tmp$dist_d
  up_d[up_d < 1] = 1
  tmp[,c("ese_up_w", "ess_up_w")] = up # absolute
  tmp[,c("ese_up_dens_w", "ess_up_dens_w")] = up/up_d
  down = tmp$intronic*tmp[,c("ese_dist", "ess_dist")] + (!tmp$intronic)*tmp[,c("ese_intron_w", "ess_intron_w")]
  down_d = tmp$intronic*tmp$dist_d + (!tmp$intronic)*tmp$intron_d
  down_d[down_d < 1] = 1
  tmp[,c("ese_down_w", "ess_down_w")] = down # absolute
  tmp[,c("ese_down_dens_w", "ess_down_dens_w")] = down/down_d
  return(tmp)
}