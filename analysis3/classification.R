# learn model to classify given pair of splice sites as de novo or decoy

library(e1071) # for SVM classification
library(rgl) # for 3-d visualization of features
source('~/Desktop/MacArthur/de-novo-splice/analysis2/auc.R')
source('~/Desktop/MacArthur/de-novo-splice/analysis2/load_features.R')
source('~/Desktop/MacArthur/de-novo-splice/analysis2/apply_weights.R')

get_cryp_train_data <- function(ntestdecoy=3000, ntraindecoy=5000, viz=F) {
  # get features
  aberrant = as.data.frame(rbind(novo_feats, cryptic_feats))
  is_cryptic = c(rep(F, nrow(novo_feats)), rep(T, nrow(cryptic_feats)))
  decoy = use
  X = rbind(decoy, aberrant)
  is_cryptic = c(rep(F, nrow(decoy)), is_cryptic)
  feats = c("ProximalESEWeighted", "ProximalESSWeighted", "ProximalESE", "ProximalESS")
  X = cbind(X[,feats], X$AltMES - X$RefMES, X$AltMES, X$intronic)
  y = c(rep(F, nrow(decoy)), rep(T, nrow(aberrant)))
  colnames(X) = c('ese', 'ess', 'ese_uw', 'ess_uw', 'mes_diff', 'alt_mes', 'intronic')
  
  # visualize features
  if (viz) {
    idx = X$intronic
    x = subset(X, idx)
    colorer <- function(i) {
      key = y[i] + is_cryptic[i] + 1
      switch(key, '#AAAAAA', '#9D1309', '#FF6103')
    }
    colors=sapply(which(idx), colorer)
    plot(x$mes_diff, x$ese, col=colors)
    plot(x$mes_diff, x$ess, col=colors)
    # plot3d(x$mes_diff, x$ese, x$ess, col=colors) 
  }
  
  # break into training and testing data
  X = cbind(X, y)
  test_decoy_idx = sample(1:nrow(decoy), ntestdecoy)
  test_decoys = X[test_decoy_idx,]
  train_decoy_idx = !(1:nrow(decoy) %in% test_decoy_idx)
  train_decoys = X[train_decoy_idx,]
  cryptic = X[is_cryptic,]
  denovo = X[!is_cryptic & y,]
  # final data set
  training_set = rbind(cryptic, sample_n(train_decoys, ntraindecoy))
  testing_set = rbind(denovo, test_decoys)
  return(list(train=training_set, test=testing_set))
}

main <- function() {
  data = get_cryp_train_data()
  training_set = data$train
  testing_set = data$test
  
  # compare different models
  alt_in_exon = T
  if (alt_in_exon) {
    train = subset(training_set, !intronic)
    test = subset(testing_set, !intronic)
  } else {
    train = subset(training_set, intronic)
    test = subset(testing_set, intronic)
  }
  roc1 = get_roc_curve(test$alt_mes, test$y)
  roc2 = get_roc_curve(test$mes_diff, test$y)
  model_a = eval_svm(train[,c('mes_diff','ese','ess','y')], test, kernel='linear')
  roc_a = get_roc_curve(model_a$score, test$y)
  model_b = eval_svm(train[,c('mes_diff','ese_uw','ess_uw','y')], test, kernel='linear')
  roc_b = get_roc_curve(model_b$score, test$y)
  model_c = eval_svm(train[,c('mes_diff','ese','ess_uw','y')], test, kernel='linear')
  roc_c = get_roc_curve(model_c$score, test$y)
  model_d = eval_svm(train[,c('mes_diff','ese_uw','ess','y')], test, kernel='linear')
  roc_d = get_roc_curve(model_d$score, test$y)
  
  # plot ROC curves
  plot(roc_a[,1], roc_a[,2])
  
  # visualize decision boundary
  plot(model_a$svm, train, ess ~ mes_diff, dataSymbol=NA)
  plot(model_a$svm, train, ese ~ mes_diff, dataSymbol=NA)
}

if (getOption('run.main', default=TRUE)) {
  main()
}
