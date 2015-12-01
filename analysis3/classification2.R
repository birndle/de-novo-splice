# train SVM model on decoys:
# 1 if downstream ss gets used
# 0 if upstream ss gets used

library(e1071) # for SVM classification
library(rgl) # for 3-d visualization of features

source('~/Desktop/MacArthur/de-novo-splice/analysis2/auc.R')
source('~/Desktop/MacArthur/de-novo-splice/analysis2/extract_features.R')
source('~/Desktop/MacArthur/de-novo-splice/analysis2/load_features.R')
source('~/Desktop/MacArthur/de-novo-splice/analysis2/apply_weights.R')

main <- function(alt_in_exon, viz=F) {
  # get training data + test decoys
  data = get_feature_matrix(use, !use$intronic)
  train = data$train
  test_decoys = data$test
  # train model
  model = svm(as.factor(y) ~ ., data=train[,c('mes_diff', 'ese_prox_w', 'ess_prox_w', 'y')], type='C', kernel='linear')
  
  # visualize decision boundary
  if (viz) {
    plot(model, data=train, ese_prox_w ~ mes_diff, dataSymbol=NA)
    plot(model, data=train, ess_prox_w ~ mes_diff, dataSymbol=NA) 
  }
  
  # add DBASS data to testing set of decoys
  test = get_feature_matrix(novo_feats, novo_feats$intronic, ntestdecoys=0)$train
  test = rbind(test, test_decoys) # combine DBASS with test decoys
  
  # before using that, test on just decoys
  model_w = eval_svm(train[,c('mes_diff', 'ese_prox_w', 'ess_prox_w', 'y')], test_decoys, kernel='linear')
  roc_w = get_roc_curve(model_w$score, test_decoys$y) 
  model_uw = eval_svm(train[,c('mes_diff', 'ese_prox', 'ess_prox', 'y')], test_decoys, kernel='linear')
  roc_uw = get_roc_curve(model_uw$score, test_decoys$y) 
  roc_mes = get_roc_curve(-test_decoys$mes_diff, test_decoys$y)
  plot(roc_w[,1], roc_w[,2], col='blue'); points(roc_uw[,1], roc_uw[,2], col='red'); points(roc_mes[,1], roc_mes[,2], col='green')
  legend("bottomright", c("svm weighted", "svm unweighted", "mes_diff"), fill=c("blue","red","green"), bty='n', cex=0.8)
  
  # visualize test decoys
  if (viz) {
    n = length(test_decoys$y)
    col = rep('red', n)
    col[test_decoys$y] = 'blue'
    plot(model, test_decoys, ese_up ~ mes_diff)
    plot(model, test_decoys, ess_prox_w ~ mes_diff) 
  }
  
  # now test on DBASS data
  cryp = get_cryp_train_data()
  cryp_train = cryp$train
  if (alt_in_exon) {
  	use_test = subset(test, !intronic)
  	use_test_cryp = use_test
  	use_test_cryp$y = !use_test_cryp$y
  	use_test_cryp$mes_diff = -use_test$mes_diff
  	use_cryp_train = subset(cryp_train, !intronic)
  } else {
  	use_test = subset(test, intronic)
  	use_test_cryp = use_test
  	use_cryp_train = subset(cryp_train, intronic)
  }
  
  # visualize testing data against decision boundary
  if (viz) {
    plot(model, data=use_test, ese ~ mes_diff)
    plot(model, data=use_test, ess ~ mes_diff) 
  }
  
  # generate models and ROC performance
  model_w = eval_svm(train[,c('mes_diff', 'ese_prox_w', 'ess_prox_w', 'y')], use_test, kernel='linear')
  roc_w = mut_roc(get_roc_curve(1 - model_w$score, use_test$y), alt_in_exon)
  model_uw = eval_svm(train[,c('mes_diff', 'ese_prox', 'ess_prox', 'y')], use_test, kernel='linear')
  roc_uw = mut_roc(get_roc_curve(1 - model_uw$score, use_test$y), alt_in_exon)
  roc_mes = mut_roc(get_roc_curve(-use_test$mes_diff, use_test$y), alt_in_exon)
  model_cryp = eval_svm(use_cryp_train[,c('mes_diff', 'ese', 'ess', 'y')], use_test_cryp, kernel='linear')
  roc_cryp = get_roc_curve(model_cryp$score, use_test_cryp$y)
  
  # plot ROC curves 
  plot(roc_w[,1], roc_w[,2], main='De Novo detection', xlab='fpr', ylab='tpr', col='red', bty='n')
  points(roc_uw[,1], roc_uw[,2], col='blue')
  points(roc_mes[,1], roc_mes[,2], col='purple')
  points(roc_cryp[,1], roc_cryp[,2], col='green')
  legend("bottomright", c("decoy unweighted", "decoy weighted", "mes diff", "cryptic weighted"),
         fill=c("blue", "red", "purple", "green"),bty='n',cex=0.85)
}

if (getOption('run.main', default=TRUE)) {
  main()
}
