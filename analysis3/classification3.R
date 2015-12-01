# train SVM model on decoys:
# 1 if downstream ss gets used
# 0 if upstream ss gets used

library(e1071)
setwd('~/Desktop/MacArthur/de-novo-splice/analysis3')
source('auc.R'); source('extract_features.R')
v_calc = 'C'; vers='both'; source('load_features.R')
source('apply_weights.R')

# get training data + test decoys
data = get_feature_matrix(use, !use$intronic)
train = data$train
test_decoys = data$test
# model = svm(as.factor(y) ~ ., data=train[,c('mes_diff', 'ese_prox_w', 'ess_prox_w', 'y')], type='C', kernel='linear')
# test performance at classifying decoys
model1 = eval_svm(train[,c('mes_diff', 'ese_prox_w', 'ess_prox_w', 'y')], test_decoys, kernel='linear')
roc1 = get_roc_curve(model1$score, test_decoys$y) 
roc_mes = get_roc_curve(-test_decoys$mes_diff, test_decoys$y)
plot(roc1[,1], roc1[,2], col='blue', bty='n'); points(roc_mes[,1], roc_mes[,2], col='green')

# get all the features!!
train = get_density_feats(train)
test_decoys = get_density_feats(test_decoys)

all_the_feats = c('mes_diff', 'ese_prox_dens', 'ess_prox_dens', 'gtrip_prox_dens', 'ese_up_dens', 'ess_up_dens', 
                  'gtrip_up_dens', 'ese_down_dens', 'ess_down_dens', 'gtrip_down_dens', 'y')
model2 = eval_svm(train[,all_the_feats], test_decoys, kernel='linear')
roc2 = get_roc_curve(model2$score, test_decoys$y) 

# plot ROC curves
plot(roc1[,1], roc1[,2], col='blue', bty='n'); points(roc2[,1], roc2[,2], col='red'); points(roc_mes[,1], roc_mes[,2], col='green')
legend("bottomright", c("svm ess/ese", "svm all the features!", "mes_diff"), 
       fill=c("blue","red","green"), bty='n', cex=0.8)

# test on DBASS data
dbass = get_feature_matrix(novo_feats, novo_feats$intronic, ntestdecoys=0)$train
dbass = get_density_feats(dbass)
test = rbind(dbass, test_decoys) # combine DBASS with test decoys

# generate models and ROC performance
alt_in_exon = T
use_test = if (alt_in_exon) subset(test, !intronic) else subset(test, intronic)
model1 = eval_svm(train[,c('mes_diff', 'ese_prox_w', 'ess_prox_w', 'y')], use_test, kernel='linear')
roc1 = mut_roc(get_roc_curve(model1$score, use_test$y), alt_in_exon)
model2 = eval_svm(train[,all_the_feats], use_test, kernel='linear')
roc2 = mut_roc(get_roc_curve(model2$score, use_test$y), alt_in_exon)
roc_mes = get_roc_curve(use_test$mes_diff, !use_test$y)
plot(roc_mes[,1], roc_mes[,2])

# plot ROC curves 
plot(roc1[,1], roc1[,2], main='De Novo detection', xlab='fpr', ylab='tpr', col='red', bty='n'); points(roc2[,1], roc2[,2], col='blue'); points(roc_mes[,1], roc_mes[,2], col='purple')
legend("bottomright", c("svm ESE/ESS", "svm all the feats!", "mes diff"),
       fill=c("red", "blue", "purple"),bty='n',cex=0.85)


