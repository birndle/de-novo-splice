library(e1071)
# mes_diff is ALT - REF
get_roc_curve <- function(score, truth, res=1000) {
  thresh = seq(range(score)[1], range(score)[2], length.out=res)
  i = 1
  roc = matrix(nrow=length(thresh),ncol=2)
  for (t in thresh) {
    yhat2 = score < t
    fpr = 1 - mean(yhat2[!truth] == FALSE)
    tpr = mean(yhat2[truth] == T)
    roc[i,] = c(fpr, tpr)
    i = i + 1 
  }
  return(roc)
}

mut_roc <- function(roc, exonic) {
  if (exonic) {
    fpr = 1 - roc[,2]
    tpr = 1 - roc[,1]
    return(cbind(fpr, tpr))
  } else {
    return(roc)
  }
}

get_auc <- function(feats, trainset, testset, alt_in_exon=T) {
  train_use = trainset[,c(feats, 'y')]
  test_use = if (alt_in_exon) subset(testset, !intronic, select=c(feats, 'y')) else subset(testset, intronic, select=c(feats, 'y'))
  model = eval_svm(train_use, test_use, kernel='linear', verbose=F)
  roc = if (alt_in_exon) get_roc_curve(1 - model$score, !test_use$y, res=100) else get_roc_curve(model$score, test_use$y, res=100)
  area = auc(roc[,1], roc[,2])
  return(area)
}

auc <- function(x, y) {
  if (length(x) != length(y)) {
    message("ERROR: x and y have unequal lengths.")
    return(1)
  }
  old_x = 0
  old_y = 0
  area = 0
  for (i in 1:length(x)) {
    dx = x[i] - old_x
    dy = y[i] - old_y
    area = area + dx*old_y + 0.5*dx*dy    
    old_x = x[i]
    old_y = y[i]
  }
  return(area)
}

eval_svm <- function(train, test, kernel='radial', verbose=T) {
  model = svm(as.factor(y) ~ ., data=train, type='C', kernel=kernel, probability=T)
  yhat = predict(model, test, probability = T)
  probs = attr(yhat, "probabilities")
  fpr = 1 - mean(yhat[!test$y] == FALSE) # fpr
  tpr = mean(yhat[test$y] == TRUE) # tpr
  if (verbose) {
    message(sprintf("TPR = %s", tpr))
    message(sprintf("FPR = %s", fpr))  
  }
  return(list(svm=model, score=probs[,1]))
}
