# implementation of SVM-RFE
library(e1071)

################################################
# Feature Ranking with SVM-RFE/R-SVM
################################################
svmrfeFeatureRanking <- function(x,y, trainingFold=1, ensembleClassifiers=1, rsvm=FALSE){
  n <- ncol(x)
  survivingFeaturesIndexes <- seq(1:n)
  featureRankedList <- vector(length=n)
  rankedFeatureIndex <- n
  if(ensembleClassifiers > 1) {
    if(trainingFold > 1) {
      warning("Both ensembleClassifiers and trainingFold larger than 1. Set trainingFold=1")
      trainingFold <- 1
    }
  }
  
  while(length(survivingFeaturesIndexes)>1){
    if(trainingFold == 1) {
      if(ensembleClassifiers > 1) {
        y_nums <- sort(table(y))
        y_sel_idx <- foreach(i = 1:ensembleClassifiers) %do% {
          sort(c(which(y==names(y_nums)[1]), sample(which(y==names(y_nums)[2]), y_nums[1])))
        }
        results <- lapply(y_sel_idx, function(sel_idx) {
          svmModel <- svm(x[sel_idx, survivingFeaturesIndexes], y[sel_idx], cost = 10, cachesize=500,  
                          scale=F, type="C-classification", kernel="linear" )
          if(rsvm) {
            md <- colMeans(x[y[sel_idx]==levels(y)[1], survivingFeaturesIndexes, drop=FALSE]) - 
              colMeans(x[y[sel_idx]==levels(y)[2], survivingFeaturesIndexes, drop=FALSE])
            rankingCriteria <- t(svmModel$coefs) %*% svmModel$SV * md
          } else {
            #compute the weight vector
            w <- t(svmModel$coefs)%*%svmModel$SV
            #compute ranking criteria
            rankingCriteria <- w * w
          }
          ranking <- sort(rankingCriteria, index.return = TRUE)$ix
        })
        ranking  <- sort(apply(sapply(results, function(x) sort(x, index.return=T)$ix), 1, mean), index=T)$ix
      } else {
        #train the support vector machine
        svmModel <- svm(x[, survivingFeaturesIndexes], y, cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" )
        if(rsvm) {
          #yy <- vector(length = length(y)) 
          #yy[y==levels(y)[1]] <- 1
          #yy[y==levels(y)[2]] <- -1
          md <- colMeans(x[y==levels(y)[1], survivingFeaturesIndexes, drop=FALSE]) - colMeans(x[y==levels(y)[2], survivingFeaturesIndexes, drop=FALSE])
          #rankingCriteria <- t(svmModel$coefs*yy[svmModel$index]) %*% svmModel$SV * md
          rankingCriteria <- t(svmModel$coefs) %*% svmModel$SV * md
        } else {
          #compute the weight vector
          w <- t(svmModel$coefs)%*%svmModel$SV
          #compute ranking criteria
          rankingCriteria <- w * w
        }
        #rank the features
        ranking <- sort(rankingCriteria, index.return = TRUE)$ix
      }
    } else {
      nSample <- length(y)
      stopifnot(trainingFold <= nSample)
      tr_folds <- rep(1:trainingFold, len=nSample)[sample(nSample)]
      tr_folds <- lapply(1:trainingFold, function(x) which(tr_folds == x))
      results <- lapply(tr_folds, function(ex_id) {
        svmModel <- svm(x[-ex_id, survivingFeaturesIndexes], y[-ex_id], cost = 10, cachesize=500,  
                       scale=F, type="C-classification", kernel="linear" )
        if(rsvm) {
          md <- colMeans(x[y[-ex_id]==levels(y)[1], survivingFeaturesIndexes, drop=FALSE]) - 
            colMeans(x[y[-ex_id]==levels(y)[2], survivingFeaturesIndexes, drop=FALSE])
          rankingCriteria <- t(svmModel$coefs) %*% svmModel$SV * md
        } else {
          #compute the weight vector
          w <- t(svmModel$coefs)%*%svmModel$SV
          #compute ranking criteria
          rankingCriteria <- w * w
        }
        ranking <- sort(rankingCriteria, index.return = TRUE)$ix
      })
##      avg.rank <- sort(apply(sapply(results, function(x) sort(x, index.return=T)$ix), 1, mean), index=T)$x
      ranking  <- sort(apply(sapply(results, function(x) sort(x, index.return=T)$ix), 1, mean), index=T)$ix
    }
    
    #update feature ranked list
    featureRankedList[rankedFeatureIndex] <- survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex <- rankedFeatureIndex - 1
    
    #eliminate the feature with smallest ranking criterion
    (survivingFeaturesIndexes <- survivingFeaturesIndexes[-ranking[1]])
  }
  featureRankedList[rankedFeatureIndex] <- survivingFeaturesIndexes
  return (featureRankedList)
}

################################################
# Cross Validation with SVM-RFE/R-SVM
################################################
svmrfePrediction <- function(test_idx, x, y, steps, trainingFold=1, ensembleClassifiers=1, rsvm=FALSE) 
# test_idx: those sample will be excluded for training, only for testing 
# x: data
# y: labels
# steps: using those numbers of top features for prediction 
{ 
  featRankList <- svmrfeFeatureRanking(x[-test_idx,], y[-test_idx], trainingFold=trainingFold, 
                                       ensembleClassifiers=ensembleClassifiers, rsvm=rsvm)
  results <- matrix(0, nrow=4, ncol=length(steps))
  rownames(results) <- c("TP","TN", "FP", "FN")
  colnames(results) <- paste("Feat", steps, sep="_")
  for(i in 1:length(steps)){
    nfeatures <- steps[i]
    if(ensembleClassifiers > 1) {
      x_train <- x[-test_idx, featRankList[1:nfeatures], drop=FALSE]
      y_train <- y[-test_idx]
      y_nums <- sort(table(y_train))
      y_sel_idx <- foreach(k = 1:ensembleClassifiers) %do% {
        sort(c(which(y_train==names(y_nums)[1]), sample(which(y_train==names(y_nums)[2]), y_nums[1])))
      }
      pred <- sapply(y_sel_idx, function(sel_idx) {
        svmModel <- svm(x_train[sel_idx, ], y_train[sel_idx], cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" ) 
        predict(svmModel,x[test_idx, featRankList[1:nfeatures], drop=FALSE])
      })
      if(length(test_idx) == 1) {
        prediction <- factor(names(sort(table(pred), decreasing=T)[1]), levels=levels(y))
      } else {
        prediction <- factor(apply(pred, 1, function(x) names(sort(table(x), decreasing=T)[1]) ), levels=levels(y))
      }
    } else {
      svmModel <- svm(x[-test_idx, featRankList[1:nfeatures]], y[-test_idx], cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" ) 
      prediction <- predict(svmModel,x[test_idx, featRankList[1:nfeatures], drop=FALSE])
    }
    for(j in 1:length(prediction)) {
      if(prediction[j] == y[test_idx][j]) {
        if(y[test_idx][j] == levels(y)[1]) results["TP", i] = results["TP", i] + 1
        else results["TN", i] = results["TN", i] + 1
      } else {
        if(y[test_idx][j] == levels(y)[1]) results["FN", i] = results["FN", i] + 1
        else results["FP", i] = results["FP", i] + 1
      }
    }
  }
  list(featRankList=featRankList, predResults=results)
}

####################################
# Performance evaluation 
####################################
svmrfePerformance <- function(x, y, fold=10, trainingFold=1, ensembleClassifiers=1, rsvm=FALSE) {
  nSample <- length(y)
  nFeat <- ncol(x)
  if(fold==0) fold = nSample
  folds <- rep(1:fold, len=nSample)[sample(nSample)]
  folds <- lapply(1:fold, function(x) which(folds == x))
  
  steps <- unique(as.integer(exp(seq(0, log(nFeat), 0.1))))
  
  #results <- lapply(folds, svmrfePrediction, x, y, steps, trainingFold=trainingFold, rsvm=rsvm)
  results <- foreach(i = 1:fold) %dopar% {
    svmrfePrediction(folds[[i]], x, y, steps, trainingFold=trainingFold, ensembleClassifiers=ensembleClassifiers, rsvm=rsvm)
  }
  featRank <- do.call(rbind, lapply(results, function(x) x$featRankList))
  
  performanceRaw <- matrix(0, nrow=4, ncol=length(steps))
  rownames(performanceRaw) <- c("TP","TN", "FP", "FN")
  colnames(performanceRaw) <- paste("Feat", steps, sep="_")
  for(i in 1:4) performanceRaw[i,] <- t(rowSums(sapply(results, function(x) x$predResults[i,])))
  Acc <- colSums(performanceRaw[1:2,]) / colSums(performanceRaw)
  Sen <- performanceRaw["TP",] / colSums(performanceRaw[c("TP","FN"),])
  Spe <- performanceRaw["TN",] / colSums(performanceRaw[c("TN","FP"),])
  list(featRank=featRank, metrics=data.frame(rbind(Acc,Sen,Spe)))
}


##################
# Testing 
##################
# initialization for parallel computing 
library(doMC)
nCores <- 3  # specify this many cores to be used in this computing 
registerDoMC(cores = nCores)

# require(caret)
# data(mdrr)
# x <- as.matrix(mdrrDescr[,-nearZeroVar(mdrrDescr)])
# y <- as.factor(mdrrClass)

key <- "log.U24"
inputfile <- paste("sample211.expr", key, "txt", sep=".")
data <- read.table(file=inputfile, header=TRUE)
x <- as.matrix(data[,-1])
y <- as.factor(data[,1])

#scale features with mean zero and standard deviation one
x <- scale(x)
#x <- 2*atan(x/2)

# Feature Ranking with SVM-RFE (fold==0 means LOO)
rsvm.res <- svmrfePerformance(x, y, fold=6, trainingFold=2, ensembleClassifiers=11, rsvm=TRUE)
svmrfe.res <- svmrfePerformance(x, y, fold=6, trainingFold=2, ensembleClassifiers=11, rsvm=FALSE)

write.table(rsvm.res$metrics, paste("sample211", key, "classify.rsvm.txt",sep="."), sep="\t",quote=FALSE)
write.table(rsvm.res$featRank, paste("sample211", key, "classify.rsvm.feat.txt",sep="."), sep="\t",quote=FALSE)
write.table(svmrfe.res$metrics, paste("sample211", key, "classify.svmrfe.txt",sep="."), sep="\t",quote=FALSE)
write.table(svmrfe.res$featRank, paste("sample211", key, "classify.svmrfe.feat.txt",sep="."), sep="\t",quote=FALSE)

#sort(apply(apply( svmrfe.res$featRank, 1,  function(x) sort(x, index.return=T)$ix),1,mean), index=T)$ix
#sort(apply(apply( rsvm.res$featRank, 1,  function(x) sort(x, index.return=T)$ix),1,mean), index=T)$ix

# TODO
# 1. imbalanced sample size 
# 2. output sensitivity and specificy and other performance cretiria
# 3. to use mSVM-RFE or frequency SVM-RFE