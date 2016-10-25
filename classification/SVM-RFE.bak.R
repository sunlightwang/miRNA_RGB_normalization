# implementation of SVM-RFE
library(e1071)

################################################
# Feature Ranking with SVM-RFE
################################################
svmrfeFeatureRanking = function(x,y, rsvm=FALSE){
  n = ncol(x)
  
  survivingFeaturesIndexes = seq(1:n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n
  
  while(length(survivingFeaturesIndexes)>0){
    #train the support vector machine
    svmModel = svm(x[, survivingFeaturesIndexes], y, cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" )
    
    if(rsvm) {
      #yy <- vector(length = length(y)) 
      #yy[y==levels(y)[1]] <- 1
      #yy[y==levels(y)[2]] <- -1
      md <- colMeans(x[y==levels(y)[1],]) - colMeans(x[y==levels(y)[2], ])
      #rankingCriteria <- t(svmModel$coefs*yy[svmModel$index]) %*% svmModel$SV * md
      rankingCriteria <- t(svmModel$coefs) %*% svmModel$SV * md
    } else {
      #compute the weight vector
      w = t(svmModel$coefs)%*%svmModel$SV
      #compute ranking criteria
      rankingCriteria = w * w
    }
    
    #rank the features
    ranking = sort(rankingCriteria, index.return = TRUE)$ix
    
    #update feature ranked list
    featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]]
    rankedFeatureIndex = rankedFeatureIndex - 1
    
    #eliminate the feature with smallest ranking criterion
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
    
  }
  
  return (featureRankedList)
}


################################################
# Feature Ranking with Average Multiclass SVM-RFE
################################################

svmrfeFeatureRankingForMulticlass = function(x,y){
  n = ncol(x)
  
  survivingFeaturesIndexes = seq(1:n)
  featureRankedList = vector(length=n)
  rankedFeatureIndex = n
  
  while(length(survivingFeaturesIndexes)>0){
    #train the support vector machine
    svmModel = svm(x[, survivingFeaturesIndexes], y, cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" )
    
    #compute the weight vector
    multiclassWeights = svm.weights(svmModel)
    
    #compute ranking criteria
    multiclassWeights = multiclassWeights * multiclassWeights
    rankingCriteria = 0
    for(i in 1:ncol(multiclassWeights))rankingCriteria[i] = mean(multiclassWeights[,i])
    
    #rank the features
    (ranking = sort(rankingCriteria, index.return = TRUE)$ix)
    
    #update feature ranked list
    (featureRankedList[rankedFeatureIndex] = survivingFeaturesIndexes[ranking[1]])
    rankedFeatureIndex = rankedFeatureIndex - 1
    
    #eliminate the feature with smallest ranking criterion
    (survivingFeaturesIndexes = survivingFeaturesIndexes[-ranking[1]])
    cat(length(survivingFeaturesIndexes),"\n")
  }
  
}

################################################
# This function gives the weights of the hiperplane
################################################
svm.weights<-function(model){
  w=0
  if(model$nclasses==2){
    w=t(model$coefs)%*%model$SV
  }else{    #when we deal with OVO svm classification
    ## compute start-index
    start <- c(1, cumsum(model$nSV)+1)
    start <- start[-length(start)]
    
    calcw <- function (i,j) {
      ## ranges for class i and j:
      ri <- start[i] : (start[i] + model$nSV[i] - 1)
      rj <- start[j] : (start[j] + model$nSV[j] - 1)
      
      ## coefs for (i,j):
      coef1 <- model$coefs[ri, j-1]
      coef2 <- model$coefs[rj, i]
      ## return w values:
      w=t(coef1)%*%model$SV[ri,]+t(coef2)%*%model$SV[rj,]
      return(w)
    }
    
    W=NULL
    for (i in 1 : (model$nclasses - 1)){
      for (j in (i + 1) : model$nclasses){
        wi=calcw(i,j)
        W=rbind(W,wi)
      }
    }
    w=W
  }
  return(w)
}


##################
# Testing 
##################
require(caret)
data(mdrr)
x <- as.matrix(mdrrDescr[,-nearZeroVar(mdrrDescr)])
y <- as.factor(mdrrClass)

#scale samples with mean zero and standard deviation one
for(i in 1:nrow(x))x[i,] = (x[i,]-mean(x[i,]))/sd(x[i,])
#scale features with mean zero and standard deviation one
for(i in 1:ncol(x))x[,i] = (x[,i]-mean(x[,i]))/sd(x[,i])
x = 2*atan(x/2)

# Feature Ranking with SVM-RFE
featureRankedList = svmrfeFeatureRanking(x,y)

# Leave One Out test using the best n ranked features of the SVM-RFE
nFeat <- length(featureRankedList)
steps <- unique(as.integer(exp(seq(0, log(nFeat), 0.1))))
nSample <- length(y)
for(nfeatures in steps){
  truePredictions = 0
  for(i in 1:nrow(x)){
    #p <- tune.svm(x[-i, featureRankedList[1:nfeatures]], y[-i], cost = 2^(-1:6), scale=F, cross=10, type="C-classification", kernel="linear",  cachesize=500)
    #cost <- as.numeric(p$best.parameters)
    svmModel = svm(x[-i, featureRankedList[1:nfeatures]], y[-i], cost = 10, cachesize=500,  scale=F, type="C-classification", kernel="linear" ) 
    prediction = predict(svmModel,matrix(x[i, featureRankedList[1:nfeatures]],nrow=1))
    if(prediction[[1]] == y[[i]]) truePredictions = truePredictions + 1
  }
  cat(nfeatures,":",truePredictions/nrow(x),"\n")
}

# TODO
# 1. imbalanced sample size 
# 2. output sensitivity and specificy and other performance cretiria
# 3. to use mSVM-RFE or frequency SVM-RFE