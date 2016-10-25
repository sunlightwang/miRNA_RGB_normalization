# normalisation using lumi for PBMC miRNA expression data
# 
# file inputs: 
dataFile <- "miRNA rehash 17-1-2011_Group_Gene_Profile_raw_bg_sub.txt"
sampleInfoFile <- "miRNA rehash 17-1-2011 raw_bg_sub Samples Table.txt"
controlData <- "miRNA rehash 17-1-2011 raw_bg_sub Control Probe Profile.txt"

require(lumi)

# data reading

expr.lumi <- lumiR.batch(dataFile, sampleInfoFile = sampleInfoFile)
expr.lumi <- addControlData2lumi(controlData, expr.lumi)
expr.data <- exprs(expr.lumi)
ctrl.data <- controlData(expr.lumi)
U49.data <- as.double(ctrl.data[52,-(1:2), drop=FALSE]) # U49: 2011, row 52

# get microarray chip information
samples <- colnames(expr.data)
arrays <- sapply(samples, function(x) {
  t <- unlist(strsplit(x, "_"))
  t[1]
})
arrays <- as.factor(as.vector(arrays))
# get label (schiz vs. control) information 
label211 <- read.table("sampleLabel211.txt", header=T)
label <- as.factor(label211[match(samples, label211[,1]),4])
label[is.na(label)] <- "NA"

# orginal data plots
nsample <- ncol(expr.data)
plot(1:nsample, colMeans(expr.data), type='l', col="red", ylim=c(0,10000), main="Sample Mean&Stdev orginal")
points(1:nsample, sqrt(colMeans(expr.data^2) - colMeans(expr.data)^2), type='l', col="blue")
points(1:nsample, U49.data/3, type='l', col="green")
legend("topright", bty="n", bg="white", legend=c("Mean", "Stdev", "U49/3"), lty=1, col=c("red","blue","green"), cex=1.0)
#
boxplot(expr.lumi, main="Sample expression boxplot: orginal", col=label)
points(1:nsample, log(U49.data), type='l', col='blue')
legend("topright", bty="n", bg="white", legend=c("Ctrl","Schiz","U49"), lty=1, col=c("black","red","blue"), cex=1.0)

#----------------------- implementing reference gene normalisaton -----------------------------#
#  
# Strategy 1: simply divided by the reference gene expression value
# step 1. divided by U49
expr.data.dividedByU49 <- t( t(expr.data) / U49.data * mean(U49.data) )
expr.lumi.dividedByU49 <- expr.lumi
exprs(expr.lumi.dividedByU49) <- as.matrix(expr.data.dividedByU49)

nsample <- ncol(expr.data.dividedByU49)
plot(1:nsample, colMeans(expr.data.dividedByU49), type='l', col="red", ylim=c(1000,5000), main="Sample Mean&Stdev after divided by U49")
points(1:nsample, sqrt(colMeans(expr.data.dividedByU49^2) - colMeans(expr.data.dividedByU49)^2), type='l', col="blue")
legend("top", bty="n", bg="white", legend=c("Mean", "Stdev"), lty=1, col=c("red","blue"), cex=1.0)
plot(1:nsample, colMeans(expr.data.dividedByU49), pch=16, col=arrays, main="Sample Mean after divided by U49")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.dividedByU49), pch=16, col=label, main="Sample Mean after divided by U49")
title(main="color for labels", line=0.5)
#
plotSampleRelation(expr.lumi.dividedByU49, method='mds', color=arrays)
title(main="Sample MDS: after divided by U49 (color for chips)", line=0.5)
boxplot(expr.lumi.dividedByU49, main="Sample expression boxplot: after divided by U49")
plotDensity(expr.lumi.dividedByU49, col=arrays, main="Sample expression density: after divided by U49")

# step 2. log transformation 
expr.data.dividedByU49.lumiT <- log(expr.data.dividedByU49 - min(expr.data.dividedByU49) + 1.01 )
expr.lumi.dividedByU49.lumiT <- expr.lumi
exprs(expr.lumi.dividedByU49.lumiT) <- expr.data.dividedByU49.lumiT 

nsample <- ncol(expr.data.dividedByU49.lumiT)
plot(1:nsample, colMeans(expr.data.dividedByU49.lumiT), type='l', col="red", ylim=c(0,10), main="Sample Mean&Stdev after divided by U49 & logT")
points(1:nsample, sqrt(colMeans(expr.data.dividedByU49.lumiT^2) - colMeans(expr.data.dividedByU49.lumiT)^2), type='l', col="blue")
legend("topright", bty="n", bg="white", legend=c("Mean", "Stdev"), lty=1, col=c("red","blue"), cex=1.0)
plot(1:nsample, colMeans(expr.data.dividedByU49.lumiT), pch=16, col=arrays, ylim=c(7,9), main="Sample Mean after divided by U49 & logT")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.dividedByU49.lumiT), pch=16, col=label, ylim=c(7,9), main="Sample Mean after divided by U49 & logT")
title(main="color for labels", line=0.5)
#
plotSampleRelation(expr.lumi.dividedByU49.lumiT, method='mds', color=arrays)
title(main="Sample MDS: after divided by U49 & logT (color for chips)", line=0.5)
boxplot(expr.lumi.dividedByU49.lumiT, main="Sample expression boxplot: after divided by U49 & logT")
plotDensity(expr.lumi.dividedByU49.lumiT, col=arrays, xlim=c(2,4), main="Sample expression density: after divided by U49 & logT")

# Strategy 2: VST + substract refgene
# step 1. VST
expr.lumi.VST <- lumiT(expr.lumi) # default - VST
expr.data.VST <- exprs(expr.lumi.VST)
trans.para <- attr(expr.lumi.VST, "vstParameter")
trans.func <- attr(expr.lumi.VST, "transformFun")
U49.data.trans <- numeric(length(U49.data))
for(i in 1:ncol(expr.data.VST)) {
  U49.data.trans[i] <- ifelse(trans.func[i] == "asinh", 
                              trans.para[i,"g"] * asinh(trans.para[i,"a"] + trans.para[i,"b"] * U49.data[i]) + trans.para[i,"Intercept"], 
                              trans.para[i,"g"] * log(trans.para[i,"a"] + trans.para[i,"b"] * U49.data[i]) + trans.para[i,"Intercept"] )
}

nsample <- ncol(expr.data.VST)
plot(1:nsample, colMeans(expr.data.VST), type='l', col="red", ylim=c(0,20), main="Sample Mean&Stdev after VST")
points(1:nsample, sqrt(colMeans(expr.data.VST^2) - colMeans(expr.data.VST)^2), type='l', col="blue")
points(1:nsample, U49.data.trans, type='l', col="green")
legend("topright", bty="n", bg="white", legend=c("Mean", "Stdev", "U49"), lty=c(1,1,1), 
       col=c("red","blue","green"), cex=1.0)
plot(1:nsample, colMeans(expr.data.VST), pch=16, col=arrays, main="Sample Mean after VST")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.VST), pch=16, col=label, main="Sample Mean after VST")
title(main="color for labels", line=0.5)
#
plotSampleRelation(expr.lumi.VST, method='mds', color=arrays)
title(main="Sample MDS: after VST (color for chips)", line=0.5)
boxplot(expr.lumi.VST, main="Sample expression boxplot: after VST")
plotDensity(expr.lumi.VST, col=arrays, main="Sample expression density: after VST")

# step 2. substract transformed U49 expression values
expr.data.VST.U49 <- t( t(expr.data.VST) / U49.data.trans * mean(U49.data.trans))

nsample <- ncol(expr.data.VST.U49)
plot(1:nsample, colMeans(expr.data.VST.U49), type='l', lwd=2, col="red", ylim=c(-4,15), main="Sample Mean&Stdev after VST & U49 norm")
points(1:nsample, sqrt(colMeans(expr.data.VST.U49^2) - colMeans(expr.data.VST.U49)^2), type='l', col="blue")
points(1:nsample, colMeans(expr.data.VST), type='l', lty=1, col="green")
points(1:nsample, U49.data.trans, type='l', col="green", lty=2)
points(1:nsample, colMeans(expr.data.dividedByU49.lumiT), type='l', col="black", lty=2)
legend("topleft", bty="n", bg="white", legend=c("Mean", "Stdev","Mean VST","U49","Mean Str_1"), lty=c(1,1,1,2,2),  
       lwd=c(2,1,1,1,1), col=c("red","blue","green","black"), cex=1.0)
plot(1:nsample, colMeans(expr.data.VST.U49), pch=16, col=arrays, main="Sample Mean after VST & U49 norm")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.VST.U49), pch=16, col=label, main="Sample Mean after VST & U49 norm")
title(main="color for labels", line=0.5)
#
expr.lumi.VST.U49 <- expr.lumi
exprs(expr.lumi.VST.U49) <- as.matrix(expr.data.VST.U49)
plotSampleRelation(expr.lumi.VST.U49, method='mds', color=arrays)
title(main="Sample MDS: after VST & U49 norm (color for chips)", line=0.5)
boxplot(expr.lumi.VST.U49, main="Sample expression boxplot: after VST & U49 norm")
plotDensity(expr.lumi.VST.U49, col=arrays, main="Sample expression density: after VST & U49 norm")

### End of normalisation by U49

######### using Erin's reference genes #############
# U24 & U49 (non U66 according to Murray)
# U66: 3792, row 11
# U24: 1992, row 9
U66.data <- as.double(ctrl.data[11,-(1:2), drop=FALSE])
U24.data <- as.double(ctrl.data[9,-(1:2), drop=FALSE])
U24U49.data <- sqrt(U49.data * U24.data)

# ------- also using two strategies as above ------------------------
# Strategy 1: simply divided by the reference gene expression value
# step 1. divided by U49
expr.data.dividedByU24U49 <- t( t(expr.data) / U24U49.data * mean(U24U49.data) )
expr.lumi.dividedByU24U49 <- expr.lumi
exprs(expr.lumi.dividedByU24U49) <- as.matrix(expr.data.dividedByU24U49)

nsample <- ncol(expr.data.dividedByU24U49)
plot(1:nsample, colMeans(expr.data.dividedByU24U49), type='l', col="red", ylim=c(0,20000), main="Sample Mean&Stdev after divided by U24U49")
points(1:nsample, sqrt(colMeans(expr.data.dividedByU24U49^2) - colMeans(expr.data.dividedByU24U49)^2), type='l', col="blue")
legend("top", bty="n", bg="white", legend=c("Mean", "Stdev"), lty=1, col=c("red","blue"), cex=1.0)
plot(1:nsample, colMeans(expr.data.dividedByU24U49), pch=16, col=arrays, main="Sample Mean after divided by U24U49")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.dividedByU24U49), pch=16, col=label, main="Sample Mean after divided by U24U49")
title(main="color for labels", line=0.5)
#
plotSampleRelation(expr.lumi.dividedByU24U49, method='mds', color=arrays)
title(main="Sample MDS: after divided by U24U49 (color for chips)", line=0.5)
boxplot(expr.lumi.dividedByU24U49, main="Sample expression boxplot: after divided by U24U49")
plotDensity(expr.lumi.dividedByU24U49, col=arrays, main="Sample expression density: after divided by U24U49")

# step 2. log transformation 
expr.data.dividedByU24U49.lumiT <- log(expr.data.dividedByU24U49 - min(expr.data.dividedByU24U49) + 1.01 )
expr.lumi.dividedByU24U49.lumiT <- expr.lumi
exprs(expr.lumi.dividedByU24U49.lumiT) <- expr.data.dividedByU24U49.lumiT 

nsample <- ncol(expr.data.dividedByU24U49.lumiT)
plot(1:nsample, colMeans(expr.data.dividedByU24U49.lumiT), type='l', col="red", ylim=c(0,10), main="Sample Mean&Stdev after divided by U24U49 & logT")
points(1:nsample, sqrt(colMeans(expr.data.dividedByU24U49.lumiT^2) - colMeans(expr.data.dividedByU24U49.lumiT)^2), type='l', col="blue")
legend("topright", bty="n", bg="white", legend=c("Mean", "Stdev"), lty=1, col=c("red","blue"), cex=1.0)
plot(1:nsample, colMeans(expr.data.dividedByU24U49.lumiT), pch=16, col=arrays, ylim=c(7,9), main="Sample Mean after divided by U24U49 & logT")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.dividedByU24U49.lumiT), pch=16, col=label, ylim=c(7,9), main="Sample Mean after divided by U24U49 & logT")
title(main="color for labels", line=0.5)
#
plotSampleRelation(expr.lumi.dividedByU24U49.lumiT, method='mds', color=arrays)
title(main="Sample MDS: after divided by U24U49 & logT (color for chips)", line=0.5)
boxplot(expr.lumi.dividedByU24U49.lumiT, main="Sample expression boxplot: after divided by U24U49 & logT")
plotDensity(expr.lumi.dividedByU24U49.lumiT, col=arrays, xlim=c(2,4), main="Sample expression density: after divided by U24U49 & logT")

# Strategy 2: VST + substract refgene
# step 1. VST
expr.lumi.VST <- lumiT(expr.lumi) # default - VST
expr.data.VST <- exprs(expr.lumi.VST)
trans.para <- attr(expr.lumi.VST, "vstParameter")
trans.func <- attr(expr.lumi.VST, "transformFun")
U24U49.data.trans <- numeric(length(U24U49.data))
for(i in 1:ncol(expr.data.VST)) {
  U24U49.data.trans[i] <- ifelse(trans.func[i] == "asinh", 
                              trans.para[i,"g"] * asinh(trans.para[i,"a"] + trans.para[i,"b"] * U24U49.data[i]) + trans.para[i,"Intercept"], 
                              trans.para[i,"g"] * log(trans.para[i,"a"] + trans.para[i,"b"] * U24U49.data[i]) + trans.para[i,"Intercept"] )
}

nsample <- ncol(expr.data.VST)
plot(1:nsample, colMeans(expr.data.VST), type='l', col="red", ylim=c(0,20), main="Sample Mean&Stdev after VST")
points(1:nsample, sqrt(colMeans(expr.data.VST^2) - colMeans(expr.data.VST)^2), type='l', col="blue")
points(1:nsample, U24U49.data.trans, type='l', col="green")
legend("topright", bty="n", bg="white", legend=c("Mean", "Stdev", "U24U49"), lty=c(1,1,1), 
       col=c("red","blue","green"), cex=1.0)
plot(1:nsample, colMeans(expr.data.VST), pch=16, col=arrays, main="Sample Mean after VST")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.VST), pch=16, col=label, main="Sample Mean after VST")
title(main="color for labels", line=0.5)
#
plotSampleRelation(expr.lumi.VST, method='mds', color=arrays)
title(main="Sample MDS: after VST (color for chips)", line=0.5)
boxplot(expr.lumi.VST, main="Sample expression boxplot: after VST")
plotDensity(expr.lumi.VST, col=arrays, main="Sample expression density: after VST")

# step 2. substract transformed U24U49 expression values
expr.data.VST.U24U49 <- t( t(expr.data.VST) / U24U49.data.trans * mean(U24U49.data.trans))

nsample <- ncol(expr.data.VST.U24U49)
plot(1:nsample, colMeans(expr.data.VST.U24U49), type='l', lwd=2, col="red", ylim=c(0,18), main="Sample Mean&Stdev after VST & U24U49 norm")
points(1:nsample, sqrt(colMeans(expr.data.VST.U24U49^2) - colMeans(expr.data.VST.U24U49)^2), type='l', col="blue")
points(1:nsample, colMeans(expr.data.VST), type='l', lty=1, col="green")
points(1:nsample, U24U49.data.trans, type='l', col="green", lty=2)
points(1:nsample, colMeans(expr.data.dividedByU24U49.lumiT), type='l', col="black", lty=2)
legend("topleft", bty="n", bg="white", legend=c("Mean", "Stdev","Mean VST","U24U49","Mean Str_1"), lty=c(1,1,1,2,2),  
       lwd=c(2,1,1,1,1), col=c("red","blue","green","black"), cex=1.0)
plot(1:nsample, colMeans(expr.data.VST.U24U49), pch=16, col=arrays, main="Sample Mean after VST & U24U49 norm")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.VST.U24U49), pch=16, col=label, main="Sample Mean after VST & U24U49 norm")
title(main="color for labels", line=0.5)
#
expr.lumi.VST.U24U49 <- expr.lumi
exprs(expr.lumi.VST.U24U49) <- as.matrix(expr.data.VST.U24U49)
plotSampleRelation(expr.lumi.VST.U24U49, method='mds', color=arrays)
title(main="Sample MDS: after VST & U24U49 norm (color for chips)", line=0.5)
boxplot(expr.lumi.VST.U24U49, main="Sample expression boxplot: after VST & U24U49 norm")
plotDensity(expr.lumi.VST.U24U49, col=arrays, main="Sample expression density: after VST & U24U49 norm")

### End of normalisation by U24 & U49


######### using U24 as the only reference gene #############
# U24 & U49 (non U66 according to Murray)
# U66: 3792, row 11
# U24: 1992, row 9
#U66.data <- as.double(ctrl.data[11,-(1:2), drop=FALSE])
U24.data <- as.double(ctrl.data[9,-(1:2), drop=FALSE])

# ------- also using two strategies as above ------------------------
# Strategy 1: simply divided by the reference gene expression value
# step 1. divided by U49
expr.data.dividedByU24 <- t( t(expr.data) / U24.data * mean(U24.data) )
expr.lumi.dividedByU24 <- expr.lumi
exprs(expr.lumi.dividedByU24) <- as.matrix(expr.data.dividedByU24)

nsample <- ncol(expr.data.dividedByU24)
plot(1:nsample, colMeans(expr.data.dividedByU24), type='l', col="red", ylim=c(0,20000), main="Sample Mean&Stdev after divided by U24")
points(1:nsample, sqrt(colMeans(expr.data.dividedByU24^2) - colMeans(expr.data.dividedByU24)^2), type='l', col="blue")
legend("top", bty="n", bg="white", legend=c("Mean", "Stdev"), lty=1, col=c("red","blue"), cex=1.0)
plot(1:nsample, colMeans(expr.data.dividedByU24), pch=16, col=arrays, main="Sample Mean after divided by U24")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.dividedByU24), pch=16, col=label, main="Sample Mean after divided by U24")
title(main="color for labels", line=0.5)
#
plotSampleRelation(expr.lumi.dividedByU24, method='mds', color=arrays)
title(main="Sample MDS: after divided by U24 (color for chips)", line=0.5)
boxplot(expr.lumi.dividedByU24, main="Sample expression boxplot: after divided by U24")
plotDensity(expr.lumi.dividedByU24, col=arrays, main="Sample expression density: after divided by U24")

# step 2. log transformation 
expr.data.dividedByU24.lumiT <- log(expr.data.dividedByU24 - min(expr.data.dividedByU24) + 1.01 )
expr.lumi.dividedByU24.lumiT <- expr.lumi
exprs(expr.lumi.dividedByU24.lumiT) <- expr.data.dividedByU24.lumiT 

nsample <- ncol(expr.data.dividedByU24.lumiT)
plot(1:nsample, colMeans(expr.data.dividedByU24.lumiT), type='l', col="red", ylim=c(0,10), main="Sample Mean&Stdev after divided by U24 & logT")
points(1:nsample, sqrt(colMeans(expr.data.dividedByU24.lumiT^2) - colMeans(expr.data.dividedByU24.lumiT)^2), type='l', col="blue")
legend("topright", bty="n", bg="white", legend=c("Mean", "Stdev"), lty=1, col=c("red","blue"), cex=1.0)
plot(1:nsample, colMeans(expr.data.dividedByU24.lumiT), pch=16, col=arrays, ylim=c(7,9), main="Sample Mean after divided by U24 & logT")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.dividedByU24.lumiT), pch=16, col=label, ylim=c(7,9), main="Sample Mean after divided by U24 & logT")
title(main="color for labels", line=0.5)
#
plotSampleRelation(expr.lumi.dividedByU24.lumiT, method='mds', color=arrays)
title(main="Sample MDS: after divided by U24 & logT (color for chips)", line=0.5)
boxplot(expr.lumi.dividedByU24.lumiT, main="Sample expression boxplot: after divided by U24 & logT")
plotDensity(expr.lumi.dividedByU24.lumiT, col=arrays, xlim=c(2,4), main="Sample expression density: after divided by U24 & logT")

# Strategy 2: VST + substract refgene
# step 1. VST
expr.lumi.VST <- lumiT(expr.lumi) # default - VST
expr.data.VST <- exprs(expr.lumi.VST)
trans.para <- attr(expr.lumi.VST, "vstParameter")
trans.func <- attr(expr.lumi.VST, "transformFun")
U24.data.trans <- numeric(length(U24.data))
for(i in 1:ncol(expr.data.VST)) {
  U24.data.trans[i] <- ifelse(trans.func[i] == "asinh", 
                                 trans.para[i,"g"] * asinh(trans.para[i,"a"] + trans.para[i,"b"] * U24.data[i]) + trans.para[i,"Intercept"], 
                                 trans.para[i,"g"] * log(trans.para[i,"a"] + trans.para[i,"b"] * U24.data[i]) + trans.para[i,"Intercept"] )
}

nsample <- ncol(expr.data.VST)
plot(1:nsample, colMeans(expr.data.VST), type='l', col="red", ylim=c(0,20), main="Sample Mean&Stdev after VST")
points(1:nsample, sqrt(colMeans(expr.data.VST^2) - colMeans(expr.data.VST)^2), type='l', col="blue")
points(1:nsample, U24.data.trans, type='l', col="green")
legend("topright", bty="n", bg="white", legend=c("Mean", "Stdev", "U24"), lty=c(1,1,1), 
       col=c("red","blue","green"), cex=1.0)
plot(1:nsample, colMeans(expr.data.VST), pch=16, col=arrays, main="Sample Mean after VST")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.VST), pch=16, col=label, main="Sample Mean after VST")
title(main="color for labels", line=0.5)
#
plotSampleRelation(expr.lumi.VST, method='mds', color=arrays)
title(main="Sample MDS: after VST (color for chips)", line=0.5)
boxplot(expr.lumi.VST, main="Sample expression boxplot: after VST")
plotDensity(expr.lumi.VST, col=arrays, main="Sample expression density: after VST")

# step 2. substract transformed U24 expression values
expr.data.VST.U24 <- t( t(expr.data.VST) / U24.data.trans * mean(U24.data.trans))

nsample <- ncol(expr.data.VST.U24)
plot(1:nsample, colMeans(expr.data.VST.U24), type='l', lwd=2, col="red", ylim=c(0,18), main="Sample Mean&Stdev after VST & U24 norm")
points(1:nsample, sqrt(colMeans(expr.data.VST.U24^2) - colMeans(expr.data.VST.U24)^2), type='l', col="blue")
points(1:nsample, colMeans(expr.data.VST), type='l', lty=1, col="green")
points(1:nsample, U24.data.trans, type='l', col="green", lty=2)
points(1:nsample, colMeans(expr.data.dividedByU24.lumiT), type='l', col="black", lty=2)
legend("topleft", bty="n", bg="white", legend=c("Mean", "Stdev","Mean VST","U24","Mean Str_1"), lty=c(1,1,1,2,2),  
       lwd=c(2,1,1,1,1), col=c("red","blue","green","black"), cex=1.0)
plot(1:nsample, colMeans(expr.data.VST.U24), pch=16, col=arrays, main="Sample Mean after VST & U24 norm")
title(main="color for chips", line=0.5)
plot(1:nsample, colMeans(expr.data.VST.U24), pch=16, col=label, main="Sample Mean after VST & U24 norm")
title(main="color for labels", line=0.5)
#
expr.lumi.VST.U24 <- expr.lumi
exprs(expr.lumi.VST.U24) <- as.matrix(expr.data.VST.U24)
plotSampleRelation(expr.lumi.VST.U24, method='mds', color=arrays)
title(main="Sample MDS: after VST & U24 norm (color for chips)", line=0.5)
boxplot(expr.lumi.VST.U24, main="Sample expression boxplot: after VST & U24 norm")
plotDensity(expr.lumi.VST.U24, col=arrays, main="Sample expression density: after VST & U24 norm")

### End of normalisation by U24 ###
