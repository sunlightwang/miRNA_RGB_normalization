# normalisation using lumi for PBMC miRNA expression data
# 
# file inputs: 
dataFile <- "miRNA rehash 17-1-2011_Group_Gene_Profile_raw_bg_sub.txt"
sampleInfoFile <- "miRNA rehash 17-1-2011 raw_bg_sub Samples Table.txt"
controlData <- "miRNA rehash 17-1-2011 raw_bg_sub Control Probe Profile.txt"

require(lumi)

# (1) lumiR: data reading
expr.lumi <- lumiR.batch(dataFile, sampleInfoFile = sampleInfoFile)
expr.lumi <- addControlData2lumi(controlData, expr.lumi)
expr.data <- exprs(expr.lumi)

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

# plot regarding control data 
plotControlData(expr.lumi, main="Control data boxplot: original")
#plotHousekeepingGene(expr.lumi)
#plotStringencyGene(expr.lumi)
#******* before VST and normalisation MDS plot *************
plotSampleRelation(expr.lumi, method='mds', color=arrays) # arrays is indicator of array numbers, from normalisation1.R
title(main = "Sample MDS: original data (color for chips)", line=0.5)
plotSampleRelation(expr.lumi, method='mds', color=label) # arrays is indicator of array numbers, from normalisation1.R
title(main = "Sample MDS: original data (color for labels)", line=0.5)
boxplot(expr.lumi, main="Sample expression boxplot: original data")
plotDensity(expr.lumi, col=arrays, main="Sample expression density: original data")

# (2) lumiB: backgroup substraction 
# expr.lumi <- lumiB(expr.lumi) # none will be done, as the bg has already been substracted 

# (3) lumiT: transformation 
expr.lumi <- lumiT(expr.lumi) # default - VST
#******* before normalisation MDS plot *************
plotSampleRelation(expr.lumi, method='mds', color=arrays)
title(main = "Sample MDS: after VS transformation (color for chips)", line=0.5)
boxplot(expr.lumi, main="Sample expression boxplot: after VS transformation")
plotDensity(expr.lumi, col=arrays, main="Sample expression density: after VST")

# (4) limuN: normaliation
expr.lumi.qn <- lumiN(expr.lumi, method='quantile') # default: quantile
#******* after normalisation MDS plot *************
plotSampleRelation(expr.lumi.qn, method='mds', color=arrays)
title(main = "Sample MDS: after quantile normalisation (color for chips)", line=0.5)
# still can observe batch effects
boxplot(expr.lumi.qn, main="Sample expression boxplot: after quantile normalisation")
plotDensity(expr.lumi.qn, col=arrays, main="Sample expression density: after quantile normalisation")

expr.lumi.rsn <- lumiN(expr.lumi, method="rsn") # Robust Spline Normalization
#******* after normalisation MDS plot *************
plotSampleRelation(expr.lumi.rsn, method='mds', color=arrays)
title(main = "Sample MDS: after Robust Spline Normalization (color for chips)", line=0.5)
# still batch effect 
boxplot(expr.lumi.rsn, main="Sample expression boxplot: after Robust Spline Normalization")
plotDensity(expr.lumi.rsn, col=arrays, main="Sample expression density: after Robust Spline Normalization")

expr.lumi.ssn <- lumiN(expr.lumi, method="ssn") # Simple Scaling Normalization
#******* after normalisation MDS plot *************
plotSampleRelation(expr.lumi.ssn, method='mds', color=arrays)
title(main = "Sample MDS: after Simple Scaling Normalization (color for chips)", line=0.5)
boxplot(expr.lumi.ssn, main="Sample expression boxplot: after Simple Scaling Normalization")
plotDensity(expr.lumi.ssn, col=arrays, xlim=c(2,6), main="Sample expression density: after Simple Scaling Normalization")

expr.lumi.vsn <- lumiN(expr.lumi, method="vsn") # Variance stabilization
#******* after normalisation MDS plot *************
plotSampleRelation(expr.lumi.vsn, method='mds', color=arrays)
title(main = "Sample MDS: after Variance Stabilization Normalisation (color for chips)", line=0.5)
boxplot(expr.lumi.vsn, main="Sample expression boxplot: after Variance Stabilization Normalisation")
plotDensity(expr.lumi.vsn, col=arrays, main="Sample expression density: after Variance Stabilization Normalisation")

expr.lumi.rin <- lumiN(expr.lumi, method="rankinvariant") # Rank Invariant Normalization
#******* after normalisation MDS plot *************
plotSampleRelation(expr.lumi.rin, method='mds', color=arrays)
title(main = "Sample MDS: after Rank Invariant Normalization (color for chips)", line=0.5)
boxplot(expr.lumi.rin, main="Sample expression boxplot: after Rank Invariant Normalization")
plotDensity(expr.lumi.rin, col=arrays, xlim=c(2,6), main="Sample expression density: after Rank Invariant Normalization")

expr.lumi.loess <- lumiN(expr.lumi, method="loess") # Loess Normalization
#******* after normalisation MDS plot *************
plotSampleRelation(expr.lumi.loess, method='mds', color=arrays)
title(main = "Sample MDS: after Loess Normalization (color for chips)", line=0.5)
boxplot(expr.lumi.loess, main="Sample expression boxplot: after Loess Normalization")
plotDensity(expr.lumi.loess, col=arrays, main="Sample expression density: after Loess Normalization")

# (4) limuQ: QC
expr.lumi <- lumiQ(expr.lumi)
# and the results by different normalisation methods 

