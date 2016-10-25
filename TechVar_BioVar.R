# file inputs: 
setwd("~/Projects/bioMarker/PBMCFromRawData")

dataFile <- "miRNA rehash 17-1-2011_Group_Gene_Profile_raw_bg_sub.txt"
sampleInfoFile <- "miRNA rehash 17-1-2011 raw_bg_sub Samples Table.txt"
controlData <- "miRNA rehash 17-1-2011 raw_bg_sub Control Probe Profile.txt"

require(lumi)

# (1) lumiR: data reading
expr.lumi <- lumiR.batch(dataFile, sampleInfoFile = sampleInfoFile)
expr.data <- exprs(expr.lumi)
se.expr.data <- se.exprs(expr.lumi)
miRNAs <- row.names(expr.data)

# (2) log transformation 
expr.lumi <- lumiT(expr.lumi, 'log2')
offset <- apply(expr.data, 2, min, na.rm = TRUE)
offset[offset <= 0] <- offset[offset <= 0] - 1.01
offset[offset > 0] <- 0
expr.data.T <- exprs(expr.lumi)

# (3) calculate tech var and bio var ## considering cal bio var after lumiN
se.expr.data.T <- se.expr.data / (expr.data - rep(1, nrow(expr.data)) %*% t(offset) )
bioSd = sqrt( rowMeans(expr.data.T ^2) - rowMeans(expr.data.T) ^ 2 )
# bioSe = sqrt( rowMeans(expr.data ^2) - rowMeans(expr.data) ^ 2 ) / sqrt(ncol(expr.data))
techSeMean = rowMeans(se.expr.data.T)
# techSeSd = sqrt( rowMeans(se.expr.data ^ 2 ) - techSeMean^2 )

plot(bioSd, techSeMean)
# points(bioSd, techSeSd, col='red')
# 
# t <- lm(techSeMean ~ bioSd) 
# abline(t$coefficients[1], t$coefficients[2], col='red')

# There are some genes/miRNAs with bioVar relatively larger than techVar, probably only these genes can be detected as DEGs. 

is.outlier <- rep(0, length(bioSd))
max.iter <- 20
distance <- 2000
for (i in 1:max.iter) {
  t <- lm(techSeMean[is.outlier==0] ~ bioSd[is.outlier==0]) 
  is.outlier <- (t$coefficients[1] + t$coefficients[2] * bioSd - 
    techSeMean) > (distance / max.iter * 2)
}
plot(bioSd, techSeMean, pch=16, col=as.numeric(is.outlier)+1)

miRNAs.dist <- as.data.frame (t$coefficients[1] + t$coefficients[2] * bioSd - techSeMean)
miRNAs.dist.sorted <- as.data.frame(sort (t$coefficients[1] + t$coefficients[2] * bioSd - techSeMean, decreasing=TRUE))
write.table(miRNAs.dist[is.outlier,,drop=FALSE], "miRNAs.bigBioVar.txt", col.names=FALSE, quote=FALSE, sep="\t")
write.table(miRNAs.dist.sorted, "miRNAs.bioVar2techVarDist.txt", col.names=FALSE, quote=FALSE, sep="\t")
miRNAs.bigBioVar <- miRNAs[is.outlier]
