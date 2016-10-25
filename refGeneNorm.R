library(lumi)
source("lumi_edited.R")

refGeneNorm <- function(x.lumi, refGeneCol, refGeneName, method=c("vst","log2")) {
  stopifnot (is(x.lumi, "eSet")) 
  method <- match.arg(method, c("vst","log2"))
  expr.data <- exprs(x.lumi)
  
  # ref gene expression 
  ctrl.data <- controlData(x.lumi)
  refGene <- as.matrix(ctrl.data[refGeneCol,-(1:2), drop=FALSE])
  if(nrow(refGene) > 1) { # using geometric mean
    refGene <- exp(colMeans(log(refGene)))
  } else {
    refGene <- as.double(refGene)
  }
    
  refGeneName <- paste(refGeneName,collapse="")
  
  # boxplot of original value
  boxplot(x.lumi, main="Sample expression boxplot: orginal")
  
  # methods: log and vst
  if(method=="log2") {  # Strategy 1: reference gene + log2
    # step 1. divided by refGeme expr
    expr.data.divided <- t( t(expr.data) / refGene * mean(refGene) )
    expr.lumi.divided <- x.lumi
    exprs(expr.lumi.divided) <- as.matrix(expr.data.divided)
    
    nsample <- ncol(expr.data.divided)
    plot(1:nsample, colMeans(expr.data.divided), type='l', col="red", ylim=c(1000,5000), main=paste("Sample Mean&Stdev after divided by",refGeneName))
    points(1:nsample, sqrt(colMeans(expr.data.divided^2) - colMeans(expr.data.divided)^2), type='l', col="blue")
    legend("top", bty="n", bg="white", legend=c("Mean", "Stdev"), lty=1, col=c("red","blue"), cex=1.0)
    #
    boxplot(expr.lumi.divided, main=paste("Sample expression boxplot: after divided by", refGeneName))
    plotDensity(expr.lumi.divided, main=paste("Sample expression density: after divided by", refGeneName))
    
    # step 2. log transformation 
    expr.data.divided.log <- log(expr.data.divided - min(expr.data.divided) + 1.01 )
    expr.lumi.divided.log <- x.lumi
    exprs(expr.lumi.divided.log) <- expr.data.divided.log 
    
    nsample <- ncol(expr.data.divided.log)
    plot(1:nsample, colMeans(expr.data.divided.log), type='l', col="red", ylim=c(0,10), main=paste("Sample Mean&Stdev after divided by",refGeneName,"& logT"))
    points(1:nsample, sqrt(colMeans(expr.data.divided.log^2) - colMeans(expr.data.divided.log)^2), type='l', col="blue")
    legend("topright", bty="n", bg="white", legend=c("Mean", "Stdev"), lty=1, col=c("red","blue"), cex=1.0)
    #
    boxplot(expr.lumi.divided.log, main=paste("Sample expression boxplot: after divided by",refGeneName, "& logT"))
    plotDensity(expr.lumi.divided.log, xlim=c(2,4), main=paste("Sample expression density: after divided by",refGeneName,"& logT"))
    
    return(expr.lumi.divided.log)
  } else {   # Strategy 2: VST + refgene
    expr.lumi.VST <- lumiT(x.lumi) 
    expr.data.VST <- exprs(expr.lumi.VST)
    trans.para <- attr(expr.lumi.VST, "vstParameter")
    trans.func <- attr(expr.lumi.VST, "transformFun")
    refGene.trans <- numeric(length(refGene))
    for(i in 1:ncol(expr.data.VST)) {
      refGene.trans[i] <- ifelse(trans.func[i] == "asinh", 
                                 trans.para[i,"g"] * asinh(trans.para[i,"a"] + trans.para[i,"b"] * refGene[i]) + trans.para[i,"Intercept"], 
                                 trans.para[i,"g"] * log(trans.para[i,"a"] + trans.para[i,"b"] * refGene[i]) + trans.para[i,"Intercept"] )
    }
    
    nsample <- ncol(expr.data.VST)
    plot(1:nsample, colMeans(expr.data.VST), type='l', col="red", ylim=c(0,20), main="Sample Mean&Stdev after VST")
    points(1:nsample, sqrt(colMeans(expr.data.VST^2) - colMeans(expr.data.VST)^2), type='l', col="blue")
    points(1:nsample, refGene.trans, type='l', col="green")
    legend("topright", bty="n", bg="white", legend=c("Mean", "Stdev", "RefGene"), lty=c(1,1,1), 
           col=c("red","blue","green"), cex=1.0)
    #
    boxplot(expr.lumi.VST, main="Sample expression boxplot: after VST")
    plotDensity(expr.lumi.VST, main="Sample expression density: after VST")
    
    # step 2. substract transformed U49 expression values
    expr.data.VST.refGene <- t( t(expr.data.VST) / refGene.trans * mean(refGene.trans))
    expr.lumi.VST.refGene <- x.lumi
    exprs(expr.lumi.VST.refGene) <- as.matrix(expr.data.VST.refGene)
    
    nsample <- ncol(expr.data.VST.refGene)
    plot(1:nsample, colMeans(expr.data.VST.refGene), type='l', lwd=2, col="red", ylim=c(-4,15), main=paste("Sample Mean&Stdev after VST &",refGeneName, "norm"))
    points(1:nsample, sqrt(colMeans(expr.data.VST.refGene^2) - colMeans(expr.data.VST.refGene)^2), type='l', col="blue")
    points(1:nsample, colMeans(expr.data.VST), type='l', lty=1, col="green")
    points(1:nsample, refGene.trans, type='l', col="black", lty=2)
    legend("topleft", bty="n", bg="white", legend=c("Mean", "Stdev","Mean VST",refGeneName), lty=c(1,1,1,2),  
           lwd=c(2,1,1,1), col=c("red","blue","green","black"), cex=1.0)
    #
    boxplot(expr.lumi.VST.refGene, main=paste("Sample expression boxplot: after VST &", refGeneName, "norm"))
    plotDensity(expr.lumi.VST.refGene, main=paste("Sample expression density: after VST &", refGeneName, "norm"))
    
    return(expr.lumi.VST.refGene)
  }
}

## main ##
dataFile <- "miRNA rehash 17-1-2011_Group_Gene_Profile_raw_bg_sub.txt"
sampleInfoFile <- "miRNA rehash 17-1-2011 raw_bg_sub Samples Table.txt"
controlDataFile <- "miRNA rehash 17-1-2011 raw_bg_sub Control Probe Profile.txt"
expr.lumi <- lumiR.batch(dataFile, sampleInfoFile = sampleInfoFile)
expr.lumi <- addControlData2lumi(controlDataFile, expr.lumi)
expr.vst.U24 <- refGeneNorm(expr.lumi, 9, "U24", method="vst")

write.table(exprs(expr.vst.U24), "PBMC.miRNA.VST_U24.norm.txt", quote=FALSE, sep="\t")

