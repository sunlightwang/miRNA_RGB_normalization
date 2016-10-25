setwd("/Users/xw364/Projects/bioMarker/PBMCFromRawData/classification")

# normaliation of PMBC data
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


# main 
require(lumi)
require(sva)
dataFile <- "../miRNA rehash 17-1-2011_Group_Gene_Profile_raw_bg_sub.txt"
sampleInfoFile <- "../miRNA rehash 17-1-2011 raw_bg_sub Samples Table.txt"
controlDataFile <- "../miRNA rehash 17-1-2011 raw_bg_sub Control Probe Profile.txt"
expr.lumi <- lumiR.batch(dataFile, sampleInfoFile = sampleInfoFile)
expr.lumi <- addControlData2lumi(controlDataFile, expr.lumi)
labelFile <- "../sampleLabel211.txt"
labelData <- read.table(labelFile, header=T)[,c(1,4)]

subtypeFile <- "../sample211Subtype.txt"
subtypeData <- read.table(subtypeFile, header=T)[,c(2,6,7)]
GoM3Type <- as.factor(subtypeData[match(samples, subtypeData[,1]),2])
cogType <- as.factor(subtypeData[match(samples, subtypeData[,1]),3])
GoM3Type <- GoM3Type[!is.na(GoM3Type)]
cogType <- cogType[!is.na(cogType)]

samples <- colnames(exprs(expr.lumi))
labels <- as.factor(labelData[match(samples, labelData[,1]),2])
labeled <- !is.na(labels)
mod <- model.matrix(~as.factor(labels), data=labels)
arrays <- sapply(samples, function(x) {
  t <- unlist(strsplit(x, "_"))
  t[1]
})
arrays <- as.factor(as.vector(arrays))

geneIDs <- rownames(exprs(expr.lumi))
geneSelected <- grepl("hsa", geneIDs)

refGene <- list(c(52),c(9),c(52,9))
refGeneName <- list(c("U49"),c("U24"),c("U49","U24"))

expr.lumi.vst <- list()
expr.lumi.log <- list()

for(i in 1:length(refGene)) {
  expr.lumi.vst[[i]] <- refGeneNorm(expr.lumi, refGene[[i]], refGeneName[[i]], method="vst") 
  expr.lumi.log[[i]] <- refGeneNorm(expr.lumi, refGene[[i]], refGeneName[[i]], method="log2") 
}
expr.vst.U49.combat <- t(ComBat(dat=exprs(expr.lumi.vst[[1]])[ geneSelected, labeled], batch=arrays[labeled], mod=mod))
expr.vst.U24.combat <- t(ComBat(dat=exprs(expr.lumi.vst[[2]])[ geneSelected, labeled], batch=arrays[labeled], mod=mod))
expr.vst.U49U24.combat <- t(ComBat(dat=exprs(expr.lumi.vst[[3]])[ geneSelected, labeled], batch=arrays[labeled], mod=mod))
expr.log.U49.combat <- t(ComBat(dat=exprs(expr.lumi.log[[1]])[ geneSelected, labeled], batch=arrays[labeled], mod=mod))
expr.log.U24.combat <- t(ComBat(dat=exprs(expr.lumi.log[[2]])[ geneSelected, labeled], batch=arrays[labeled], mod=mod))
expr.log.U49U24.combat <- t(ComBat(dat=exprs(expr.lumi.log[[3]])[ geneSelected, labeled], batch=arrays[labeled], mod=mod))

# global normalisation methods 
# (3) lumiT: transformation 
expr.lumi.vst <- lumiT(expr.lumi) # default - VST

# (4) limuN: normaliation
expr.lumi.qn <- lumiN(expr.lumi.vst, method='quantile') # default: quantile
expr.lumi.rsn <- lumiN(expr.lumi.vst, method="rsn") # Robust Spline Normalization
expr.lumi.ssn <- lumiN(expr.lumi.vst, method="ssn") # Simple Scaling Normalization
expr.lumi.vsn <- lumiN(expr.lumi.vst, method="vsn") # Variance stabilization
expr.lumi.rin <- lumiN(expr.lumi.vst, method="rankinvariant") # Rank Invariant Normalization
expr.lumi.loess <- lumiN(expr.lumi.vst, method="loess") # Loess Normalization

expr.qn.combat <- t(ComBat(dat=exprs(expr.lumi.qn)[ geneSelected, labeled], batch=arrays[labeled], mod=mod))
expr.rsn.combat <- t(ComBat(dat=exprs(expr.lumi.rsn)[ geneSelected, labeled], batch=arrays[labeled], mod=mod))
expr.ssn.combat <- t(ComBat(dat=exprs(expr.lumi.ssn)[ geneSelected, labeled], batch=arrays[labeled], mod=mod))
expr.vsn.combat <- t(ComBat(dat=exprs(expr.lumi.vsn)[ geneSelected, labeled], batch=arrays[labeled], mod=mod))
expr.rin.combat <- t(ComBat(dat=exprs(expr.lumi.rin)[ geneSelected, labeled], batch=arrays[labeled], mod=mod))
expr.loess.combat <- t(ComBat(dat=exprs(expr.lumi.loess)[ geneSelected, labeled], batch=arrays[labeled], mod=mod))


cog1 <- which(cogType == 1)
cog2 <- which(cogType == 2)
GoM1 <- which(GoM3Type == 1)
GoM2 <- which(GoM3Type == 2)
GoM3 <- which(GoM3Type == 3)
ctrl <- which(labels[labeled] == "CTRL")


## cog 1v2
outprefix <- "cog1v2.expr"
sub_labels <- rep(NA_character_, sum(labeled))
sub_labels[cog1] <- "cog1"
sub_labels[cog2] <- "cog2"
sub_labeled <- !is.na(sub_labels)
  
write.table(data.frame(label=sub_labels[sub_labeled], expr.vst.U49.combat[sub_labeled,]), 
            file=paste(outprefix, "vst.U49.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.vst.U24.combat[sub_labeled,]), 
            file=paste(outprefix, "vst.U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.vst.U49U24.combat[sub_labeled,]), 
            file=paste(outprefix, "vst.U49U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.log.U49.combat[sub_labeled,]), 
            file=paste(outprefix, "log.U49.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.log.U24.combat[sub_labeled,]), 
            file=paste(outprefix, "log.U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.log.U49U24.combat[sub_labeled,]), 
            file=paste(outprefix, "log.U49U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)

write.table(data.frame(label=sub_labels[sub_labeled], expr.qn.combat[sub_labeled,]), 
            file=paste(outprefix, "qn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.rsn.combat[sub_labeled,]), 
            file=paste(outprefix, "rsn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.ssn.combat[sub_labeled,]), 
            file=paste(outprefix, "ssn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.vsn.combat[sub_labeled,]), 
            file=paste(outprefix, "vsn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.rin.combat[sub_labeled,]), 
            file=paste(outprefix, "rin.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.loess.combat[sub_labeled,]), 
            file=paste(outprefix, "loess.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)

### cog 1vc
outprefix <- "cog1vc.expr"
sub_labels <- rep(NA_character_, sum(labeled))
sub_labels[cog1] <- "cog1"
sub_labels[ctrl] <- "ctrl"
sub_labeled <- !is.na(sub_labels)

write.table(data.frame(label=sub_labels[sub_labeled], expr.vst.U49.combat[sub_labeled,]), 
            file=paste(outprefix, "vst.U49.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.vst.U24.combat[sub_labeled,]), 
            file=paste(outprefix, "vst.U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.vst.U49U24.combat[sub_labeled,]), 
            file=paste(outprefix, "vst.U49U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.log.U49.combat[sub_labeled,]), 
            file=paste(outprefix, "log.U49.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.log.U24.combat[sub_labeled,]), 
            file=paste(outprefix, "log.U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.log.U49U24.combat[sub_labeled,]), 
            file=paste(outprefix, "log.U49U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)

write.table(data.frame(label=sub_labels[sub_labeled], expr.qn.combat[sub_labeled,]), 
            file=paste(outprefix, "qn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.rsn.combat[sub_labeled,]), 
            file=paste(outprefix, "rsn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.ssn.combat[sub_labeled,]), 
            file=paste(outprefix, "ssn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.vsn.combat[sub_labeled,]), 
            file=paste(outprefix, "vsn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.rin.combat[sub_labeled,]), 
            file=paste(outprefix, "rin.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.loess.combat[sub_labeled,]), 
            file=paste(outprefix, "loess.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)


### cog 2vc
outprefix <- "cog2vc.expr"
sub_labels <- rep(NA_character_, sum(labeled))
sub_labels[cog2] <- "cog2"
sub_labels[ctrl] <- "ctrl"
sub_labeled <- !is.na(sub_labels)

write.table(data.frame(label=sub_labels[sub_labeled], expr.vst.U49.combat[sub_labeled,]), 
            file=paste(outprefix, "vst.U49.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.vst.U24.combat[sub_labeled,]), 
            file=paste(outprefix, "vst.U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.vst.U49U24.combat[sub_labeled,]), 
            file=paste(outprefix, "vst.U49U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.log.U49.combat[sub_labeled,]), 
            file=paste(outprefix, "log.U49.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.log.U24.combat[sub_labeled,]), 
            file=paste(outprefix, "log.U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.log.U49U24.combat[sub_labeled,]), 
            file=paste(outprefix, "log.U49U24.combat.txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)

write.table(data.frame(label=sub_labels[sub_labeled], expr.qn.combat[sub_labeled,]), 
            file=paste(outprefix, "qn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.rsn.combat[sub_labeled,]), 
            file=paste(outprefix, "rsn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.ssn.combat[sub_labeled,]), 
            file=paste(outprefix, "ssn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.vsn.combat[sub_labeled,]), 
            file=paste(outprefix, "vsn.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.rin.combat[sub_labeled,]), 
            file=paste(outprefix, "rin.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(data.frame(label=sub_labels[sub_labeled], expr.loess.combat[sub_labeled,]), 
            file=paste(outprefix, "loess.combat", "txt", sep="."), sep="\t", 
            quote=FALSE, row.names=TRUE, col.names=TRUE)

## END