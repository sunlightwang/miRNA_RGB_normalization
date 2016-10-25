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

limmaDEG <- function(x.lumi, mod, labeled=rep(TRUE, ncol(expr.lumi)), resultFile, comBat=FALSE, batch) {
  require(limma)
  require(sva)
  stopifnot(is(x.lumi, "eSet"))
  
  expr.data <- exprs(x.lumi)
  if(comBat) {
    expr.data <- ComBat(dat=expr.data[,labeled], batch=batch[labeled], mod=mod)
    expr.data.fit <- lmFit(expr.data, mod)
  } else {
    expr.data.fit <- lmFit(expr.data[,labeled], mod)
  }
  expr.data.fit <- eBayes(expr.data.fit)
  ngenes <- nrow(expr.data)
  expr.data.fit.res <- topTable(expr.data.fit, coef=colnames(mod)[2], adjust='BH', number=ngenes)
  write.table(expr.data.fit.res, resultFile, quote=FALSE, sep='\t', row.names=FALSE)
}

ftestDEG <- function(x.lumi, mod, mod0, labeled=rep(TRUE, ncol(expr.lumi)), resultFile, comBat=FALSE, batch) {
  require(sva)
  stopifnot(is(x.lumi, "eSet"))
  
  expr.data <- exprs(x.lumi)
  if(comBat) {
    expr.data <- ComBat(dat=expr.data[,labeled], batch=batch[labeled], mod=mod)
    pval <- f.pvalue(dat=expr.data, mod=mod, mod0=mod0)
  } else 
  {
    pval <- f.pvalue(dat=expr.data[,labeled], mod=mod, mod0=mod0)
  }
  res <- data.frame(cbind(geneID=rownames(expr.data), pval))
  res <- res[order(res[,2]),]
  write.table(res, resultFile, quote=FALSE, sep='\t', row.names=FALSE)
}


samrDEG <- function(x.lumi, mod, labeled=rep(TRUE, ncol(expr.lumi)), resultFile, comBat=FALSE, batch) {
  require(sva)
  require(samr)
  stopifnot(is(x.lumi, "eSet"))
  
  expr.data <- exprs(x.lumi)
  geneIDs <- rownames(expr.data)
  if(comBat) {
    expr.data <- ComBat(dat=expr.data[,labeled], batch=batch[labeled], mod=mod)
    expr.data.samr.data <- list(x=expr.data, y=as.factor(mod[,2]), geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
  } else {
    expr.data.samr.data <- list(x=expr.data[,labeled], y=as.factor(mod[,2]), geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
  }
  delta <- 0.25
  expr.data.samr.obj <- samr(expr.data.samr.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
  expr.data.samr.delta.table <- samr.compute.delta.table(expr.data.samr.obj, nvals=50)
  expr.data.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.samr.obj, delta, expr.data.samr.data, expr.data.samr.delta.table, all.genes=TRUE)
  res <- rbind(expr.data.samr.siggenes.table$genes.up, expr.data.samr.siggenes.table$genes.lo)
  res <- res[ order(as.double(res[,8]), -abs(log(as.double(res[,7])))), ]
  write.table(res[,-c(1,2)], resultFile, quote=FALSE, sep='\t', row.names=FALSE)
}

filename <- function(outFilePrefix, norm, comBat, DE) { 
  if(comBat) {
    return( paste(outFilePrefix, norm, "combat", DE, "rst", sep=".") )
  } else {
    return( paste(outFilePrefix, norm, DE, "rst", sep=".") ) 
  }
}

arrayAnalysis <- function(dataFile, sampleInfoFile, controlDataFile, labels, refGene=list(c(52),c(9),c(52,9)), 
                          refGeneName=list(c("U49"),c("U24"),c("U49","U24")), outFilePrefix="ArrayAnalysisResult", comBat=FALSE, ...) {
# dataFile: expression raw file
# sampleInfoFile: sample info file
# controlDataFile: control data
# label: two columns, fisrt col matchs sample IDs in expression file, second col contains the label
# refGene: the row number of the ref gene expression values in controlData
  require(lumi)
  source("lumi_edited.R")
  stopifnot(length(refGene) == length(refGeneName))

  ########### lumi normalization ##########################################
  # (0) lumiB: data loading
  expr.lumi <- lumiR.batch(dataFile, sampleInfoFile = sampleInfoFile)
  expr.lumi <- addControlData2lumi(controlDataFile, expr.lumi)
  expr.data <- exprs(expr.lumi)
  
  # (1) generate plots
  # plot regarding control data 
  plotControlData(expr.lumi, main="Control data boxplot: original")
  boxplot(expr.lumi, main="Sample expression boxplot: original data")
  plotDensity(expr.lumi, main="Sample expression density: original data")
  
  # (2) lumiB: backgroup substraction 
  # expr.lumi <- lumiB(expr.lumi) # none will be done, as the bg has already been substracted 
  
  # (3) lumiT: transformation 
  expr.lumi.vst <- lumiT(expr.lumi) # default - VST
  boxplot(expr.lumi.vst, main="Sample expression boxplot: after VS transformation")
  plotDensity(expr.lumi.vst, main="Sample expression density: after VST")
  
  # (4) limuN: normaliation
  expr.lumi.qn <- lumiN(expr.lumi.vst, method='quantile') # default: quantile
  boxplot(expr.lumi.qn, main="Sample expression boxplot: after quantile normalisation")
  plotDensity(expr.lumi.qn, main="Sample expression density: after quantile normalisation")
  
  expr.lumi.rsn <- lumiN(expr.lumi.vst, method="rsn") # Robust Spline Normalization
  boxplot(expr.lumi.rsn, main="Sample expression boxplot: after Robust Spline Normalization")
  plotDensity(expr.lumi.rsn, main="Sample expression density: after Robust Spline Normalization")
  
  expr.lumi.ssn <- lumiN(expr.lumi.vst, method="ssn") # Simple Scaling Normalization
  boxplot(expr.lumi.ssn, main="Sample expression boxplot: after Simple Scaling Normalization")
  plotDensity(expr.lumi.ssn, xlim=c(2,6), main="Sample expression density: after Simple Scaling Normalization")
  
  expr.lumi.vsn <- lumiN(expr.lumi.vst, method="vsn") # Variance stabilization
  boxplot(expr.lumi.vsn, main="Sample expression boxplot: after Variance Stabilization Normalisation")
  plotDensity(expr.lumi.vsn, main="Sample expression density: after Variance Stabilization Normalisation")
  
  expr.lumi.rin <- lumiN(expr.lumi.vst, method="rankinvariant") # Rank Invariant Normalization
  boxplot(expr.lumi.rin, main="Sample expression boxplot: after Rank Invariant Normalization")
  plotDensity(expr.lumi.rin, xlim=c(2,6), main="Sample expression density: after Rank Invariant Normalization")
  
  expr.lumi.loess <- lumiN(expr.lumi.vst, method="loess") # Loess Normalization
  boxplot(expr.lumi.loess, main="Sample expression boxplot: after Loess Normalization")
  plotDensity(expr.lumi.loess, main="Sample expression density: after Loess Normalization")
  
  # (4) limuQ: QC
  #expr.lumi <- lumiQ(expr.lumi)
  # and the results by different normalisation methods 
  
  ########### END lumi normalization ######################################
  
  ########### refGene normalization #######################################
  expr.refGene.vst.norm <- list()
  expr.refGene.log.norm <- list()
  for(i in 1:length(refGene)) {
    expr.refGene.vst.norm[[i]] <- refGeneNorm(expr.lumi, refGene[[i]], refGeneName[[i]], method="vst") 
    expr.refGene.log.norm[[i]] <- refGeneNorm(expr.lumi, refGene[[i]], refGeneName[[i]], method="log2") 
  }
  ########### END refGene normalization ###################################
  
  samples <- colnames(expr.data)
  arrays <- sapply(samples, function(x) {
    t <- unlist(strsplit(x, "_"))
    t[1]
  })
  arrays <- as.factor(as.vector(arrays))
  # get label (schiz vs. control) information 
  labels <- as.factor(labels[match(samples, labels[,1]),2])
  labeled <- !is.na(labels)
  mod <- model.matrix(~as.factor(labels), data=labels)
  colnames(mod) <- c('All', 'SZ')
  mod0 <- model.matrix(~1, data=labels[labeled])
  ngenes <- nrow(expr.data)
  geneIDs <- rownames(expr.data)
  
  ########### DEG detection with limma ####################################
  limmaDEG(expr.lumi.qn, mod=mod, labeled=labeled,  resultFile=filename(outFilePrefix, "vst+qn", comBat, "limma"), comBat=comBat, batch=arrays)
  limmaDEG(expr.lumi.rsn, mod=mod, labeled=labeled,  resultFile=filename(outFilePrefix, "vst+rsn", comBat, "limma"), comBat=comBat, batch=arrays)
  limmaDEG(expr.lumi.ssn, mod=mod, labeled=labeled,  resultFile=filename(outFilePrefix, "vst+ssn", comBat, "limma"), comBat=comBat, batch=arrays)
  limmaDEG(expr.lumi.vsn, mod=mod, labeled=labeled,  resultFile=filename(outFilePrefix, "vst+vsn", comBat, "limma"), comBat=comBat, batch=arrays)
  limmaDEG(expr.lumi.rin, mod=mod, labeled=labeled,  resultFile=filename(outFilePrefix, "vst+rin", comBat, "limma"), comBat=comBat, batch=arrays)
  limmaDEG(expr.lumi.loess, mod=mod, labeled=labeled,  resultFile=filename(outFilePrefix, "vst+loess", comBat, "limma"), comBat=comBat, batch=arrays)
  for(i in 1:length(refGene)) {
    limmaDEG(expr.refGene.vst.norm[[i]], mod=mod, labeled=labeled,  resultFile=filename(outFilePrefix, paste("vst", paste(refGeneName[[i]], collapse=""), sep="."), comBat, "limma"), comBat=comBat, batch=arrays)
    limmaDEG(expr.refGene.log.norm[[i]], mod=mod, labeled=labeled,  resultFile=filename(outFilePrefix, paste("log", paste(refGeneName[[i]], collapse=""), sep="."), comBat, "limma"), comBat=comBat, batch=arrays)
  }
  ########### END DEG detection with limma ################################
  
  ########### DEG detection with ftest ####################################
  ftestDEG(expr.lumi.qn, mod=mod, mod0=mod0, labeled=labeled, resultFile=filename(outFilePrefix, "vst+qn", comBat, "ftest"), comBat=comBat, batch=arrays) 
  ftestDEG(expr.lumi.rsn, mod=mod, mod0=mod0, labeled=labeled, resultFile=filename(outFilePrefix, "vst+rsn", comBat, "ftest"), comBat=comBat, batch=arrays) 
  ftestDEG(expr.lumi.ssn, mod=mod, mod0=mod0, labeled=labeled, resultFile=filename(outFilePrefix, "vst+ssn", comBat, "ftest"), comBat=comBat, batch=arrays) 
  ftestDEG(expr.lumi.vsn, mod=mod, mod0=mod0, labeled=labeled, resultFile=filename(outFilePrefix, "vst+vsn", comBat, "ftest"), comBat=comBat, batch=arrays) 
  ftestDEG(expr.lumi.rin, mod=mod, mod0=mod0, labeled=labeled, resultFile=filename(outFilePrefix, "vst+rin", comBat, "ftest"), comBat=comBat, batch=arrays) 
  ftestDEG(expr.lumi.loess, mod=mod, mod0=mod0, labeled=labeled, resultFile=filename(outFilePrefix, "vst+loess", comBat, "ftest"), comBat=comBat, batch=arrays) 
  for(i in 1:length(refGene)) {
    ftestDEG(expr.refGene.vst.norm[[i]], mod=mod, mod0=mod0,  labeled=labeled,  resultFile=filename(outFilePrefix, paste("vst", paste(refGeneName[[i]], collapse=""), sep="."), comBat, "ftest"), comBat=comBat, batch=arrays)
    ftestDEG(expr.refGene.log.norm[[i]], mod=mod, mod0=mod0,  labeled=labeled,  resultFile=filename(outFilePrefix, paste("log", paste(refGeneName[[i]], collapse=""), sep="."), comBat, "ftest"), comBat=comBat, batch=arrays)
  }
  ########### END DEG detection with ftest ################################
    
  ########### DEG detection with samr #####################################
  samrDEG(expr.lumi.qn, mod=mod, labeled=labeled, resultFile=filename(outFilePrefix, "vst+qn", comBat, "samr"), comBat=comBat, batch=arrays)
  samrDEG(expr.lumi.rsn, mod=mod, labeled=labeled, resultFile=filename(outFilePrefix, "vst+rsn", comBat, "samr"), comBat=comBat, batch=arrays)
  samrDEG(expr.lumi.ssn, mod=mod, labeled=labeled, resultFile=filename(outFilePrefix, "vst+ssn", comBat, "samr"), comBat=comBat, batch=arrays)
  samrDEG(expr.lumi.vsn, mod=mod, labeled=labeled, resultFile=filename(outFilePrefix, "vst+vsn", comBat, "samr"), comBat=comBat, batch=arrays)
  samrDEG(expr.lumi.rin, mod=mod, labeled=labeled, resultFile=filename(outFilePrefix, "vst+rin", comBat, "samr"), comBat=comBat, batch=arrays)
  samrDEG(expr.lumi.loess, mod=mod, labeled=labeled, resultFile=filename(outFilePrefix, "vst+loess", comBat, "samr"), comBat=comBat, batch=arrays)
  for(i in 1:length(refGene)) {
    samrDEG(expr.refGene.vst.norm[[i]], mod=mod, labeled=labeled,  resultFile=filename(outFilePrefix, paste("vst", paste(refGeneName[[i]], collapse=""), sep="."), comBat, "samr"), comBat=comBat, batch=arrays)
    samrDEG(expr.refGene.log.norm[[i]], mod=mod, labeled=labeled,  resultFile=filename(outFilePrefix, paste("log", paste(refGeneName[[i]], collapse=""), sep="."), comBat, "samr"), comBat=comBat, batch=arrays)
  }
  ########### DEG detection with samr #####################################
}

################################################
################################################

DEresultAnalysis <- function(dir, keys, norm.method, DEG.method="limma", top=100) {
  require(Vennerable)
  res <- list()
  n <- top
  tops <- list()
  for(i in 1:length(keys)) {
    res[[i]] <- read.table(paste(paste(dir, keys[i], sep=""), norm.method, DEG.method, "rst", sep="."), header=T, sep="\t")
    tops[[i]] <- res[[i]][1:n,1]
  }
  v <- Venn(tops)
  x <- v@IndicatorWeight
  VennScoring(x)
#   v <- venn(list(A019=res[[1]][1:n,1], 
#                    A032=res[[2]][1:n,1], 
#                    A035=res[[3]][1:n,1], 
#                    A081=res[[4]][1:n,1]
#   ))
#   title(main=paste(norm.method, DEG.method, "with top", top, sep=" "))
  # score <- 1 * overlap by two + 2 * overlap by 3 + 4 * overlap by 4
#   v["0011",1] + v["0101",1] + v["0110",1] + v["1001",1] + v["1010",1] + v["1100",1] + 
#     2 * (v["0111",1] + v["1011",1] + v["1101",1] + v["1110",1]) +
#     4 * v["1111", 1]
}

VennScoring <- function(Venn.IndicatorWeight) {
  x <- Venn.IndicatorWeight 
  stopifnot(2^(ncol(x)-1) == nrow(x))
  y <- rowSums(x[,-ncol(x)])
  sum(x[y>1,ncol(x)] * 2 ^ (y[y>1] - 2))
}

####################################################
### DE analysis of array data: main
####################################################
### Set parameters here
####################################################

args <- commandArgs(TRUE)
if(length(args) < 2 || length(args) > 2)
{
  cat ("USAGE: Rscript arrayPipe.R <dir:string> <combat:logical>\n")
  quit("no")
}

dir <- paste(args[1], "/", sep="")
combat <- as.logical(args[2])
if( is.na(combat)) {
  cat ("combat:logical must be LOGICAL\n")
  cat ("USAGE: Rscript arrayPipe.R <dir:string> <combat:logical>\n")
  quit("no")
}
keys <- c("A001","A002","A003","A004")
labelFile <- "sampleLabel211.txt"
labels <- read.table(labelFile, header=T)[,c(1,4)]
methods <- c("vst+qn", "vst+rsn", "vst+ssn", "vst+vsn","vst+rin", "vst+loess", "log.U49", "log.U24", "log.U49U24", "vst.U49", "vst.U24", "vst.U49U24")
out.prefix <- paste(dir, "conscore.", sep="")
if(combat) {
  methods <- paste(methods, "combat", sep=".")
  out.prefix <- paste(out.prefix, "combat.", sep="")
}
top <- c(5, 10,20,50,100,150,200,250)



for(key in keys) {
  dataFile <- paste(dir, key, "_expr.txt", sep="")
  sampleInfoFile <- paste(dir, key, "_sample.txt", sep="") 
  controlDataFile <-  paste(dir, key, "_ctrl.txt", sep="")
  outFilePrefix <-  paste(dir, key, sep="")
  arrayAnalysis(dataFile=dataFile, sampleInfoFile=sampleInfoFile, controlDataFile=controlDataFile, labels=labels, outFilePrefix=outFilePrefix, comBat=combat)
}
# analysis of the results 
scores.limma <- matrix(0, nrow=length(methods), ncol=length(top))
rownames(scores.limma) <- methods
colnames(scores.limma) <- top
for(i in 1:length(methods)) {
  for(j in 1:length(top)) {
    scores.limma[i,j] <- DEresultAnalysis(dir, keys, methods[i], top=top[j])
  }
}
scores.ftest <- matrix(0, nrow=length(methods), ncol=length(top))
rownames(scores.ftest) <- methods
colnames(scores.ftest) <- top
for(i in 1:length(methods)) {
  for(j in 1:length(top)) {
    scores.ftest[i,j] <- DEresultAnalysis(dir, keys, methods[i], DEG.method="ftest", top=top[j])
  }
}
scores.samr <- matrix(0, nrow=length(methods), ncol=length(top))
rownames(scores.samr) <- methods
colnames(scores.samr) <- top
for(i in 1:length(methods)) {
  for(j in 1:length(top)) {
    scores.samr[i,j] <- DEresultAnalysis(dir, keys, methods[i], DEG.method="samr", top=top[j])
  }
}

write.table(scores.limma, paste(out.prefix, "limma", sep=""), col.names=T, row.names=T, quote=F, sep="\t")
write.table(scores.ftest, paste(out.prefix, "ftest", sep=""), col.names=T, row.names=T, quote=F, sep="\t")
write.table(scores.samr, paste(out.prefix, "samr", sep=""), col.names=T, row.names=T, quote=F, sep="\t")

# end 
