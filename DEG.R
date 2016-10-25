# Expr data transformation and normalisation done with "normalisation_lumi.R" and "normalisation_refgene.R"
setwd("~/Projects/bioMarker/PBMCFromRawData")
source("normalisation_lumi.R")
source("normalisation_refgene.R")

################## Task here #######################
# DEG detection (w/o ComBat) for PBMC data 

# batch and mod data 
batch <- arrays
labeled <- !is.na(label) # some samples were not labeled 
mod <- model.matrix(~as.factor(label), data=label)
colnames(mod) <- c('All', 'SZ')
mod0 <- model.matrix(~1, data=label[labeled])
ngenes <- nrow(exprs(expr.lumi))
geneIDs <- rownames(expr.data)

expr.data.qn <- exprs(expr.lumi.qn)
expr.data.rsn <- exprs(expr.lumi.rsn)
expr.data.ssn <- exprs(expr.lumi.ssn)
expr.data.vsn <- exprs(expr.lumi.vsn)
expr.data.rin <- exprs(expr.lumi.rin)
expr.data.loess <- exprs(expr.lumi.loess)

# 1. UNDO combat
# 2. DEG detetction
# 2.1 using limma
require(limma)
expr.data.qn.fit <- lmFit(expr.data.qn[,labeled], mod)
expr.data.qn.fit <- eBayes(expr.data.qn.fit)
expr.data.qn.limma.res <- topTable(expr.data.qn.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.rsn.fit <- lmFit(expr.data.rsn[,labeled], mod)
expr.data.rsn.fit <- eBayes(expr.data.rsn.fit)
expr.data.rsn.limma.res <- topTable(expr.data.rsn.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.ssn.fit <- lmFit(expr.data.ssn[,labeled], mod)
expr.data.ssn.fit <- eBayes(expr.data.ssn.fit)
expr.data.ssn.limma.res <- topTable(expr.data.ssn.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.vsn.fit <- lmFit(expr.data.vsn[,labeled], mod)
expr.data.vsn.fit <- eBayes(expr.data.vsn.fit)
expr.data.vsn.limma.res <- topTable(expr.data.vsn.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.rin.fit <- lmFit(expr.data.rin[,labeled], mod)
expr.data.rin.fit <- eBayes(expr.data.rin.fit)
expr.data.rin.limma.res <- topTable(expr.data.rin.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.loess.fit <- lmFit(expr.data.loess[,labeled], mod)
expr.data.loess.fit <- eBayes(expr.data.loess.fit)
expr.data.loess.limma.res <- topTable(expr.data.loess.fit, coef='SZ', adjust='BH', number=ngenes)
# 
expr.data.dividedByU49.lumiT.fit <- lmFit(expr.data.dividedByU49.lumiT[,labeled], mod)
expr.data.dividedByU49.lumiT.fit <- eBayes(expr.data.dividedByU49.lumiT.fit)
expr.data.dividedByU49.lumiT.limma.res <- topTable(expr.data.dividedByU49.lumiT.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.VST.U49.fit <- lmFit(expr.data.VST.U49[,labeled], mod)
expr.data.VST.U49.fit <- eBayes(expr.data.VST.U49.fit)
expr.data.VST.U49.limma.res <- topTable(expr.data.VST.U49.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.dividedByU24U49.lumiT.fit <- lmFit(expr.data.dividedByU24U49.lumiT[,labeled], mod)
expr.data.dividedByU24U49.lumiT.fit <- eBayes(expr.data.dividedByU24U49.lumiT.fit)
expr.data.dividedByU24U49.lumiT.limma.res <- topTable(expr.data.dividedByU24U49.lumiT.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.VST.U24U49.fit <- lmFit(expr.data.VST.U24U49[,labeled], mod)
expr.data.VST.U24U49.fit <- eBayes(expr.data.VST.U24U49.fit)
expr.data.VST.U24U49.limma.res <- topTable(expr.data.VST.U24U49.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.dividedByU24.lumiT.fit <- lmFit(expr.data.dividedByU24.lumiT[,labeled], mod)
expr.data.dividedByU24.lumiT.fit <- eBayes(expr.data.dividedByU24.lumiT.fit)
expr.data.dividedByU24.lumiT.limma.res <- topTable(expr.data.dividedByU24.lumiT.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.VST.U24.fit <- lmFit(expr.data.VST.U24[,labeled], mod)
expr.data.VST.U24.fit <- eBayes(expr.data.VST.U24.fit)
expr.data.VST.U24.limma.res <- topTable(expr.data.VST.U24.fit, coef='SZ', adjust='BH', number=ngenes)

# End of 2.1 DEG using limma 

# 2.2 using f.pvalue in sva package 
require(sva)
expr.data.qn.f.padj <- p.adjust(f.pvalue(dat=expr.data.qn[,labeled], mod=mod, mod0=mod0))
expr.data.rsn.f.padj <- p.adjust(f.pvalue(dat=expr.data.rsn[,labeled], mod=mod, mod0=mod0))
expr.data.ssn.f.padj <- p.adjust(f.pvalue(dat=expr.data.ssn[,labeled], mod=mod, mod0=mod0))
expr.data.vsn.f.padj <- p.adjust(f.pvalue(dat=expr.data.vsn[,labeled], mod=mod, mod0=mod0))
expr.data.rin.f.padj <- p.adjust(f.pvalue(dat=expr.data.rin[,labeled], mod=mod, mod0=mod0))
expr.data.loess.f.padj <- p.adjust(f.pvalue(dat=expr.data.loess[,labeled], mod=mod, mod0=mod0))
#
expr.data.dividedByU49.lumiT.f.padj <- p.adjust(f.pvalue(dat=expr.data.dividedByU49.lumiT[,labeled], mod=mod, mod0=mod0))
expr.data.VST.U49.f.padj <- p.adjust(f.pvalue(dat=expr.data.VST.U49[,labeled], mod=mod, mod0=mod0))
expr.data.dividedByU24U49.lumiT.f.padj <- p.adjust(f.pvalue(dat=expr.data.dividedByU24U49.lumiT[,labeled], mod=mod, mod0=mod0))
expr.data.VST.U24U49.f.padj <- p.adjust(f.pvalue(dat=expr.data.VST.U24U49[,labeled], mod=mod, mod0=mod0))
expr.data.dividedByU24.lumiT.f.padj <- p.adjust(f.pvalue(dat=expr.data.dividedByU24.lumiT[,labeled], mod=mod, mod0=mod0))
expr.data.VST.U24.f.padj <- p.adjust(f.pvalue(dat=expr.data.VST.U24[,labeled], mod=mod, mod0=mod0))

# 2.3 using SAM 
require(samr)
delta <- 0.25
expr.data.qn.samr.data <- list(x=expr.data.qn[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.qn.samr.obj <- samr(expr.data.qn.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.qn.samr.delta.table <- samr.compute.delta.table(expr.data.qn.samr.obj, nvals=50)
expr.data.qn.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.qn.samr.obj, delta, expr.data.qn.samr.data, expr.data.qn.samr.delta.table, all.genes=TRUE)

expr.data.rsn.samr.data <- list(x=expr.data.rsn[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.rsn.samr.obj <- samr(expr.data.rsn.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.rsn.samr.delta.table <- samr.compute.delta.table(expr.data.rsn.samr.obj, nvals=50)
expr.data.rsn.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.rsn.samr.obj, delta, expr.data.rsn.samr.data, expr.data.rsn.samr.delta.table, all.genes=TRUE)

expr.data.ssn.samr.data <- list(x=expr.data.ssn[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.ssn.samr.obj <- samr(expr.data.ssn.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.ssn.samr.delta.table <- samr.compute.delta.table(expr.data.ssn.samr.obj, nvals=50)
expr.data.ssn.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.ssn.samr.obj, delta, expr.data.ssn.samr.data, expr.data.ssn.samr.delta.table, all.genes=TRUE)

expr.data.vsn.samr.data <- list(x=expr.data.vsn[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.vsn.samr.obj <- samr(expr.data.vsn.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.vsn.samr.delta.table <- samr.compute.delta.table(expr.data.vsn.samr.obj, nvals=50)
expr.data.vsn.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.vsn.samr.obj, delta, expr.data.vsn.samr.data, expr.data.vsn.samr.delta.table, all.genes=TRUE)

expr.data.rin.samr.data <- list(x=expr.data.rin[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.rin.samr.obj <- samr(expr.data.rin.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.rin.samr.delta.table <- samr.compute.delta.table(expr.data.rin.samr.obj, nvals=50)
expr.data.rin.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.rin.samr.obj, delta, expr.data.rin.samr.data, expr.data.rin.samr.delta.table, all.genes=TRUE)

expr.data.loess.samr.data <- list(x=expr.data.loess[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.loess.samr.obj <- samr(expr.data.loess.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.loess.samr.delta.table <- samr.compute.delta.table(expr.data.loess.samr.obj, nvals=50)
expr.data.loess.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.loess.samr.obj, delta, expr.data.loess.samr.data, expr.data.loess.samr.delta.table, all.genes=TRUE)

expr.data.dividedByU49.lumiT.samr.data <- list(x=expr.data.dividedByU49.lumiT[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU49.lumiT.samr.obj <- samr(expr.data.dividedByU49.lumiT.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.dividedByU49.lumiT.samr.delta.table <- samr.compute.delta.table(expr.data.dividedByU49.lumiT.samr.obj, nvals=50)
expr.data.dividedByU49.lumiT.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU49.lumiT.samr.obj, delta, expr.data.dividedByU49.lumiT.samr.data, expr.data.dividedByU49.lumiT.samr.delta.table, all.genes=TRUE)

expr.data.VST.U49.samr.data <- list(x=expr.data.VST.U49[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U49.samr.obj <- samr(expr.data.VST.U49.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.VST.U49.samr.delta.table <- samr.compute.delta.table(expr.data.VST.U49.samr.obj, nvals=50)
expr.data.VST.U49.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U49.samr.obj, delta, expr.data.VST.U49.samr.data, expr.data.VST.U49.samr.delta.table, all.genes=TRUE)

expr.data.dividedByU24U49.lumiT.samr.data <- list(x=expr.data.dividedByU24U49.lumiT[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU24U49.lumiT.samr.obj <- samr(expr.data.dividedByU24U49.lumiT.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.dividedByU24U49.lumiT.samr.delta.table <- samr.compute.delta.table(expr.data.dividedByU24U49.lumiT.samr.obj, nvals=50)
expr.data.dividedByU24U49.lumiT.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU24U49.lumiT.samr.obj, delta, expr.data.dividedByU24U49.lumiT.samr.data, expr.data.dividedByU24U49.lumiT.samr.delta.table, all.genes=TRUE)

expr.data.VST.U24U49.samr.data <- list(x=expr.data.VST.U24U49[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U24U49.samr.obj <- samr(expr.data.VST.U24U49.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.VST.U24U49.samr.delta.table <- samr.compute.delta.table(expr.data.VST.U24U49.samr.obj, nvals=50)
expr.data.VST.U24U49.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U24U49.samr.obj, delta, expr.data.VST.U24U49.samr.data, expr.data.VST.U24U49.samr.delta.table, all.genes=TRUE)

expr.data.dividedByU24.lumiT.samr.data <- list(x=expr.data.dividedByU24.lumiT[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU24.lumiT.samr.obj <- samr(expr.data.dividedByU24.lumiT.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.dividedByU24.lumiT.samr.delta.table <- samr.compute.delta.table(expr.data.dividedByU24.lumiT.samr.obj, nvals=50)
expr.data.dividedByU24.lumiT.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU24.lumiT.samr.obj, delta, expr.data.dividedByU24.lumiT.samr.data, expr.data.dividedByU24.lumiT.samr.delta.table, all.genes=TRUE)

expr.data.VST.U24.samr.data <- list(x=expr.data.VST.U24[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U24.samr.obj <- samr(expr.data.VST.U24.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.VST.U24.samr.delta.table <- samr.compute.delta.table(expr.data.VST.U24.samr.obj, nvals=50)
expr.data.VST.U24.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U24.samr.obj, delta, expr.data.VST.U24.samr.data, expr.data.VST.U24.samr.delta.table, all.genes=TRUE)

# End of 2.3 DEG using samr 

# 2.3 using SAM w/ wilcoxon stat 
require(samr)
delta <- 0.25
expr.data.qn.samr.wxc.data <- list(x=expr.data.qn[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.qn.samr.wxc.obj <- samr(expr.data.qn.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.qn.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.qn.samr.wxc.obj, nvals=50)
expr.data.qn.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.qn.samr.wxc.obj, delta, expr.data.qn.samr.wxc.data, expr.data.qn.samr.wxc.delta.table, all.genes=TRUE)

expr.data.rsn.samr.wxc.data <- list(x=expr.data.rsn[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.rsn.samr.wxc.obj <- samr(expr.data.rsn.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.rsn.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.rsn.samr.wxc.obj, nvals=50)
expr.data.rsn.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.rsn.samr.wxc.obj, delta, expr.data.rsn.samr.wxc.data, expr.data.rsn.samr.wxc.delta.table, all.genes=TRUE)

expr.data.ssn.samr.wxc.data <- list(x=expr.data.ssn[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.ssn.samr.wxc.obj <- samr(expr.data.ssn.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.ssn.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.ssn.samr.wxc.obj, nvals=50)
expr.data.ssn.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.ssn.samr.wxc.obj, delta, expr.data.ssn.samr.wxc.data, expr.data.ssn.samr.wxc.delta.table, all.genes=TRUE)

expr.data.vsn.samr.wxc.data <- list(x=expr.data.vsn[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.vsn.samr.wxc.obj <- samr(expr.data.vsn.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.vsn.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.vsn.samr.wxc.obj, nvals=50)
expr.data.vsn.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.vsn.samr.wxc.obj, delta, expr.data.vsn.samr.wxc.data, expr.data.vsn.samr.wxc.delta.table, all.genes=TRUE)

expr.data.rin.samr.wxc.data <- list(x=expr.data.rin[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.rin.samr.wxc.obj <- samr(expr.data.rin.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.rin.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.rin.samr.wxc.obj, nvals=50)
expr.data.rin.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.rin.samr.wxc.obj, delta, expr.data.rin.samr.wxc.data, expr.data.rin.samr.wxc.delta.table, all.genes=TRUE)

expr.data.loess.samr.wxc.data <- list(x=expr.data.loess[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.loess.samr.wxc.obj <- samr(expr.data.loess.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.loess.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.loess.samr.wxc.obj, nvals=50)
expr.data.loess.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.loess.samr.wxc.obj, delta, expr.data.loess.samr.wxc.data, expr.data.loess.samr.wxc.delta.table, all.genes=TRUE)

expr.data.dividedByU49.lumiT.samr.wxc.data <- list(x=expr.data.dividedByU49.lumiT[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU49.lumiT.samr.wxc.obj <- samr(expr.data.dividedByU49.lumiT.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.dividedByU49.lumiT.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.dividedByU49.lumiT.samr.wxc.obj, nvals=50)
expr.data.dividedByU49.lumiT.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU49.lumiT.samr.wxc.obj, delta, expr.data.dividedByU49.lumiT.samr.wxc.data, expr.data.dividedByU49.lumiT.samr.wxc.delta.table, all.genes=TRUE)

expr.data.VST.U49.samr.wxc.data <- list(x=expr.data.VST.U49[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U49.samr.wxc.obj <- samr(expr.data.VST.U49.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.VST.U49.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.VST.U49.samr.wxc.obj, nvals=50)
expr.data.VST.U49.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U49.samr.wxc.obj, delta, expr.data.VST.U49.samr.wxc.data, expr.data.VST.U49.samr.wxc.delta.table, all.genes=TRUE)

expr.data.dividedByU24U49.lumiT.samr.wxc.data <- list(x=expr.data.dividedByU24U49.lumiT[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU24U49.lumiT.samr.wxc.obj <- samr(expr.data.dividedByU24U49.lumiT.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.dividedByU24U49.lumiT.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.dividedByU24U49.lumiT.samr.wxc.obj, nvals=50)
expr.data.dividedByU24U49.lumiT.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU24U49.lumiT.samr.wxc.obj, delta, expr.data.dividedByU24U49.lumiT.samr.wxc.data, expr.data.dividedByU24U49.lumiT.samr.wxc.delta.table, all.genes=TRUE)

expr.data.VST.U24U49.samr.wxc.data <- list(x=expr.data.VST.U24U49[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U24U49.samr.wxc.obj <- samr(expr.data.VST.U24U49.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.VST.U24U49.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.VST.U24U49.samr.wxc.obj, nvals=50)
expr.data.VST.U24U49.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U24U49.samr.wxc.obj, delta, expr.data.VST.U24U49.samr.wxc.data, expr.data.VST.U24U49.samr.wxc.delta.table, all.genes=TRUE)

expr.data.dividedByU24.lumiT.samr.wxc.data <- list(x=expr.data.dividedByU24.lumiT[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU24.lumiT.samr.wxc.obj <- samr(expr.data.dividedByU24.lumiT.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.dividedByU24.lumiT.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.dividedByU24.lumiT.samr.wxc.obj, nvals=50)
expr.data.dividedByU24.lumiT.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU24.lumiT.samr.wxc.obj, delta, expr.data.dividedByU24.lumiT.samr.wxc.data, expr.data.dividedByU24.lumiT.samr.wxc.delta.table, all.genes=TRUE)

expr.data.VST.U24.samr.wxc.data <- list(x=expr.data.VST.U24[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U24.samr.wxc.obj <- samr(expr.data.VST.U24.samr.wxc.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.VST.U24.samr.wxc.delta.table <- samr.compute.delta.table(expr.data.VST.U24.samr.wxc.obj, nvals=50)
expr.data.VST.U24.samr.wxc.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U24.samr.wxc.obj, delta, expr.data.VST.U24.samr.wxc.data, expr.data.VST.U24.samr.wxc.delta.table, all.genes=TRUE)

# End of 2.4 DEG using samr w/ wilcoxon
#================================

# write results to txt files 
write.table(expr.data.qn.limma.res, "PBMC.211samp.qn.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.rsn.limma.res, "PBMC.211samp.rsn.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.ssn.limma.res, "PBMC.211samp.ssn.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.vsn.limma.res, "PBMC.211samp.vsn.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.rin.limma.res, "PBMC.211samp.rin.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.loess.limma.res, "PBMC.211samp.loess.limma.res.txt", quote=FALSE, sep='\t')
#
write.table(expr.data.dividedByU49.lumiT.limma.res, "PBMC.211samp.U49_log.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U49.limma.res, "PBMC.211samp.VST_U49.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.dividedByU24U49.lumiT.limma.res, "PBMC.211samp.U24U49_log.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U24U49.limma.res, "PBMC.211samp.VST_U24U49.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.dividedByU24.lumiT.limma.res, "PBMC.211samp.U24_log.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U24.limma.res, "PBMC.211samp.VST_U24.limma.res.txt", quote=FALSE, sep='\t')
###
write.table(expr.data.qn.f.padj, "PBMC.211samp.qn.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.rsn.f.padj, "PBMC.211samp.rsn.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.ssn.f.padj, "PBMC.211samp.ssn.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.vsn.f.padj, "PBMC.211samp.vsn.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.rin.f.padj, "PBMC.211samp.rin.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.loess.f.padj, "PBMC.211samp.loess.f.padj.txt", quote=FALSE, sep='\t')
#
write.table(expr.data.dividedByU49.lumiT.f.padj, "PBMC.211samp.U49_log.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U49.f.padj, "PBMC.211samp.VST_U49.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.dividedByU24U49.lumiT.f.padj, "PBMC.211samp.U24U49_log.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U24U49.f.padj, "PBMC.211samp.VST_U24U49.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.dividedByU24.lumiT.f.padj, "PBMC.211samp.U24_log.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U24.f.padj, "PBMC.211samp.VST_U24.f.padj.txt", quote=FALSE, sep='\t')
#####
samr.qval.table <- data.frame(qn=rep(NA_real_, ngenes), 
                              rsn=rep(NA_real_, ngenes), 
                              ssn=rep(NA_real_, ngenes), 
                              vsn=rep(NA_real_, ngenes), 
                              rin=rep(NA_real_, ngenes), 
                              loess=rep(NA_real_, ngenes), 
                              U49_log=rep(NA_real_, ngenes), 
                              VST_U49=rep(NA_real_, ngenes), 
                              U24U49_log=rep(NA_real_, ngenes), 
                              VST_U24U49=rep(NA_real_, ngenes), 
                              U24_log=rep(NA_real_, ngenes), 
                              VST_U24=rep(NA_real_, ngenes),
                              row.names=geneIDs)
samr.qval.table$qn[match(expr.data.qn.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.qn.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$qn[match(expr.data.qn.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.qn.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$rsn[match(expr.data.rsn.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.rsn.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$rsn[match(expr.data.rsn.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.rsn.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$ssn[match(expr.data.ssn.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.ssn.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$ssn[match(expr.data.ssn.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.ssn.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$vsn[match(expr.data.vsn.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.vsn.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$vsn[match(expr.data.vsn.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.vsn.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$rin[match(expr.data.rin.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.rin.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$rin[match(expr.data.rin.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.rin.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$loess[match(expr.data.loess.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.loess.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$loess[match(expr.data.loess.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.loess.samr.siggenes.table$genes.lo[,8]) / 100.0
#
samr.qval.table$U49_log[match(expr.data.dividedByU49.lumiT.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU49.lumiT.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$U49_log[match(expr.data.dividedByU49.lumiT.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU49.lumiT.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$VST_U49[match(expr.data.VST.U49.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U49.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$VST_U49[match(expr.data.VST.U49.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U49.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$U24U49_log[match(expr.data.dividedByU24U49.lumiT.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24U49.lumiT.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$U24U49_log[match(expr.data.dividedByU24U49.lumiT.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24U49.lumiT.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$VST_U24U49[match(expr.data.VST.U24U49.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U24U49.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$VST_U24U49[match(expr.data.VST.U24U49.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U24U49.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$U24_log[match(expr.data.dividedByU24.lumiT.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24.lumiT.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$U24_log[match(expr.data.dividedByU24.lumiT.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24.lumiT.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$VST_U24[match(expr.data.VST.U24.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U24.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$VST_U24[match(expr.data.VST.U24.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U24.samr.siggenes.table$genes.lo[,8]) / 100.0
#
write.table(samr.qval.table, "PBMC.211samp.nonComBat.samr.qval.txt", quote=FALSE, sep='\t')

### plot SAM objects
# pdf("SAM.ttest.woComBat.pdf")
samr.plot(expr.data.qn.samr.obj, delta); title(main="VST+QN")
samr.plot(expr.data.rsn.samr.obj, delta); title(main="VST+RSN")
samr.plot(expr.data.ssn.samr.obj, delta); title(main="VST+SSN")
samr.plot(expr.data.vsn.samr.obj, delta); title(main="VST+VSN")
samr.plot(expr.data.rin.samr.obj, delta); title(main="VST+RIN")
samr.plot(expr.data.loess.samr.obj, delta); title(main="VST+LOESS")
samr.plot(expr.data.dividedByU49.lumiT.samr.obj, delta); title(main="U49+log")
samr.plot(expr.data.VST.U49.samr.obj, delta); title(main="VST+U49")
samr.plot(expr.data.dividedByU24U49.lumiT.samr.obj, delta); title(main="U24U49+log")
samr.plot(expr.data.VST.U24U49.samr.obj, delta); title(main="VST+U24U49")
samr.plot(expr.data.dividedByU24.lumiT.samr.obj, delta); title(main="U24+log")
samr.plot(expr.data.VST.U24.samr.obj, delta); title(main="VST+U24")
# dev.off()
####
# pdf("SAM.wcx.woComBat.pdf")
samr.plot(expr.data.qn.samr.wxc.obj, delta); title(main="VST+QN")
samr.plot(expr.data.rsn.samr.wxc.obj, delta); title(main="VST+RSN")
samr.plot(expr.data.ssn.samr.wxc.obj, delta); title(main="VST+SSN")
samr.plot(expr.data.vsn.samr.wxc.obj, delta); title(main="VST+VSN")
samr.plot(expr.data.rin.samr.wxc.obj, delta); title(main="VST+RIN")
samr.plot(expr.data.loess.samr.wxc.obj, delta); title(main="VST+LOESS")
samr.plot(expr.data.dividedByU49.lumiT.samr.wxc.obj, delta); title(main="U49+log")
samr.plot(expr.data.VST.U49.samr.wxc.obj, delta); title(main="VST+U49")
samr.plot(expr.data.dividedByU24U49.lumiT.samr.wxc.obj, delta); title(main="U24U49+log")
samr.plot(expr.data.VST.U24U49.samr.wxc.obj, delta); title(main="VST+U24U49")
samr.plot(expr.data.dividedByU24.lumiT.samr.wxc.obj, delta); title(main="U24+log")
samr.plot(expr.data.VST.U24.samr.wxc.obj, delta); title(main="VST+U24")
# dev.off()

### plot limma fit
volcanoplot2 <- function(fit, coef=2, adjP=0.05, ...) {
  xlab = "Log Fold Change"
  ylab = "Log Odds"
  pch = 16
  cex = 0.5
  x <- as.matrix(fit$coef)[, coef]
  y <- as.matrix(fit$lods)[, coef]
  plot(x, y, xlab = xlab, ylab = ylab, pch = pch+as.numeric(p.adjust(fit$p.value[, coef], method="BH")<adjP), 
       cex = cex, col=1+as.numeric(p.adjust(fit$p.value[, coef], method="BH")<adjP), ...)
  #invisible()
}
# pdf("limma.volcano.woCombat.pdf")
volcanoplot2(expr.data.qn.fit, main="VST+QN")
volcanoplot2(expr.data.rsn.fit, main="VST+RSN")
volcanoplot2(expr.data.ssn.fit, main="VST+SSN")
volcanoplot2(expr.data.vsn.fit, main="VST+VSN")
volcanoplot2(expr.data.rin.fit, main="VST+RIN")
volcanoplot2(expr.data.loess.fit, main="VST+LOESS")
volcanoplot2(expr.data.dividedByU49.lumiT.fit, main="U49+log")
volcanoplot2(expr.data.VST.U49.fit, main="VST+U49")
volcanoplot2(expr.data.dividedByU24U49.lumiT.fit, main="U24U49+log")
volcanoplot2(expr.data.VST.U24U49.fit, main="VST+U24U49")
volcanoplot2(expr.data.dividedByU24.lumiT.fit, main="U24+log")
volcanoplot2(expr.data.VST.U24.fit, main="VST+U24")
# dev.off()
# END 
