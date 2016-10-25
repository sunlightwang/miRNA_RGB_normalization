# Expr data transformation and normalisation done with "normalisation_lumi.R" and "normalisation_refgene.R"
setwd("~/Projects/bioMarker/PBMCFromRawData")
source("normalisation_lumi.R")
source("normalisation_refgene.R")

################## Task here #######################
# Batch effect adjust and DEG detection for PBMC data 

# 1. Batch effect adjustment using ComBat in sva bioconductor package
require(sva)

# batch and mod data 
batch <- arrays
labeled <- !is.na(label) # some samples were not labeled 
mod <- model.matrix(~as.factor(label), data=label)
colnames(mod) <- c('All', 'SZ')
mod0 <- model.matrix(~1, data=label[labeled])
ngenes <- nrow(exprs(expr.lumi))
geneIDs <- rownames(expr.data)

# for diff normalisation methods from lumi package
expr.data.qn <- exprs(expr.lumi.qn)
expr.data.qn.combat <- ComBat(dat=expr.data.qn[,labeled], batch=batch[labeled], mod=mod)
expr.data.rsn <- exprs(expr.lumi.rsn)
expr.data.rsn.combat <- ComBat(dat=expr.data.rsn[,labeled], batch=batch[labeled], mod=mod)
expr.data.ssn <- exprs(expr.lumi.ssn)
expr.data.ssn.combat <- ComBat(dat=expr.data.ssn[,labeled], batch=batch[labeled], mod=mod)
expr.data.vsn <- exprs(expr.lumi.vsn)
expr.data.vsn.combat <- ComBat(dat=expr.data.vsn[,labeled], batch=batch[labeled], mod=mod)
expr.data.rin <- exprs(expr.lumi.rin)
expr.data.rin.combat <- ComBat(dat=expr.data.rin[,labeled], batch=batch[labeled], mod=mod)
expr.data.loess <- exprs(expr.lumi.loess)
expr.data.loess.combat <- ComBat(dat=expr.data.loess[,labeled], batch=batch[labeled], mod=mod)
#
expr.data.dividedByU49.lumiT <- exprs(expr.lumi.dividedByU49.lumiT)
expr.data.dividedByU49.lumiT.combat <- ComBat(dat=expr.data.dividedByU49.lumiT[,labeled], batch=batch[labeled], mod=mod)
expr.data.VST.U49 <- exprs(expr.lumi.VST.U49)
expr.data.VST.U49.combat <- ComBat(dat=expr.data.VST.U49[,labeled], batch=batch[labeled], mod=mod)
expr.data.dividedByU24U49.lumiT <- exprs(expr.lumi.dividedByU24U49.lumiT)
expr.data.dividedByU24U49.lumiT.combat <- ComBat(dat=expr.data.dividedByU24U49.lumiT[,labeled], batch=batch[labeled], mod=mod)
expr.data.VST.U24U49 <- exprs(expr.lumi.VST.U24U49)
expr.data.VST.U24U49.combat <- ComBat(dat=expr.data.VST.U24U49[,labeled], batch=batch[labeled], mod=mod)
expr.data.dividedByU24.lumiT <- exprs(expr.lumi.dividedByU24.lumiT)
expr.data.dividedByU24.lumiT.combat <- ComBat(dat=expr.data.dividedByU24.lumiT[,labeled], batch=batch[labeled], mod=mod)
expr.data.VST.U24 <- exprs(expr.lumi.VST.U24)
expr.data.VST.U24.combat <- ComBat(dat=expr.data.VST.U24[,labeled], batch=batch[labeled], mod=mod)
# End of 1. ComBat

# 2. DEG detetction
# 2.1 using limma
require(limma)
expr.data.qn.combat.fit <- lmFit(expr.data.qn.combat, mod)
expr.data.qn.combat.fit <- eBayes(expr.data.qn.combat.fit)
expr.data.qn.combat.limma.res <- topTable(expr.data.qn.combat.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.rsn.combat.fit <- lmFit(expr.data.rsn.combat, mod)
expr.data.rsn.combat.fit <- eBayes(expr.data.rsn.combat.fit)
expr.data.rsn.combat.limma.res <- topTable(expr.data.rsn.combat.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.ssn.combat.fit <- lmFit(expr.data.ssn.combat, mod)
expr.data.ssn.combat.fit <- eBayes(expr.data.ssn.combat.fit)
expr.data.ssn.combat.limma.res <- topTable(expr.data.ssn.combat.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.vsn.combat.fit <- lmFit(expr.data.vsn.combat, mod)
expr.data.vsn.combat.fit <- eBayes(expr.data.vsn.combat.fit)
expr.data.vsn.combat.limma.res <- topTable(expr.data.vsn.combat.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.rin.combat.fit <- lmFit(expr.data.rin.combat, mod)
expr.data.rin.combat.fit <- eBayes(expr.data.rin.combat.fit)
expr.data.rin.combat.limma.res <- topTable(expr.data.rin.combat.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.loess.combat.fit <- lmFit(expr.data.loess.combat, mod)
expr.data.loess.combat.fit <- eBayes(expr.data.loess.combat.fit)
expr.data.loess.combat.limma.res <- topTable(expr.data.loess.combat.fit, coef='SZ', adjust='BH', number=ngenes)
# 
expr.data.dividedByU49.lumiT.combat.fit <- lmFit(expr.data.dividedByU49.lumiT.combat, mod)
expr.data.dividedByU49.lumiT.combat.fit <- eBayes(expr.data.dividedByU49.lumiT.combat.fit)
expr.data.dividedByU49.lumiT.combat.limma.res <- topTable(expr.data.dividedByU49.lumiT.combat.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.VST.U49.combat.fit <- lmFit(expr.data.VST.U49.combat, mod)
expr.data.VST.U49.combat.fit <- eBayes(expr.data.VST.U49.combat.fit)
expr.data.VST.U49.combat.limma.res <- topTable(expr.data.VST.U49.combat.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.dividedByU24U49.lumiT.combat.fit <- lmFit(expr.data.dividedByU24U49.lumiT.combat, mod)
expr.data.dividedByU24U49.lumiT.combat.fit <- eBayes(expr.data.dividedByU24U49.lumiT.combat.fit)
expr.data.dividedByU24U49.lumiT.combat.limma.res <- topTable(expr.data.dividedByU24U49.lumiT.combat.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.VST.U24U49.combat.fit <- lmFit(expr.data.VST.U24U49.combat, mod)
expr.data.VST.U24U49.combat.fit <- eBayes(expr.data.VST.U24U49.combat.fit)
expr.data.VST.U24U49.combat.limma.res <- topTable(expr.data.VST.U24U49.combat.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.dividedByU24.lumiT.combat.fit <- lmFit(expr.data.dividedByU24.lumiT.combat, mod)
expr.data.dividedByU24.lumiT.combat.fit <- eBayes(expr.data.dividedByU24.lumiT.combat.fit)
expr.data.dividedByU24.lumiT.combat.limma.res <- topTable(expr.data.dividedByU24.lumiT.combat.fit, coef='SZ', adjust='BH', number=ngenes)
expr.data.VST.U24.combat.fit <- lmFit(expr.data.VST.U24.combat, mod)
expr.data.VST.U24.combat.fit <- eBayes(expr.data.VST.U24.combat.fit)
expr.data.VST.U24.combat.limma.res <- topTable(expr.data.VST.U24.combat.fit, coef='SZ', adjust='BH', number=ngenes)

# End of 2.1 DEG using limma 

# 2.2 using f.pvalue in sva package 
require(sva)
expr.data.qn.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.qn.combat, mod=mod, mod0=mod0))
expr.data.rsn.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.rsn.combat, mod=mod, mod0=mod0))
expr.data.ssn.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.ssn.combat, mod=mod, mod0=mod0))
expr.data.vsn.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.vsn.combat, mod=mod, mod0=mod0))
expr.data.rin.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.rin.combat, mod=mod, mod0=mod0))
expr.data.loess.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.loess.combat, mod=mod, mod0=mod0))
#
expr.data.dividedByU49.lumiT.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.dividedByU49.lumiT.combat, mod=mod, mod0=mod0))
expr.data.VST.U49.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.VST.U49.combat, mod=mod, mod0=mod0))
expr.data.dividedByU24U49.lumiT.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.dividedByU24U49.lumiT.combat, mod=mod, mod0=mod0))
expr.data.VST.U24U49.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.VST.U24U49.combat, mod=mod, mod0=mod0))
expr.data.dividedByU24.lumiT.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.dividedByU24.lumiT.combat, mod=mod, mod0=mod0))
expr.data.VST.U24.combat.f.padj <- p.adjust(f.pvalue(dat=expr.data.VST.U24.combat, mod=mod, mod0=mod0))

# 2.3 using SAM 
require(samr)
delta <- 0.25
expr.data.qn.combat.samr.data <- list(x=expr.data.qn.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.qn.combat.samr.obj <- samr(expr.data.qn.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.qn.combat.samr.delta.table <- samr.compute.delta.table(expr.data.qn.combat.samr.obj, nvals=50)
expr.data.qn.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.qn.combat.samr.obj, delta, expr.data.qn.combat.samr.data, expr.data.qn.combat.samr.delta.table, all.genes=TRUE)

expr.data.rsn.combat.samr.data <- list(x=expr.data.rsn.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.rsn.combat.samr.obj <- samr(expr.data.rsn.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.rsn.combat.samr.delta.table <- samr.compute.delta.table(expr.data.rsn.combat.samr.obj, nvals=50)
expr.data.rsn.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.rsn.combat.samr.obj, delta, expr.data.rsn.combat.samr.data, expr.data.rsn.combat.samr.delta.table, all.genes=TRUE)

expr.data.ssn.combat.samr.data <- list(x=expr.data.ssn.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.ssn.combat.samr.obj <- samr(expr.data.ssn.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.ssn.combat.samr.delta.table <- samr.compute.delta.table(expr.data.ssn.combat.samr.obj, nvals=50)
expr.data.ssn.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.ssn.combat.samr.obj, delta, expr.data.ssn.combat.samr.data, expr.data.ssn.combat.samr.delta.table, all.genes=TRUE)

expr.data.vsn.combat.samr.data <- list(x=expr.data.vsn.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.vsn.combat.samr.obj <- samr(expr.data.vsn.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.vsn.combat.samr.delta.table <- samr.compute.delta.table(expr.data.vsn.combat.samr.obj, nvals=50)
expr.data.vsn.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.vsn.combat.samr.obj, delta, expr.data.vsn.combat.samr.data, expr.data.vsn.combat.samr.delta.table, all.genes=TRUE)

expr.data.rin.combat.samr.data <- list(x=expr.data.rin.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.rin.combat.samr.obj <- samr(expr.data.rin.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.rin.combat.samr.delta.table <- samr.compute.delta.table(expr.data.rin.combat.samr.obj, nvals=50)
expr.data.rin.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.rin.combat.samr.obj, delta, expr.data.rin.combat.samr.data, expr.data.rin.combat.samr.delta.table, all.genes=TRUE)

expr.data.loess.combat.samr.data <- list(x=expr.data.loess.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.loess.combat.samr.obj <- samr(expr.data.loess.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.loess.combat.samr.delta.table <- samr.compute.delta.table(expr.data.loess.combat.samr.obj, nvals=50)
expr.data.loess.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.loess.combat.samr.obj, delta, expr.data.loess.combat.samr.data, expr.data.loess.combat.samr.delta.table, all.genes=TRUE)

expr.data.dividedByU49.lumiT.combat.samr.data <- list(x=expr.data.dividedByU49.lumiT.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU49.lumiT.combat.samr.obj <- samr(expr.data.dividedByU49.lumiT.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.dividedByU49.lumiT.combat.samr.delta.table <- samr.compute.delta.table(expr.data.dividedByU49.lumiT.combat.samr.obj, nvals=50)
expr.data.dividedByU49.lumiT.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU49.lumiT.combat.samr.obj, delta, expr.data.dividedByU49.lumiT.combat.samr.data, expr.data.dividedByU49.lumiT.combat.samr.delta.table, all.genes=TRUE)

expr.data.VST.U49.combat.samr.data <- list(x=expr.data.VST.U49.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U49.combat.samr.obj <- samr(expr.data.VST.U49.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.VST.U49.combat.samr.delta.table <- samr.compute.delta.table(expr.data.VST.U49.combat.samr.obj, nvals=50)
expr.data.VST.U49.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U49.combat.samr.obj, delta, expr.data.VST.U49.combat.samr.data, expr.data.VST.U49.combat.samr.delta.table, all.genes=TRUE)

expr.data.dividedByU24U49.lumiT.combat.samr.data <- list(x=expr.data.dividedByU24U49.lumiT.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU24U49.lumiT.combat.samr.obj <- samr(expr.data.dividedByU24U49.lumiT.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.dividedByU24U49.lumiT.combat.samr.delta.table <- samr.compute.delta.table(expr.data.dividedByU24U49.lumiT.combat.samr.obj, nvals=50)
expr.data.dividedByU24U49.lumiT.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU24U49.lumiT.combat.samr.obj, delta, expr.data.dividedByU24U49.lumiT.combat.samr.data, expr.data.dividedByU24U49.lumiT.combat.samr.delta.table, all.genes=TRUE)

expr.data.VST.U24U49.combat.samr.data <- list(x=expr.data.VST.U24U49.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U24U49.combat.samr.obj <- samr(expr.data.VST.U24U49.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.VST.U24U49.combat.samr.delta.table <- samr.compute.delta.table(expr.data.VST.U24U49.combat.samr.obj, nvals=50)
expr.data.VST.U24U49.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U24U49.combat.samr.obj, delta, expr.data.VST.U24U49.combat.samr.data, expr.data.VST.U24U49.combat.samr.delta.table, all.genes=TRUE)

expr.data.dividedByU24.lumiT.combat.samr.data <- list(x=expr.data.dividedByU24.lumiT.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU24.lumiT.combat.samr.obj <- samr(expr.data.dividedByU24.lumiT.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.dividedByU24.lumiT.combat.samr.delta.table <- samr.compute.delta.table(expr.data.dividedByU24.lumiT.combat.samr.obj, nvals=50)
expr.data.dividedByU24.lumiT.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU24.lumiT.combat.samr.obj, delta, expr.data.dividedByU24.lumiT.combat.samr.data, expr.data.dividedByU24.lumiT.combat.samr.delta.table, all.genes=TRUE)

expr.data.VST.U24.combat.samr.data <- list(x=expr.data.VST.U24.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U24.combat.samr.obj <- samr(expr.data.VST.U24.combat.samr.data, resp.type="Two class unpaired", nperms=1000)
expr.data.VST.U24.combat.samr.delta.table <- samr.compute.delta.table(expr.data.VST.U24.combat.samr.obj, nvals=50)
expr.data.VST.U24.combat.samr.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U24.combat.samr.obj, delta, expr.data.VST.U24.combat.samr.data, expr.data.VST.U24.combat.samr.delta.table, all.genes=TRUE)

# End of 2.3 DEG using samr 

# 2.4 using SAM w/ wilcoxon statistic
require(samr)
delta <- 0.25
expr.data.qn.combat.samr.wcx.data <- list(x=expr.data.qn.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.qn.combat.samr.wcx.obj <- samr(expr.data.qn.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.qn.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.qn.combat.samr.wcx.obj, nvals=50)
expr.data.qn.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.qn.combat.samr.wcx.obj, delta, expr.data.qn.combat.samr.wcx.data, expr.data.qn.combat.samr.wcx.delta.table, all.genes=TRUE)

expr.data.rsn.combat.samr.wcx.data <- list(x=expr.data.rsn.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.rsn.combat.samr.wcx.obj <- samr(expr.data.rsn.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.rsn.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.rsn.combat.samr.wcx.obj, nvals=50)
expr.data.rsn.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.rsn.combat.samr.wcx.obj, delta, expr.data.rsn.combat.samr.wcx.data, expr.data.rsn.combat.samr.wcx.delta.table, all.genes=TRUE)

expr.data.ssn.combat.samr.wcx.data <- list(x=expr.data.ssn.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.ssn.combat.samr.wcx.obj <- samr(expr.data.ssn.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.ssn.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.ssn.combat.samr.wcx.obj, nvals=50)
expr.data.ssn.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.ssn.combat.samr.wcx.obj, delta, expr.data.ssn.combat.samr.wcx.data, expr.data.ssn.combat.samr.wcx.delta.table, all.genes=TRUE)

expr.data.vsn.combat.samr.wcx.data <- list(x=expr.data.vsn.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.vsn.combat.samr.wcx.obj <- samr(expr.data.vsn.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.vsn.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.vsn.combat.samr.wcx.obj, nvals=50)
expr.data.vsn.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.vsn.combat.samr.wcx.obj, delta, expr.data.vsn.combat.samr.wcx.data, expr.data.vsn.combat.samr.wcx.delta.table, all.genes=TRUE)

expr.data.rin.combat.samr.wcx.data <- list(x=expr.data.rin.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.rin.combat.samr.wcx.obj <- samr(expr.data.rin.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.rin.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.rin.combat.samr.wcx.obj, nvals=50)
expr.data.rin.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.rin.combat.samr.wcx.obj, delta, expr.data.rin.combat.samr.wcx.data, expr.data.rin.combat.samr.wcx.delta.table, all.genes=TRUE)

expr.data.loess.combat.samr.wcx.data <- list(x=expr.data.loess.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.loess.combat.samr.wcx.obj <- samr(expr.data.loess.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.loess.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.loess.combat.samr.wcx.obj, nvals=50)
expr.data.loess.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.loess.combat.samr.wcx.obj, delta, expr.data.loess.combat.samr.wcx.data, expr.data.loess.combat.samr.wcx.delta.table, all.genes=TRUE)

expr.data.dividedByU49.lumiT.combat.samr.wcx.data <- list(x=expr.data.dividedByU49.lumiT.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU49.lumiT.combat.samr.wcx.obj <- samr(expr.data.dividedByU49.lumiT.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.dividedByU49.lumiT.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.dividedByU49.lumiT.combat.samr.wcx.obj, nvals=50)
expr.data.dividedByU49.lumiT.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU49.lumiT.combat.samr.wcx.obj, delta, expr.data.dividedByU49.lumiT.combat.samr.wcx.data, expr.data.dividedByU49.lumiT.combat.samr.wcx.delta.table, all.genes=TRUE)

expr.data.VST.U49.combat.samr.wcx.data <- list(x=expr.data.VST.U49.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U49.combat.samr.wcx.obj <- samr(expr.data.VST.U49.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.VST.U49.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.VST.U49.combat.samr.wcx.obj, nvals=50)
expr.data.VST.U49.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U49.combat.samr.wcx.obj, delta, expr.data.VST.U49.combat.samr.wcx.data, expr.data.VST.U49.combat.samr.wcx.delta.table, all.genes=TRUE)

expr.data.dividedByU24U49.lumiT.combat.samr.wcx.data <- list(x=expr.data.dividedByU24U49.lumiT.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU24U49.lumiT.combat.samr.wcx.obj <- samr(expr.data.dividedByU24U49.lumiT.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.dividedByU24U49.lumiT.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.dividedByU24U49.lumiT.combat.samr.wcx.obj, nvals=50)
expr.data.dividedByU24U49.lumiT.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU24U49.lumiT.combat.samr.wcx.obj, delta, expr.data.dividedByU24U49.lumiT.combat.samr.wcx.data, expr.data.dividedByU24U49.lumiT.combat.samr.wcx.delta.table, all.genes=TRUE)

expr.data.VST.U24U49.combat.samr.wcx.data <- list(x=expr.data.VST.U24U49.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U24U49.combat.samr.wcx.obj <- samr(expr.data.VST.U24U49.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.VST.U24U49.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.VST.U24U49.combat.samr.wcx.obj, nvals=50)
expr.data.VST.U24U49.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U24U49.combat.samr.wcx.obj, delta, expr.data.VST.U24U49.combat.samr.wcx.data, expr.data.VST.U24U49.combat.samr.wcx.delta.table, all.genes=TRUE)

expr.data.dividedByU24.lumiT.combat.samr.wcx.data <- list(x=expr.data.dividedByU24.lumiT.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU24.lumiT.combat.samr.wcx.obj <- samr(expr.data.dividedByU24.lumiT.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.dividedByU24.lumiT.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.dividedByU24.lumiT.combat.samr.wcx.obj, nvals=50)
expr.data.dividedByU24.lumiT.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.dividedByU24.lumiT.combat.samr.wcx.obj, delta, expr.data.dividedByU24.lumiT.combat.samr.wcx.data, expr.data.dividedByU24.lumiT.combat.samr.wcx.delta.table, all.genes=TRUE)

expr.data.VST.U24.combat.samr.wcx.data <- list(x=expr.data.VST.U24.combat, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.VST.U24.combat.samr.wcx.obj <- samr(expr.data.VST.U24.combat.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=1000)
expr.data.VST.U24.combat.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.VST.U24.combat.samr.wcx.obj, nvals=50)
expr.data.VST.U24.combat.samr.wcx.siggenes.table<-samr.compute.siggenes.table(expr.data.VST.U24.combat.samr.wcx.obj, delta, expr.data.VST.U24.combat.samr.wcx.data, expr.data.VST.U24.combat.samr.wcx.delta.table, all.genes=TRUE)

# End of 2.4 DEG using samr with wilcoxon

# write results to txt files 
write.table(expr.data.qn.combat.limma.res, "PBMC.211samp.qn.combat.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.rsn.combat.limma.res, "PBMC.211samp.rsn.combat.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.ssn.combat.limma.res, "PBMC.211samp.ssn.combat.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.vsn.combat.limma.res, "PBMC.211samp.vsn.combat.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.rin.combat.limma.res, "PBMC.211samp.rin.combat.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.loess.combat.limma.res, "PBMC.211samp.loess.combat.limma.res.txt", quote=FALSE, sep='\t')
#
write.table(expr.data.dividedByU49.lumiT.combat.limma.res, "PBMC.211samp.U49_log.combat.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U49.combat.limma.res, "PBMC.211samp.VST_U49.combat.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.dividedByU24U49.lumiT.combat.limma.res, "PBMC.211samp.U24U49_log.combat.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U24U49.combat.limma.res, "PBMC.211samp.VST_U24U49.combat.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.dividedByU24.lumiT.combat.limma.res, "PBMC.211samp.U24_log.combat.limma.res.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U24.combat.limma.res, "PBMC.211samp.VST_U24.combat.limma.res.txt", quote=FALSE, sep='\t')
###
write.table(expr.data.qn.combat.f.padj, "PBMC.211samp.qn.combat.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.rsn.combat.f.padj, "PBMC.211samp.rsn.combat.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.ssn.combat.f.padj, "PBMC.211samp.ssn.combat.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.vsn.combat.f.padj, "PBMC.211samp.vsn.combat.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.rin.combat.f.padj, "PBMC.211samp.rin.combat.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.loess.combat.f.padj, "PBMC.211samp.loess.combat.f.padj.txt", quote=FALSE, sep='\t')
#
write.table(expr.data.dividedByU49.lumiT.combat.f.padj, "PBMC.211samp.U49_log.combat.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U49.combat.f.padj, "PBMC.211samp.VST_U49.combat.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.dividedByU24U49.lumiT.combat.f.padj, "PBMC.211samp.U24U49_log.combat.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U24U49.combat.f.padj, "PBMC.211samp.VST_U24U49.combat.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.dividedByU24.lumiT.combat.f.padj, "PBMC.211samp.U24_log.combat.f.padj.txt", quote=FALSE, sep='\t')
write.table(expr.data.VST.U24.combat.f.padj, "PBMC.211samp.VST_U24.combat.f.padj.txt", quote=FALSE, sep='\t')
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
samr.qval.table$qn[match(expr.data.qn.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.qn.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$qn[match(expr.data.qn.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.qn.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$rsn[match(expr.data.rsn.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.rsn.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$rsn[match(expr.data.rsn.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.rsn.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$ssn[match(expr.data.ssn.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.ssn.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$ssn[match(expr.data.ssn.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.ssn.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$vsn[match(expr.data.vsn.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.vsn.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$vsn[match(expr.data.vsn.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.vsn.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$rin[match(expr.data.rin.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.rin.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$rin[match(expr.data.rin.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.rin.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$loess[match(expr.data.loess.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.loess.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$loess[match(expr.data.loess.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.loess.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
#
samr.qval.table$U49_log[match(expr.data.dividedByU49.lumiT.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU49.lumiT.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$U49_log[match(expr.data.dividedByU49.lumiT.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU49.lumiT.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$VST_U49[match(expr.data.VST.U49.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U49.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$VST_U49[match(expr.data.VST.U49.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U49.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$U24U49_log[match(expr.data.dividedByU24U49.lumiT.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24U49.lumiT.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$U24U49_log[match(expr.data.dividedByU24U49.lumiT.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24U49.lumiT.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$VST_U24U49[match(expr.data.VST.U24U49.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U24U49.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$VST_U24U49[match(expr.data.VST.U24U49.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U24U49.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$U24_log[match(expr.data.dividedByU24.lumiT.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24.lumiT.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$U24_log[match(expr.data.dividedByU24.lumiT.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24.lumiT.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
samr.qval.table$VST_U24[match(expr.data.VST.U24.combat.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U24.combat.samr.siggenes.table$genes.up[,8]) / 100.0
samr.qval.table$VST_U24[match(expr.data.VST.U24.combat.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.VST.U24.combat.samr.siggenes.table$genes.lo[,8]) / 100.0
#
write.table(samr.qval.table, "PBMC.211samp.samr.qval.txt", quote=FALSE, sep='\t')

####
### plot SAM objects
# pdf("SAM.ttest.wComBat.pdf")
samr.plot(expr.data.qn.combat.samr.obj, delta); title(main="VST+QN+ComBat")
samr.plot(expr.data.rsn.combat.samr.obj, delta); title(main="VST+RSN+ComBat")
samr.plot(expr.data.ssn.combat.samr.obj, delta); title(main="VST+SSN+ComBat")
samr.plot(expr.data.vsn.combat.samr.obj, delta); title(main="VST+VSN+ComBat")
samr.plot(expr.data.rin.combat.samr.obj, delta); title(main="VST+RIN+ComBat")
samr.plot(expr.data.loess.combat.samr.obj, delta); title(main="VST+LOESS+ComBat")
samr.plot(expr.data.dividedByU49.lumiT.combat.samr.obj, delta); title(main="U49+log+ComBat")
samr.plot(expr.data.VST.U49.combat.samr.obj, delta); title(main="VST+U49+ComBat")
samr.plot(expr.data.dividedByU24U49.lumiT.combat.samr.obj, delta); title(main="U24U49+log+ComBat")
samr.plot(expr.data.VST.U24U49.combat.samr.obj, delta); title(main="VST+U24U49+ComBat")
samr.plot(expr.data.dividedByU24.lumiT.combat.samr.obj, delta); title(main="U24+log+ComBat")
samr.plot(expr.data.VST.U24.combat.samr.obj, delta); title(main="VST+U24+ComBat")
# dev.off()
####
# pdf("SAM.wcx.wComBat.pdf")
samr.plot(expr.data.qn.combat.samr.wcx.obj, delta); title(main="VST+QN+ComBat")
samr.plot(expr.data.rsn.combat.samr.wcx.obj, delta); title(main="VST+RSN+ComBat")
samr.plot(expr.data.ssn.combat.samr.wcx.obj, delta); title(main="VST+SSN+ComBat")
samr.plot(expr.data.vsn.combat.samr.wcx.obj, delta); title(main="VST+VSN+ComBat")
samr.plot(expr.data.rin.combat.samr.wcx.obj, delta); title(main="VST+RIN+ComBat")
samr.plot(expr.data.loess.combat.samr.wcx.obj, delta); title(main="VST+LOESS+ComBat")
samr.plot(expr.data.dividedByU49.lumiT.combat.samr.wcx.obj, delta); title(main="U49+log+ComBat")
samr.plot(expr.data.VST.U49.combat.samr.wcx.obj, delta); title(main="VST+U49+ComBat")
samr.plot(expr.data.dividedByU24U49.lumiT.combat.samr.wcx.obj, delta); title(main="U24U49+log+ComBat")
samr.plot(expr.data.VST.U24U49.combat.samr.wcx.obj, delta); title(main="VST+U24U49+ComBat")
samr.plot(expr.data.dividedByU24.lumiT.combat.samr.wcx.obj, delta); title(main="U24+log+ComBat")
samr.plot(expr.data.VST.U24.combat.samr.wcx.obj, delta); title(main="VST+U24+ComBat")
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
# pdf("limma.volcano.wCombat.pdf")
volcanoplot2(expr.data.qn.combat.fit, main="VST+QN")
volcanoplot2(expr.data.rsn.combat.fit, main="VST+RSN")
volcanoplot2(expr.data.ssn.combat.fit, main="VST+SSN")
volcanoplot2(expr.data.vsn.combat.fit, main="VST+VSN")
volcanoplot2(expr.data.rin.combat.fit, main="VST+RIN")
volcanoplot2(expr.data.loess.combat.fit, main="VST+LOESS")
volcanoplot2(expr.data.dividedByU49.lumiT.combat.fit, main="U49+log")
volcanoplot2(expr.data.VST.U49.combat.fit, main="VST+U49")
volcanoplot2(expr.data.dividedByU24U49.lumiT.combat.fit, main="U24U49+log")
volcanoplot2(expr.data.VST.U24U49.combat.fit, main="VST+U24U49")
volcanoplot2(expr.data.dividedByU24.lumiT.combat.fit, main="U24+log")
volcanoplot2(expr.data.VST.U24.combat.fit, main="VST+U24")
# dev.off()

# END 

