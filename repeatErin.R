# This is only to repeat Erin's analysis 

# used miRNA expression values 
# expr.data.dividedByU24U49
# or 
# expr.data.dividedByU49

delta <- 0.25

# ref gene normalised unlogged values 
expr.data.dividedByU24U49.samr.data <- list(x=expr.data.dividedByU24U49[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=FALSE)
expr.data.dividedByU24U49.samr.obj <- samr(expr.data.dividedByU24U49.samr.data, resp.type="Two class unpaired", nperms=5000)
expr.data.dividedByU24U49.samr.delta.table <- samr.compute.delta.table(expr.data.dividedByU24U49.samr.obj, nvals=50)
expr.data.dividedByU24U49.samr.siggenes.table <- samr.compute.siggenes.table(expr.data.dividedByU24U49.samr.obj, delta, expr.data.dividedByU24U49.samr.data, expr.data.dividedByU24U49.samr.delta.table)

expr.data.dividedByU49.samr.data <- list(x=expr.data.dividedByU49[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=FALSE)
expr.data.dividedByU49.samr.obj <- samr(expr.data.dividedByU49.samr.data, resp.type="Two class unpaired", nperms=5000)
expr.data.dividedByU49.samr.delta.table <- samr.compute.delta.table(expr.data.dividedByU49.samr.obj, nvals=50)
expr.data.dividedByU49.samr.siggenes.table <- samr.compute.siggenes.table(expr.data.dividedByU49.samr.obj, delta, expr.data.dividedByU49.samr.data, expr.data.dividedByU49.samr.delta.table)

# ref gene normalised logged un-combatted values 
expr.data.dividedByU24U49.lumiT.samr.data <- list(x=expr.data.dividedByU24U49.lumiT[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU24U49.lumiT.samr.obj <- samr(expr.data.dividedByU24U49.lumiT.samr.data, resp.type="Two class unpaired", nperms=5000)
expr.data.dividedByU24U49.lumiT.samr.delta.table <- samr.compute.delta.table(expr.data.dividedByU24U49.lumiT.samr.obj, nvals=50)
expr.data.dividedByU24U49.lumiT.samr.siggenes.table <- samr.compute.siggenes.table(expr.data.dividedByU24U49.lumiT.samr.obj, delta, expr.data.dividedByU24U49.lumiT.samr.data, expr.data.dividedByU24U49.lumiT.samr.delta.table)

expr.data.dividedByU49.lumiT.samr.data <- list(x=expr.data.dividedByU49.lumiT[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU49.lumiT.samr.obj <- samr(expr.data.dividedByU49.lumiT.samr.data, resp.type="Two class unpaired", nperms=5000)
expr.data.dividedByU49.lumiT.samr.delta.table <- samr.compute.delta.table(expr.data.dividedByU49.lumiT.samr.obj, nvals=50)
expr.data.dividedByU49.lumiT.samr.siggenes.table <- samr.compute.siggenes.table(expr.data.dividedByU49.lumiT.samr.obj, delta, expr.data.dividedByU49.lumiT.samr.data, expr.data.dividedByU49.lumiT.samr.delta.table)


# exp VST+QN+COMBAT value
expr.data.qn.combat.exp <- exp(expr.data.qn.combat)
expr.data.qn.combat.exp.samr.data <- list(x=expr.data.qn.combat.exp, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=FALSE)
expr.data.qn.combat.exp.samr.obj <- samr(expr.data.qn.combat.exp.samr.data, resp.type="Two class unpaired", nperms=5000)
expr.data.qn.combat.exp.samr.delta.table <- samr.compute.delta.table(expr.data.qn.combat.exp.samr.obj, nvals=50)
expr.data.qn.combat.exp.samr.siggenes.table <- samr.compute.siggenes.table(expr.data.qn.combat.exp.samr.obj, delta, expr.data.qn.combat.exp.samr.data, expr.data.qn.combat.exp.samr.delta.table)

# exp VST+QN value
expr.data.qn.exp <- exp(expr.data.qn[,labeled])
expr.data.qn.exp.samr.data <- list(x=expr.data.qn.exp, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=FALSE)
expr.data.qn.exp.samr.obj <- samr(expr.data.qn.exp.samr.data, resp.type="Two class unpaired", nperms=5000)
expr.data.qn.exp.samr.delta.table <- samr.compute.delta.table(expr.data.qn.exp.samr.obj, nvals=50)
expr.data.qn.exp.samr.siggenes.table <- samr.compute.siggenes.table(expr.data.qn.exp.samr.obj, delta, expr.data.qn.exp.samr.data, expr.data.qn.exp.samr.delta.table)

# VST+QN value
expr.data.qn.samr.data <- list(x=expr.data.qn[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=FALSE)
expr.data.qn.samr.obj <- samr(expr.data.qn.samr.data, resp.type="Two class unpaired", nperms=5000)
expr.data.qn.samr.delta.table <- samr.compute.delta.table(expr.data.qn.samr.obj, nvals=50)
expr.data.qn.samr.siggenes.table <- samr.compute.siggenes.table(expr.data.qn.samr.obj, delta, expr.data.qn.samr.data, expr.data.qn.samr.delta.table)

###
# output 
trial.samr.qval.table <- data.frame(U24U49=rep(NA_real_, ngenes), 
                                    U49=rep(NA_real_, ngenes), 
                                    U24U49.log=rep(NA_real_, ngenes), 
                                    U49.log=rep(NA_real_, ngenes), 
                                    VST.QN.ComBat.exp=rep(NA_real_, ngenes), 
                                    VST.QN.exp=rep(NA_real_, ngenes), 
                                    VST.QN=rep(NA_real_, ngenes), 
                                    row.names=geneIDs)
trial.samr.qval.table$U24U49[match(expr.data.dividedByU24U49.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24U49.samr.siggenes.table$genes.up[,8]) / 100.0
trial.samr.qval.table$U24U49[match(expr.data.dividedByU24U49.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24U49.samr.siggenes.table$genes.lo[,8]) / 100.0
trial.samr.qval.table$U49[match(expr.data.dividedByU49.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU49.samr.siggenes.table$genes.up[,8]) / 100.0
trial.samr.qval.table$U49[match(expr.data.dividedByU49.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU49.samr.siggenes.table$genes.lo[,8]) / 100.0
#
trial.samr.qval.table$U24U49.log[match(expr.data.dividedByU24U49.lumiT.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24U49.lumiT.samr.siggenes.table$genes.up[,8]) / 100.0
trial.samr.qval.table$U24U49.log[match(expr.data.dividedByU24U49.lumiT.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU24U49.lumiT.samr.siggenes.table$genes.lo[,8]) / 100.0
trial.samr.qval.table$U49.log[match(expr.data.dividedByU49.lumiT.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU49.lumiT.samr.siggenes.table$genes.up[,8]) / 100.0
trial.samr.qval.table$U49.log[match(expr.data.dividedByU49.lumiT.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.dividedByU49.lumiT.samr.siggenes.table$genes.lo[,8]) / 100.0
#
trial.samr.qval.table$VST.QN.ComBat.exp[match(expr.data.qn.combat.exp.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.qn.combat.exp.samr.siggenes.table$genes.up[,8]) / 100.0
trial.samr.qval.table$VST.QN.ComBat.exp[match(expr.data.qn.combat.exp.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.qn.combat.exp.samr.siggenes.table$genes.lo[,8]) / 100.0
trial.samr.qval.table$VST.QN.exp[match(expr.data.qn.exp.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.qn.exp.samr.siggenes.table$genes.up[,8]) / 100.0
trial.samr.qval.table$VST.QN.exp[match(expr.data.qn.exp.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.qn.exp.samr.siggenes.table$genes.lo[,8]) / 100.0
trial.samr.qval.table$VST.QN[match(expr.data.qn.samr.siggenes.table$genes.up[,2] ,geneIDs)] <- as.numeric(expr.data.qn.samr.siggenes.table$genes.up[,8]) / 100.0
trial.samr.qval.table$VST.QN[match(expr.data.qn.samr.siggenes.table$genes.lo[,2] ,geneIDs)] <- as.numeric(expr.data.qn.samr.siggenes.table$genes.lo[,8]) / 100.0
##
write.table(trial.samr.qval.table, "PBMC.211samp.trial.samr.qval.txt", quote=FALSE, sep='\t')

# END ##

# using wilcox option (the same input as above taking the same order)

# ref gene normalised unlogged values 
expr.data.dividedByU24U49.samr.wcx.data <- list(x=expr.data.dividedByU24U49[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=FALSE)
expr.data.dividedByU24U49.samr.wcx.obj <- samr(expr.data.dividedByU24U49.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=5000)
expr.data.dividedByU24U49.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.dividedByU24U49.samr.wcx.obj, nvals=50)
expr.data.dividedByU24U49.samr.wcx.siggenes.table <- samr.compute.siggenes.table(expr.data.dividedByU24U49.samr.wcx.obj, delta, expr.data.dividedByU24U49.samr.wcx.data, expr.data.dividedByU24U49.samr.wcx.delta.table)

expr.data.dividedByU49.samr.wcx.data <- list(x=expr.data.dividedByU49[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=FALSE)
expr.data.dividedByU49.samr.wcx.obj <- samr(expr.data.dividedByU49.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=5000)
expr.data.dividedByU49.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.dividedByU49.samr.wcx.obj, nvals=50)
expr.data.dividedByU49.samr.wcx.siggenes.table <- samr.compute.siggenes.table(expr.data.dividedByU49.samr.wcx.obj, delta, expr.data.dividedByU49.samr.wcx.data, expr.data.dividedByU49.samr.wcx.delta.table)

# ref gene normalised logged un-combatted values 
expr.data.dividedByU24U49.lumiT.samr.wcx.data <- list(x=expr.data.dividedByU24U49.lumiT[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU24U49.lumiT.samr.wcx.obj <- samr(expr.data.dividedByU24U49.lumiT.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=5000)
expr.data.dividedByU24U49.lumiT.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.dividedByU24U49.lumiT.samr.wcx.obj, nvals=50)
expr.data.dividedByU24U49.lumiT.samr.wcx.siggenes.table <- samr.compute.siggenes.table(expr.data.dividedByU24U49.lumiT.samr.wcx.obj, delta, expr.data.dividedByU24U49.lumiT.samr.wcx.data, expr.data.dividedByU24U49.lumiT.samr.wcx.delta.table)

expr.data.dividedByU49.lumiT.samr.wcx.data <- list(x=expr.data.dividedByU49.lumiT[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=TRUE)
expr.data.dividedByU49.lumiT.samr.wcx.obj <- samr(expr.data.dividedByU49.lumiT.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=5000)
expr.data.dividedByU49.lumiT.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.dividedByU49.lumiT.samr.wcx.obj, nvals=50)
expr.data.dividedByU49.lumiT.samr.wcx.siggenes.table <- samr.compute.siggenes.table(expr.data.dividedByU49.lumiT.samr.wcx.obj, delta, expr.data.dividedByU49.lumiT.samr.wcx.data, expr.data.dividedByU49.lumiT.samr.wcx.delta.table)


# exp VST+QN+COMBAT value
expr.data.qn.combat.exp <- exp(expr.data.qn.combat)
expr.data.qn.combat.exp.samr.wcx.data <- list(x=expr.data.qn.combat.exp, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=FALSE)
expr.data.qn.combat.exp.samr.wcx.obj <- samr(expr.data.qn.combat.exp.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=5000)
expr.data.qn.combat.exp.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.qn.combat.exp.samr.wcx.obj, nvals=50)
expr.data.qn.combat.exp.samr.wcx.siggenes.table <- samr.compute.siggenes.table(expr.data.qn.combat.exp.samr.wcx.obj, delta, expr.data.qn.combat.exp.samr.wcx.data, expr.data.qn.combat.exp.samr.wcx.delta.table)

# exp VST+QN value
expr.data.qn.exp <- exp(expr.data.qn[,labeled])
expr.data.qn.exp.samr.wcx.data <- list(x=expr.data.qn.exp, y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=FALSE)
expr.data.qn.exp.samr.wcx.obj <- samr(expr.data.qn.exp.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=5000)
expr.data.qn.exp.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.qn.exp.samr.wcx.obj, nvals=50)
expr.data.qn.exp.samr.wcx.siggenes.table <- samr.compute.siggenes.table(expr.data.qn.exp.samr.wcx.obj, delta, expr.data.qn.exp.samr.wcx.data, expr.data.qn.exp.samr.wcx.delta.table)

# VST+QN value
expr.data.qn.samr.wcx.data <- list(x=expr.data.qn[,labeled], y=label[labeled], geneid=geneIDs, genenames=geneIDs, logged2=FALSE)
expr.data.qn.samr.wcx.obj <- samr(expr.data.qn.samr.wcx.data, resp.type="Two class unpaired", testStatistic="wilcoxon", nperms=5000)
expr.data.qn.samr.wcx.delta.table <- samr.compute.delta.table(expr.data.qn.samr.wcx.obj, nvals=50)
expr.data.qn.samr.wcx.siggenes.table <- samr.compute.siggenes.table(expr.data.qn.samr.wcx.obj, delta, expr.data.qn.samr.wcx.data, expr.data.qn.samr.wcx.delta.table)
