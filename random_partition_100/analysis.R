i <- 0:99
dir <- paste("RUN", i, sep="")

# for conscore.limma
x <- array(0, dim=c(100,12,8))
for(i in 1:100) {
  x[i,,] <- as.matrix(read.table(paste(dir[i],"/conscore.limma", sep="")))
}
write.table(apply(x, c(2,3), mean), "conscore.limma.mean", row.names=F, col.names=F, quote=F, sep="\t")
write.table(sqrt(apply(x, c(2,3), var)/100), "conscore.limma.se", row.names=F, col.names=F, quote=F, sep="\t")

# conscore.ftest
x <- array(0, dim=c(100,12,8))
for(i in 1:100) {
  x[i,,] <- as.matrix(read.table(paste(dir[i],"/conscore.ftest", sep="")))
}
write.table(apply(x, c(2,3), mean), "conscore.ftest.mean", row.names=F, col.names=F, quote=F, sep="\t")
write.table(sqrt(apply(x, c(2,3), var)/100), "conscore.ftest.se", row.names=F, col.names=F, quote=F, sep="\t")

# conscore.samr
x <- array(0, dim=c(100,12,8))
for(i in 1:100) {
  x[i,,] <- as.matrix(read.table(paste(dir[i],"/conscore.samr", sep="")))
}
write.table(apply(x, c(2,3), mean), "conscore.samr.mean", row.names=F, col.names=F, quote=F, sep="\t")
write.table(sqrt(apply(x, c(2,3), var)/100), "conscore.samr.se", row.names=F, col.names=F, quote=F, sep="\t")

##########################################################
# for conscore.combat.limma
x <- array(0, dim=c(100,12,8))
for(i in 1:100) {
  x[i,,] <- as.matrix(read.table(paste(dir[i],"/conscore.combat.limma", sep="")))
}
write.table(apply(x, c(2,3), mean), "conscore.combat.limma.mean", row.names=F, col.names=F, quote=F, sep="\t")
write.table(sqrt(apply(x, c(2,3), var)/100), "conscore.combat.limma.se", row.names=F, col.names=F, quote=F, sep="\t")

# conscore.combat.ftest
x <- array(0, dim=c(100,12,8))
for(i in 1:100) {
  x[i,,] <- as.matrix(read.table(paste(dir[i],"/conscore.combat.ftest", sep="")))
}
write.table(apply(x, c(2,3), mean), "conscore.combat.ftest.mean", row.names=F, col.names=F, quote=F, sep="\t")
write.table(sqrt(apply(x, c(2,3), var)/100), "conscore.combat.ftest.se", row.names=F, col.names=F, quote=F, sep="\t")

# conscore.combat.samr
x <- array(0, dim=c(100,12,8))
for(i in 1:100) {
  x[i,,] <- as.matrix(read.table(paste(dir[i],"/conscore.combat.samr", sep="")))
}
write.table(apply(x, c(2,3), mean), "conscore.combat.samr.mean", row.names=F, col.names=F, quote=F, sep="\t")
write.table(sqrt(apply(x, c(2,3), var)/100), "conscore.combat.samr.se", row.names=F, col.names=F, quote=F, sep="\t")
