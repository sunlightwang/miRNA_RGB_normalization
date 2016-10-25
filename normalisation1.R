# analysis PBMC miRNA data 
# the raw data were exported from GenomeStudio 

data <- read.table("~/Projects/miRNA-mRNA/miRNA_array/PBMC/original/miRNA rehash 17-1-2011_Group_Gene_Profile_raw_bg_sub.txt", 
                   skip=7, header=TRUE)
nc <- ncol(data)
expr <- data[,seq(3,nc,8)]
nsample <- ncol(expr)
ngene <- nrow(expr)
rownames(expr) <- data[,1]
samples <- colnames(expr)
samples <- sapply(samples, function(x) {
  t <- unlist(strsplit(x,"\\."))
  t[2]
})
arrays <- sapply(samples, function(x) {
  t <- unlist(strsplit(x, "_"))
  t[1]
})
arrays <- as.factor(as.vector(arrays))
sample.mean <- colMeans(expr)
sample.stdev <- sqrt( (colMeans(expr^2) - (sample.mean)^2) * ngene / (ngene - 1) )

data <- read.table("~/Projects/miRNA-mRNA/miRNA_array/PBMC/original/miRNA rehash 17-1-2011 raw_bg_sub Control Probe Profile.txt", 
                      sep="\t", header=TRUE)
ctrl <- as.matrix(data[,-(1:2)])
ctrl.meta <- data[,1:2]
samples.ctrl <- colnames(ctrl)
samples.ctrl <- sapply(samples.ctrl, function(x) {
  t <- unlist(strsplit(x,"\\."))
  t <- unlist(strsplit(t[1],"X"))
  t[2]
})
# for ctrl data
# U49: 2011, row 52
# U66: 3792, row 11
# U24: 1992, row 9

pdf("controls.pdf")
plot(1:nsample, sample.mean * 30 , type="l", col=1, xlab="samples", ylab="expr", main="all samples")
points(1:nsample, rep(0,nsample), col=arrays, pch=20)
points(1:nsample, sample.stdev * 10 , type="l", col=2)
points(1:nsample, ctrl[52,], type="l", col=3)
points(1:nsample, ctrl[11,], type="l", col=4)
points(1:nsample, ctrl[9,], type="l", col=5)
points(1:nsample, sqrt(ctrl[11,] * ctrl[9,])*2, type="l", col=6)
leg.txt <- c("expr mean", "sample stdev","U49", "U66", "U24","geo mean(U66, U24)")
legend("topright", bty="n", bg = "white", legend=leg.txt, 
       lty = rep(1,6), lwd = rep(2,6), col = 1:6, cex = 0.9)
dev.off()

# dealing with sample labels
##label188 <- read.table("sampleLabel188.txt")
label211 <- read.table("sampleLabel211.txt", header=T)
##match <- label188[match(label211[,1], label188[,1]),4] == label211[,4]
##match[is.na(match)] <- FALSE
##stopifnot( sum(match) == 188)
# using label211 below
label <- label211[match(samples, label211[,1]),4]
label[is.na(label)] <- "NA"
SZ <- label == "SZ"
CTRL <- label == "CTRL"
SZ.expr <- expr[,SZ] 
SZ.ctrl <- ctrl[,SZ]
CTRL.expr <- expr[,CTRL] 
CTRL.ctrl <- ctrl[,CTRL] 
nSZ <- sum(SZ)
nCTRL <- sum(CTRL)

pdf("SZ.controls.pdf")
plot(1:nSZ, sample.mean[SZ] * 30 , type="l", col=1, xlab="samples", ylab="expr", ylim=range(sample.mean*30), 
     main="Diagnosis=SZ")
points(1:nSZ, rep(0,nSZ), col=arrays[SZ], pch=20)
points(1:nSZ, sample.stdev[SZ] * 10 , type="l", col=2)
points(1:nSZ, ctrl[52,SZ], type="l", col=3)
points(1:nSZ, ctrl[11,SZ], type="l", col=4)
points(1:nSZ, ctrl[9,SZ], type="l", col=5)
points(1:nSZ, sqrt(ctrl[11,SZ] * ctrl[9,SZ])*2, type="l", col=6)
leg.txt <- c("expr mean", "sample stdev","U49", "U66", "U24","geo mean(U66, U24)")
legend("topright", bty="n", bg = "white", legend=leg.txt, 
       lty = rep(1,6), lwd = rep(2,6), col = 1:6, cex = 0.9)
dev.off()

pdf("CTRL.controls.pdf")
plot(1:nCTRL, sample.mean[CTRL] * 30 , type="l", col=1, xlab="samples", ylab="expr", ylim=range(sample.mean*30), 
     main="Diagnosis=CTRL")
points(1:nCTRL, rep(0,nCTRL), col=arrays[CTRL], pch=20)
points(1:nCTRL, sample.stdev[CTRL] * 10 , type="l", col=2)
points(1:nCTRL, ctrl[52,CTRL], type="l", col=3)
points(1:nCTRL, ctrl[11,CTRL], type="l", col=4)
points(1:nCTRL, ctrl[9,CTRL], type="l", col=5)
points(1:nCTRL, sqrt(ctrl[11,CTRL] * ctrl[9,CTRL])*2, type="l", col=6)
leg.txt <- c("expr mean", "sample stdev","U49", "U66", "U24","geo mean(U66, U24)")
legend("topright", bty="n", bg = "white", legend=leg.txt, 
       lty = rep(1,6), lwd = rep(2,6), col = 1:6, cex = 0.9)
dev.off()


