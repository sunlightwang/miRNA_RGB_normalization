args <- commandArgs(TRUE)
if(length(args) < 1 || length(args) > 1)
{
  cat ("USAGE: Rscript split_samples.R <dir>\n")
  quit("no")
}

dir <- paste(args[1], "/", sep="")

dir.create(dir)

expr <- read.table(file="./miRNA rehash 17-1-2011_Group_Gene_Profile_raw_bg_sub.txt", skip=7, header=T)
ctrl <- read.table(file="./miRNA rehash 17-1-2011 raw_bg_sub Control Probe Profile.txt", header=T, sep="\t")
sample <- read.table("./miRNA rehash 17-1-2011 raw_bg_sub Samples Table.txt", header=T, sep="\t")

spliter <- seq(2, ncol(expr), by=8)
sampler <- sample(1:length(spliter), length(spliter))
spliter <- spliter[sampler]

write.table(matrix(sampler[seq(1,length(spliter),by=4)], nrow=1), paste(dir, "sample.matrix.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.name=FALSE)
write.table(matrix(sampler[seq(2,length(spliter),by=4)], nrow=1), paste(dir, "sample.matrix.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.name=FALSE, append=TRUE)
write.table(matrix(sampler[seq(3,length(spliter),by=4)], nrow=1), paste(dir, "sample.matrix.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.name=FALSE, append=TRUE)
write.table(matrix(sampler[seq(4,length(spliter),by=4)], nrow=1), paste(dir, "sample.matrix.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE, col.name=FALSE, append=TRUE)

cols1 <- c(sapply(spliter[seq(1,length(spliter),by=4)], function(x) x:(x+7)))
cols2 <- c(sapply(spliter[seq(2,length(spliter),by=4)], function(x) x:(x+7)))
cols3 <- c(sapply(spliter[seq(3,length(spliter),by=4)], function(x) x:(x+7)))
cols4 <- c(sapply(spliter[seq(4,length(spliter),by=4)], function(x) x:(x+7)))
expr1 <- expr[,c(1,cols1)]
expr2 <- expr[,c(1,cols2)]
expr3 <- expr[,c(1,cols3)]
expr4 <- expr[,c(1,cols4)]

ctrl1 <- ctrl[,c(1,2,2+sampler[seq(1,length(spliter),by=4)])]
ctrl2 <- ctrl[,c(1,2,2+sampler[seq(2,length(spliter),by=4)])]
ctrl3 <- ctrl[,c(1,2,2+sampler[seq(3,length(spliter),by=4)])]
ctrl4 <- ctrl[,c(1,2,2+sampler[seq(4,length(spliter),by=4)])]

sample1 <- sample[sampler[seq(1,length(spliter),by=4)],]
sample2 <- sample[sampler[seq(2,length(spliter),by=4)],]
sample3 <- sample[sampler[seq(3,length(spliter),by=4)],]
sample4 <- sample[sampler[seq(4,length(spliter),by=4)],]

write.table(expr1, paste(dir, "A001_expr.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
write.table(expr2, paste(dir, "A002_expr.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
write.table(expr3, paste(dir, "A003_expr.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
write.table(expr4, paste(dir, "A004_expr.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)

write.table(ctrl1, paste(dir, "A001_ctrl.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
write.table(ctrl2, paste(dir, "A002_ctrl.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
write.table(ctrl3, paste(dir, "A003_ctrl.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
write.table(ctrl4, paste(dir, "A004_ctrl.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)

write.table(sample1, paste(dir, "A001_sample.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
write.table(sample2, paste(dir, "A002_sample.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
write.table(sample3, paste(dir, "A003_sample.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
write.table(sample4, paste(dir, "A004_sample.txt", sep=""), quote=FALSE, sep="\t", row.names=FALSE)
