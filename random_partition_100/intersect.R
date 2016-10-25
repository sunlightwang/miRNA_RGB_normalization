i <- 0:99
dirs <- paste("RUN", i, "/", sep="")
keys <- c("A001","A002","A003","A004")
outfile <- "limma.intersect.250"

intersection <- function(x, y, ...){
 if (missing(...)) intersect(x, y)
 else intersect(x, intersection(y, ...))
}

DEresultAnalysis <- function(dir, keys, norm.method="vst.U24", DEG.method="limma", top=250) {
  require(Vennerable)
  res <- list()
  n <- top
  tops <- list()
  for(i in 1:length(keys)) {
    res[[i]] <- read.table(paste(paste(dir, keys[i], sep=""), norm.method, DEG.method, "rst", sep="."), header=T, sep="\t")
    tops[[i]] <- res[[i]][1:n,1]
  }
  intersection(tops[[1]], tops[[2]], tops[[3]], tops[[4]])
}

if(file.exists(outfile)) {file.remove(outfile)} 
for(i in 1:length(dirs)) {
  write.table(DEresultAnalysis(dirs[i], keys), outfile, append=T, quote=FALSE, col.names=F, row.names=F, sep="\t")
  write("\n", outfile, append=T)
}
