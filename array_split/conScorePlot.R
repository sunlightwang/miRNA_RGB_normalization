require(ggplot2)
theme_set(theme_grey(20)) # set theme with grey background and white gridlines, base size = 18

consist_score_errorbar_plot <- function(mean, se) {
  stopifnot( all(dim(mean) == dim(se)) ) 
  x <- 1 : ncol(mean)
  c <- 1 : nrow(mean)
  x_lab <- c("5","10","20","50","100","150","200","250")
  c_lab <- c("VST+QN","VST+RSN","VST+SSN","VST+VSN","VST+RIN","VST+Loess","Log+U49","Log+U24","Log+U49U24","VST+U49","VST+U24","VST+U49U24")
  
  m <- rep(as.factor(c), each=length(x))
  t <- rep(x, nrow(mean))
  mean_vec <- c(t(mean))
  se_vec <- c(t(se))
  data <- data.frame(methods=m, cbind(topN=t, score=mean_vec, stde=se_vec))
  
  pd <- position_dodge(.6) # move them .05 to the left and right
  
  ggplot(data, aes(x=topN, y=score, colour=methods, group=methods)) + 
    geom_errorbar(aes(ymin=score-stde, ymax=score+stde), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size=3, shape=22, fill="white") +
    scale_x_continuous(breaks = x, labels = x_lab) + 
    ylab("Consistency Score") + 
    xlab("Top N Genes") + 
    scale_colour_hue(name="Norm Methods", breaks=c, labels=c_lab, l=40) + 
    ggtitle("Consistency Score Comparison") + 
    theme_bw(base_size = 18) +
    theme(legend.justification=c(0,1), legend.position=c(0,1)) # Position legend in bottom right
}


consist_score_errorbar_barplot <- function(mean, se, main="") {
  stopifnot( all(dim(mean) == dim(se)) ) 
  x <- 1 : ncol(mean)
  c <- 1 : nrow(mean)
  x_lab <- c("5","10","20","50","100","150","200","250")
  c_lab <- c("VST+QN","VST+RSN","VST+SSN","VST+VSN","VST+RIN","VST+Loess","Log+U49","Log+U24","Log+U49U24","VST+U49","VST+U24","VST+U49U24")
  
  m <- rep(factor(c_lab, levels=c_lab, ordered =TRUE), each=length(x))
  t <- rep(x, nrow(mean))
  mean_vec <- c(t(mean))
  se_vec <- c(t(se))
  data <- data.frame(Methods=m, cbind(topN=t, score=mean_vec, stde=se_vec))
  
  pd <- position_dodge(0.9)
  cols <- c("VST+QN"="darkolivegreen1","VST+RSN"="darkolivegreen2","VST+SSN"="darkolivegreen3",
            "VST+VSN"="springgreen1","VST+RIN"="springgreen2","VST+Loess"="springgreen3",
            "Log+U49"="royalblue1","Log+U24"="royalblue2","Log+U49U24"="royalblue3",
            "VST+U49"="purple1","VST+U24"="purple3","VST+U49U24"="purple4")
  
  ggplot(data, aes(x=topN, y=score, group=Methods)) + 
    geom_bar(position=pd,  stat="identity", colour="black", aes(fill=Methods)) + 
    #     scale_fill_grey(start=0.8,end=0.2) + 
    scale_fill_manual(values = cols) + 
    geom_errorbar(aes(ymin=score-stde, ymax=score+stde), colour="black", width=.45, position=pd) +
    #     geom_line(position=pd) +
    scale_x_continuous(breaks = x, labels = x_lab) + 
    scale_y_sqrt(breaks=c(1,10,50,100,200,300)) + 
    ylab("Consistency Score") + 
    xlab("Top N DE Genes") + 
    #     scale_fill_hue(name="Norm Methods", breaks=c, labels=c_lab, l=40) + 
    ggtitle(main) + 
    theme_bw(base_size = 16) +
    theme(legend.justification=c(0,1), legend.position=c(0,1)) # Position legend in bottom right
}

setwd("~/Projects/bioMarker/PBMCFromRawData/array_split/")
se <- 0
# limma
mean <- read.table("conscore.limma")
consist_score_errorbar_barplot(mean, se, main="limma")
# ftest
mean <- read.table("conscore.ftest")
consist_score_errorbar_barplot(mean, se, main="F-test")
# samr
mean <- read.table("conscore.samr")
consist_score_errorbar_barplot(mean, se, main="SAM")

##### combat #####
# limma
mean <- read.table("conscore.combat.limma")
consist_score_errorbar_barplot(mean, se, main="limma with Combat")
# ftest
mean <- read.table("conscore.combat.ftest")
consist_score_errorbar_barplot(mean, se, main="F-test with Combat")
# samr
mean <- read.table("conscore.combat.samr")
consist_score_errorbar_barplot(mean, se, main="SAM with Combat")
