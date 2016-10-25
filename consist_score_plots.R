require(gplots)
require(ggplot2)
theme_set(theme_grey(20)) # set theme with grey background and white gridlines, base size = 18

consist_score_plot <- function(data) {
  x <- 1 : ncol(data)
  x_lab <- c("5","10","20","50","100","150","200","250")
  q <- qplot(ylab="Consistency score", xlab="Number of top genes") #+ scale_y_log10()
  q <- q + 
    geom_point(aes(x=x, y=data[1,]), color="yellow", pch=15, size=3) + 
    geom_point(aes(x=x, y=data[2,]), color="yellow", pch=16, size=3) + 
    geom_point(aes(x=x, y=data[3,]), color="yellow", pch=17, size=3) + 
    geom_point(aes(x=x, y=data[4,]), color="green", pch=15, size=3) + 
    geom_point(aes(x=x, y=data[5,]), color="green", pch=16, size=3) + 
    geom_point(aes(x=x, y=data[6,]), color="green", pch=17, size=3) + 
    geom_point(aes(x=x, y=data[7,]), color="red", pch=15, size=3) + 
    geom_point(aes(x=x, y=data[8,]), color="red", pch=16, size=3) + 
    geom_point(aes(x=x, y=data[9,]), color="red", pch=17, size=3) + 
    geom_point(aes(x=x, y=data[10,]), color="purple", pch=15, size=3) + 
    geom_point(aes(x=x, y=data[11,]), color="purple", pch=16, size=3) + 
    geom_point(aes(x=x, y=data[12,]), color="purple", pch=17, size=3) 
  q <- q + 
    geom_line(aes(x=x, y=data[1,]), color="yellow", linetype=3, size=1) + 
    geom_line(aes(x=x, y=data[2,]), color="yellow", linetype=2, size=1) + 
    geom_line(aes(x=x, y=data[3,]), color="yellow", linetype=4, size=1) + 
    geom_line(aes(x=x, y=data[4,]), color="green", linetype=3, size=1) + 
    geom_line(aes(x=x, y=data[5,]), color="green", linetype=2, size=1) + 
    geom_line(aes(x=x, y=data[6,]), color="green", linetype=4, size=1) + 
    geom_line(aes(x=x, y=data[7,]), color="red", linetype=3, size=1) + 
    geom_line(aes(x=x, y=data[8,]), color="red", linetype=2, size=1) + 
    geom_line(aes(x=x, y=data[9,]), color="red", linetype=4, size=1) + 
    geom_line(aes(x=x, y=data[10,]), color="purple", linetype=3, size=1) + 
    geom_line(aes(x=x, y=data[11,]), color="purple", linetype=2, size=1) + 
    geom_line(aes(x=x, y=data[12,]), color="purple", linetype=4, size=1) 
  q <- q + scale_x_discrete(breaks = x, labels = x_lab)
  q
}

setwd("~/Projects/bioMarker/PBMCFromRawData")
#1
four_batch_data <- as.matrix(read.table("4batch_cs.txt", header=T, row.names=1))
consist_score_plot(four_batch_data)
#2
two_comb_four_batch_data <- as.matrix(read.table("2comb_4batch_cs.txt", header=T, row.names=1))
consist_score_plot(two_comb_four_batch_data)
#3
two_comb_four_batch_data_combat <- as.matrix(read.table("2comb_4batch_combat_cs.txt", header=T, row.names=1))
consist_score_plot(two_comb_four_batch_data_combat)
#4
three_comb_four_batch_data <- as.matrix(read.table("3comb_4batch_cs.txt", header=T, row.names=1))
consist_score_plot(three_comb_four_batch_data)
#5
three_comb_four_batch_data_combat <- as.matrix(read.table("3comb_4batch_combat_cs.txt", header=T, row.names=1))
consist_score_plot(three_comb_four_batch_data_combat)
#6
random_half_data <- as.matrix(read.table("random_half_cs.txt", header=T, row.names=1)) 
consist_score_plot(random_half_data)
#7 
random_half_data_combat <- as.matrix(read.table("random_half_combat_cs.txt", header=T, row.names=1)) 
consist_score_plot(random_half_data_combat)
#8
random_quater_data <- as.matrix(read.table("random_quater_cs.txt", header=T, row.names=1)) 
consist_score_plot(random_quater_data)
#9 
random_quater_data_combat <- as.matrix(read.table("random_quater_combat_cs.txt", header=T, row.names=1)) 
consist_score_plot(random_quater_data_combat)


legend.txt <- c("QN", "RSN", "SSN", "VSN", "RIN", "Loess", "Log+U49", "Log+U24", "Log+U49U24", "VST+U49", "VST+U24", "VST+U49U24")
legend.col <- c("yellow", "yellow", "yellow", "green", "green", "green", "red", "red", "red", "purple", "purple","purple")
legend.pch <- rep(c(15, 16, 17), 4)
legend.lty <- rep(c(3,2,4), 4)
plot.new()
legend("center",legend.txt, lty=legend.lty, pch=legend.pch, col=legend.col, bty="o", bg="grey92")