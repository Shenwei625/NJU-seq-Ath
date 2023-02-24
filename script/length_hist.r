#! usr/bin/R
# Rscript length_hist.r dir PREFIX

args <- commandArgs(T)
library("ggplot2")
FILE <- paste(sep="", args[1], "/filter/", args[2], "/plot.tsv")
LAB <- paste(sep="", args[2], "_", args[1])

ggplot(DATA, aes(reads_length, number, fill=group)) +
geom_bar(stat = "identity")+
scale_x_continuous(breaks = seq(10, 50, 5))+
labs(title = LAB)+
theme(text=element_text(face = "bold", size = 8), axis.text=element_text(face = "bold", size = 8), plot.title = element_text(hjust=0.5))

DIR <- paste(sep="", "../visuliaztion/plot/", args[2], "/", LAB, ".png")
ggsave(file=DIR, width=5, height=3.5, dpi=300)
