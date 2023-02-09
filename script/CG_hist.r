#! usr/bin/R
# Rscript CG_hist.r dir PREFIX

args <- commandArgs(T)
library("ggplot2")
FILE <- paste(sep="", args[1], "/filter/", args[2], "/plot_CG.tsv")
DATA <- read.table(FILE, header = TRUE, sep = "\t")
LAB <- paste(sep="", args[2], "_", args[1], "_CG")

ggplot(DATA, aes(CG_content, number, colour=group)) +
geom_line()+
xlab("CG_content(%)")+
scale_x_continuous(breaks = seq(25, 75, 5))+
labs(title = LAB)+
theme(text=element_text(face = "bold"), axis.text=element_text(face = "bold"), plot.title = element_text(hjust=0.5))

DIR <- paste(sep="", "../visuliaztion/plot/", args[2], "/", LAB, ".png")
ggsave(file=DIR, width=7.5, height=5, dpi=300)
