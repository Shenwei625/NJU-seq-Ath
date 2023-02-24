#! usr/bin/R
# Rscript end_count_hist.r ${TISSUE} rrna_5-8s/18s/25s/mrna

args <- commandArgs(T)
library("ggplot2")
library("patchwork")

FILE1 <- paste(sep="", "visuliaztion/End_count_hist/", args[1], "_", args[2], "_NC.tsv")
FILE2 <- paste(sep="", "visuliaztion/End_count_hist/", args[1], "_", args[2], "_TR.tsv")
TITLE1 <- paste(sep=" ", args[1], args[2], "NC")
TITLE2 <- paste(sep=" ", args[1], args[2], "TR")

DATA1 <- read.table(FILE1, header = TRUE, sep = "\t")
DATA2 <- read.table(FILE2, header = TRUE, sep = "\t")

P1 <- ggplot(DATA1, aes(x = Site, y = End_count)) +
geom_bar(stat = "identity") +
labs(title = TITLE1) +
theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank())

P2 <- ggplot(DATA2, aes(x = Site, y = End_count)) +
labs(title = TITLE2) +
geom_bar(stat = "identity")

P1/P2
DIR <- paste(sep="", "visuliaztion/End_count_hist/", args[1], "_", args[2], ".png")
ggsave(file=DIR, width=20, height=10, dpi=300)
