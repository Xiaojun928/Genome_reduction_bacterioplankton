setwd("E:/CHUG_reduction/MA/synthesis")
library(dplyr)
# Vertical box plot
df <- read.table("Roseo_gnmsize.tsv",sep="\t",header = T)
df <- df[-c(1:4),]
qdf <- df %>% filter(Completeness >= 90 & Contamination <= 5)
qdf$Msize <- as.numeric(qdf$estimated.size) / 1000000
boxplot(qdf$Msize, col = "white")

