disp <-read.table("liver.disp.1642962.summary")

dispnasc <- na.omit(disp)

dispnascmat <- dist(dispnasc, method ="euclidean")
dispnascmat <- na.omit(dispnascmat)

fit <- hclust(dispnascmat, method="ward")

library(ggplot2)
#install.packages("ggdendro")
library(ggdendro)

fit$labels <- dispnasc$V1
gf<-ggdendrogram(fit)


pdf("SRR1642962.pdf", paper="a4r", width=15, height=4)
plot(gf)
dev.off()