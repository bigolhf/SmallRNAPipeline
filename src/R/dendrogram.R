disp <-read.table("liver.disp.summary")

dispnasc <- na.omit(disp)

dispnascmat <- dist(dispnasc, method ="euclidean")
dispnascmat <- na.omit(dispnascmat)
print("start clustering")
fit <- hclust(dispnascmat, method="ward")

library(ggplot2)
#install.packages("ggdendro")
library(ggdendro)

fit$labels <- dispnasc$V1
print("generate dendrogram")
gf<-ggdendrogram(fit)

print("plotting")
pdf("liver.pdf", width=500, height=20)
plot(gf)
dev.off()
print("done")