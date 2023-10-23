dat <- read.table("tamakiarray_hclust.txt", header=T, row.names=1, sep="\t")
data <- dat[, -which(colnames(dat) %in% c("empty.vector2","empty.vector4") )]
distance <- dist(t(data), method="euclidean")
hc <- hclust(distance, method="ward.D2")

pdf("hclust.pdf", width=4, height=5)
plot(hc, xlab="", ylab="distance (intensity)")
dev.off()
