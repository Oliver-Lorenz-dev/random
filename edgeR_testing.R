library(edgeR)
library(readr)
library(dplyr)
library(magrittr)
library(tximport)
library(statmod)

x <- read.delim("rna-seq-chris-analysis/input_data/salmon.merged.gene_counts.txt")

drops <- c("gene_id")
x = x[ , !(names(x) %in% drops)]

# just keep 2 conditions for testing purposes
x <- x[ -c( 2:21, 26:28) ]

names = x[,1]
rownames(x) = make.names(names, unique = TRUE)

x_final = data.matrix(x)
x_final = x_final[, colnames(x_final) != "Symbol"]
group <- factor(c(1,1,2,2))
y <- DGEList(counts=x_final,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

logcpm <- cpm(y, log=TRUE)
plotMDS(y)

Cell <- factor(c(3,4,3,4))
Toxin <- factor(c("Y","Y","N","N"))
data.frame(Sample=colnames(y),Cell,Toxin)
design <- model.matrix(~Cell+Toxin)
rownames(design) <- colnames(y)

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)

plotMD(lrt)
abline(h=c(-1, 1), col="blue")
