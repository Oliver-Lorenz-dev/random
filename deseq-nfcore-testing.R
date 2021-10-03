# DESeq2 differential gene expression analysis code

library(BiocManager)
library(DESeq2)
library(readr)
library(dplyr)
library(magrittr)
library(tximport)

run_table = read_csv('./demo_run_table.txt')

# read in read count data using tximport
count_data = read_tsv('./salmon.merged.gene_counts.tsv')

count_data = data.frame(count_data)

count_data = count_data[,-1]
rownames(count_data) = count_data[,1]

count_data_final = subset(count_data, select = -c(1))

count_data_dds <- mutate_all(count_data_final, function(x) as.integer(as.character(x)))

deseq_data <- DESeqDataSetFromMatrix(countData = count_data_dds,
                              colData = run_table,
                              design= ~ Condition)

# run DESeq2 analysis on data
dds = DESeq(deseq_data)

# get results
dds_results = results(dds , contrast = c("Condition","WT","UNINDUCED"))

# check each gene for differential expression
dds_results$dif_exp = dds_results$padj < 0.05 & abs(dds_results$log2FoldChange) > 1
dds_dif_exp_results = data.frame(dds_results)

# use dplyr to filter dataframe for differentially expressed genes only
dds_dif_exp_results = filter(dds_dif_exp_results, padj < 0.05)
dds_dif_exp_results = filter(dds_dif_exp_results, abs(log2FoldChange) > 1)

# MA plot
plotMA(dds , alpha = 0.01, main = "MA plot")

# PCA plot
dds_var_transform = varianceStabilizingTransformation(deseq_data)
plotPCA(dds_var_transform , intgroup = 'Condition')
