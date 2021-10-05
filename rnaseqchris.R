# DESeq2 differential gene expression analysis code for nfcore-rnaseq results

library(BiocManager)
library(DESeq2)
library(readr)
library(dplyr)
library(magrittr)
library(tximport)

run_table = read_csv('./run_table_chris.txt')

# extra filtering (if needed)
run_table_ed = filter(run_table, p_value == '35')

# filter run table for values you want
run_table_1 = filter(run_table_ed, toxin == 'CTX')
run_table_2 = filter(run_table_ed, toxin == 'ctrl')


run_table_final = full_join(run_table_1,run_table_2)

# read in read count data
count_data = read_tsv('./salmon.merged.gene_counts.tsv')

count_data = data.frame(count_data)
names = count_data[,2]

rownames(count_data) = make.names(names, unique = TRUE)
# drop unnecessary columns and set row names to gene ids
count_data_final = subset(count_data, select = -c(1,2))

# convert values to integer
count_data_dds <- mutate_all(count_data_final, function(x) as.integer(as.character(x)))

# filter count data to match run table
sample_list = run_table_final
sample_list = data.frame(sample_list)
rownames(sample_list) = sample_list[,1]

final_samples = rownames(sample_list)
count_data_filter = count_data_dds[final_samples]



deseq_data <- DESeqDataSetFromMatrix(countData = count_data_filter,
                                     colData = run_table_final,
                                     design= ~ toxin)

# run DESeq2 analysis on data
dds = DESeq(deseq_data)

# get results
dds_results = results(dds , contrast = c("toxin","CTX","ctrl"))

# check each gene for differential expression
dds_results$dif_exp = dds_results$padj < 0.05 & abs(dds_results$log2FoldChange) > 1
dds_dif_exp_results = data.frame(dds_results)

# use dplyr to filter dataframe for differentially expressed genes only
dds_dif_exp_results = filter(dds_dif_exp_results, padj < 0.05)
dds_dif_exp_results = filter(dds_dif_exp_results, abs(log2FoldChange) > 1)

# MA plot
plotMA(dds , alpha = 0.01, main = "P35 - Toxin: CTX v ctrl")

# PCA plot
dds_var_transform = varianceStabilizingTransformation(deseq_data)
plotPCA(dds_var_transform , intgroup = 'toxin')
