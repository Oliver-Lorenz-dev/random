library(edgeR)
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


Toxin <- factor(c("Y","Y","N","N"))
data.frame(Sample=colnames(y),Toxin)
design <- model.matrix(~Toxin)
rownames(design) <- colnames(y)

y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)

fit <- glmFit(y, design)
lrt <- glmLRT(fit)
topTags(lrt)

plotMD(lrt)
abline(h=c(-1, 1), col="blue")


# pathway analysis
expression_table = lrt[["table"]]
pathway_df = expression_table[ -c( 2,3) ]
pathway_df <- tibble::rownames_to_column(pathway_df, "Gene.symbol")
colnames(pathway_df)[colnames(pathway_df) == 'PValue'] <- 'adj.P.val'

library(pathfindR)

# process inputs
RA_processed <- input_processing(input = pathway_df, # results from edgeR as input (in form of dataframe)
                                 p_val_threshold = 0.05, # p value threshold to filter significant genes
                                 pin_name_path  = "Biogrid", # the name of the PIN to use for active subnetwork search
                                 convert2alias = TRUE) # boolean indicating whether or not to convert missing symbols to alias symbols in the PIN

# The available gene sets in pathfindR are “KEGG”, “Reactome”, “BioCarta”, “GO-All”, “GO-BP”, “GO-CC” and “GO-MF”. 
biocarta_list <- fetch_gene_set(gene_sets = "BioCarta",
                                min_gset_size = 10,
                                max_gset_size = 300)
biocarta_gsets <- biocarta_list[[1]]
biocarta_descriptions <- biocarta_list[[2]]

# number of iterations
n_iter <- 10

# to store the result of each iteration
combined_res <- NULL

for (i in 1:n_iter) {
  
  # Active Subnetwork Search
  snws_file <- paste0("active_snws_", i) 
  active_snws <- active_snw_search(input_for_search = RA_processed, 
                                   pin_name_path = "Biogrid", 
                                   snws_file = snws_file,
                                   score_quan_thr = 0.8, # arg can be changed
                                   sig_gene_thr = 0.02, # arg can be changed
                                   search_method = "GR")
  
  # Enrichment Analyses
  current_res <- enrichment_analyses(snws = active_snws,
                                     sig_genes_vec = RA_processed$GENE,
                                     pin_name_path = "Biogrid", 
                                     genes_by_term = biocarta_gsets,
                                     term_descriptions = biocarta_descriptions,
                                     adj_method = "bonferroni",
                                     enrichment_threshold = 0.05,
                                     list_active_snw_genes = TRUE) 
  
  # Combine results 
  combined_res <- rbind(combined_res, current_res)
}

# Summarize Combined Enrichment Results
summarized_df <- summarize_enrichment_results(combined_res, 
                                              list_active_snw_genes = TRUE)

# Annotate Affected Genes Involved in Each Enriched Term
final_res <- annotate_term_genes(result_df = summarized_df, 
                                 input_processed = RA_processed, 
                                 genes_by_term = biocarta_gsets)

visualize_terms(result_df = final_res, 
                hsa_KEGG = FALSE, # true if using human KEGG gene set
                pin_name_path = "Biogrid")

enrichment_chart(final_res[1:10, ])

## biocarta is the quickest pathway analysis to run as it is the
## smallest gene set, this does take a while to run
