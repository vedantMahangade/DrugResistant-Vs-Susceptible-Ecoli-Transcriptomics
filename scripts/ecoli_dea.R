library(ggplot2)
library(pheatmap)
library(tidyverse)


# Set working directory
setwd("~/Projects/DrugResistant-Vs-Susceptible-Ecoli-Transcriptomics")

# Load expression matrix
# GEO file includes metadata prefixed with !. 
# comment.char = "!" skips those lines
geo_data <- read.delim(
  "data/GSE59408_series_matrix.txt", 
  comment.char = "!", 
  header = TRUE, 
  row.names = 1, 
  check.names = FALSE
)

head(geo_data)
str(geo_data)

# Converting to matrix
expr_matrix <- as.matrix(geo_data)
# rownames(expr_matrix) <- geo_data$ID_REF

# Create Sample Metadata
sample_names <- c(
  "Parent_strain_rep1", "Parent_strain_rep2",
  "Cefoperazone_CPZ_Line1", "Cefoperazone_CPZ_Line2", "Cefoperazone_CPZ_Line3", "Cefoperazone_CPZ_Line4",
  "Cefixime_CFIX_Line1", "Cefixime_CFIX_Line2", "Cefixime_CFIX_Line3", "Cefixime_CFIX_Line4",
  "Amikacin_AMK_Line1", "Amikacin_AMK_Line2", "Amikacin_AMK_Line3", "Amikacin_AMK_Line4",
  "Neomycin_NM_Line1", "Neomycin_NM_Line2", "Neomycin_NM_Line3", "Neomycin_NM_Line4",
  "Doxycycline_DOXY_Line1", "Doxycycline_DOXY_Line2", "Doxycycline_DOXY_Line3", "Doxycycline_DOXY_Line4",
  "Chloramphenicol_CP_Line1", "Chloramphenicol_CP_Line2", "Chloramphenicol_CP_Line3", "Chloramphenicol_CP_Line4",
  "Azithromycin_AZM_Line1", "Azithromycin_AZM_Line2", "Azithromycin_AZM_Line3", "Azithromycin_AZM_Line4",
  "Trimethoprim_TP_Line1", "Trimethoprim_TP_Line2", "Trimethoprim_TP_Line3", "Trimethoprim_TP_Line4",
  "Enoxacin_ENX_Line1", "Enoxacin_ENX_Line2", "Enoxacin_ENX_Line3", "Enoxacin_ENX_Line4",
  "Ciprofloxacin_CPFX_Line1", "Ciprofloxacin_CPFX_Line2", "Ciprofloxacin_CPFX_Line3", "Ciprofloxacin_CPFX_Line4"
)


# Creating sample metadata
metadata <- data.frame(
  # SampleID = colnames(expr_matrix),
  SampleName = sample_names,
  IsResistant = c("No", "No", rep("Yes", 40)),
  ResistanceType = c(
    "No", "No", 
    rep("CPZ", 4), rep("CFIX", 4), rep("AMK", 4), rep("NM", 4), rep("DOXY", 4),
    rep("CP", 4), rep("AZM", 4), rep("TP", 4), rep("ENX", 4), rep("CPFX", 4)
  )
)

# Alignment of metadata and matrix
rownames(metadata) <- colnames(expr_matrix)

# # Extract metadata manually by reading the file as plain text
# file_metadata <- readLines('data/GSE59408_series_matrix.txt')
# file_metadata <- file_metadata[grepl("^!", file_metadata)]
# head(file_metadata) # View the metadata lines

# Preprocessing Function
preprocess_data <- function(expr_matrix, metadata, min_expression = 300, min_sd = 0.5) {
  
  # Remove rows with missing values
  expr_matrix <- na.omit(expr_matrix)  
  
  # Filter low-expression genes
  expr_matrix <- expr_matrix[rowSums(expr_matrix >= min_expression) > 0, ] 
  
  # Log2 transformation
  expr_matrix <- log2(expr_matrix + 1)
  
  # Calculate standard deviation
  gene_sd <- apply(expr_matrix, 1, sd)
  
  # Filter low-variability genes
  expr_matrix <- expr_matrix[gene_sd >= min_sd, ] 
  
  return(expr_matrix)
}


# PCA Function
perform_pca <- function(expr_matrix, metadata) {
  pca <- prcomp(t(expr_matrix), scale. = TRUE)
  explained_variance <- round((pca$sdev^2 / sum(pca$sdev^2)) * 100, 1)
  pca_data <- as.data.frame(pca$x)
  pca_data$IsResistant <- metadata$IsResistant
  pca_data$ResistanceType <- metadata$ResistanceType
  list(pca_data = pca_data, explained_variance = explained_variance)
}


# Volcano Plot Function
create_volcano_plot <- function(tweedie_results, title_prefix, path) {
  
  # Prepare the data for the volcano plot
  tweedie_results$logFC <- tweedie_results$coef
  tweedie_results$negLog10qval <- -log10(tweedie_results$qval)
  
  volcano_plot <- EnhancedVolcano(tweedie_results,
                                  lab = tweedie_results$feature,
                                  title = paste0(title_prefix,' vs No resistance'),
                                  x = 'logFC',
                                  y = 'pval',
                                  xlim = c(-1.5, 1.5),
                                  ylim = c(0, 7),
                                  pCutoff = 0.25,
                                  FCcutoff = 0.25,
                                  pointSize = 3.0,
                                  labSize = 3.0)
  
  # Save the volcano plot
  ggsave(
    filename = paste0(path, "/volcano_plot.png"),
    plot = volcano_plot,
    width = 8, height = 6, dpi = 300
  )
}


# Heatmap Function
create_heatmap <- function(expr_matrix, significant_genes, title_prefix, path) {
  sig_expr_matrix <- expr_matrix[rownames(expr_matrix) %in% significant_genes, ]
  # Save heatmap as a PNG
  png(filename = paste0(path, "/heatmap.png"), width = 800, height = 600)
  pheatmap(
    sig_expr_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    display_numbers = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(50),
    main = paste0(title_prefix, " Heatmap of Differentially Expressed Genes"),
    fontsize = 8
  )
  dev.off()  # Close the device
}



# Antibiotics list
antibiotics <- c("CPZ", "CFIX", "AMK", "NM", "DOXY", "CP", "AZM", "TP", "ENX", "CPFX")
i<-1

results_list <- list()
significant_genes <- list()

# Loop through each antibiotic
for (ab in antibiotics) {
  
  # Subset data for the current antibiotic
  ab_data <- expr_matrix[, metadata$ResistanceType %in% c('No', ab)]
  ab_metadata <- metadata[metadata$ResistanceType %in% c('No', ab), ]
  
  path = paste0("analysis/", as.character(i), "_", ab)
  
  dir.create(path)
  
  # Preprocess the data
  ab_data <- preprocess_data(ab_data, ab_metadata)
  
  # PCA Analysis
  pca_results <- perform_pca(ab_data, ab_metadata)
  pca_data <- pca_results$pca_data
  explained_variance <- pca_results$explained_variance
  
  # Plot PCA Results
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = IsResistant)) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      title = paste0("PCA of Gene Expression Data (", ab, " Resistance)"),
      x = paste0("PC1 (", explained_variance[1], "% variance)"),
      y = paste0("PC2 (", explained_variance[2], "% variance)"),
      color = "Resistance"
    )
  # save plot
  ggsave(
    filename = paste0(path,"/PCA", ".png"),
    plot = pca_plot,
    width = 8, height = 6, dpi = 300
  )
  
  # Run Tweedieverse
  Tweedieverse(
    input_features = as.data.frame(ab_data),
    input_metadata = metadata,
    output = paste0(path, "/tweedieverse/"),
    base_model = "CPLM",
    fixed_effects = 'ResistanceType',
    max_significance = 0.05,
    correction = "BH",
    plot_heatmap = TRUE,
    plot_scatter = TRUE,
    heatmap_first_n = 50
  )
  
  tweedie_results <- read.delim(paste0(path,"/tweedieverse/all_results.tsv"), 
                                header = TRUE,
                                sep = '\t')
  
  # Volcano Plot
  create_volcano_plot(tweedie_results, ab, path)
  
  # Heatmap for Significant Genes
  create_heatmap(as.data.frame(ab_data),
                 tweedie_results$feature[tweedie_results$qval < 0.05], ab, path)
  
  # Filter significant genes (q-value < 0.05)
  sig_genes <- tweedie_results$feature[tweedie_results$qval < 0.05]
  significant_genes[[ab]] <- sig_genes
  
  # Store full results for meta-analysis
  tweedie_results$Antibiotic <- ab
  results_list[[ab]] <- tweedie_results
  
  i <- i+1
}

# Merge results into a single data frame and save
combined_results <- bind_rows(results_list)


top_genes <- combined_results %>%
  filter(qval < 0.05) %>%
  group_by(Antibiotic) %>%
  arrange(qval) %>%
  slice_head(n = 20) %>%
  select(Antibiotic, feature, coef, pval, qval)

write.csv(combined_results, file.path('analysis/summary', "aggregated_results.csv"), row.names = FALSE)
write.csv(top_genes, file.path('analysis/summary', "top_genes.csv"), row.names = FALSE)


# Get shared and unique gene sets
shared_genes <- Reduce(intersect, significant_genes)
unique_genes <- lapply(names(significant_genes), function(ab) {
  setdiff(significant_genes[[ab]], unlist(significant_genes[names(significant_genes) != ab]))
})
# Assign antibiotic names to the list
names(unique_genes) <- antibiotics  


# Create a Heatmap of Top Genes
heatmap_data <- combined_results %>%
  select(feature, Antibiotic, coef) %>%
  pivot_wider(names_from = Antibiotic, values_from = coef, values_fill = 0)

heatmap_matrix <- as.matrix(heatmap_data[, -1])  # Remove 'feature' column
rownames(heatmap_matrix) <- heatmap_data$feature

# Generate heatmap
png(file.path('analysis/summary', "heatmap_top_genes.png"), width = 800, height = 600)
pheatmap(
  heatmap_matrix,
  # cluster_rows = TRUE,
  # cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of Coefficients for Top Differentially Expressed Genes",
  fontsize = 8
)
dev.off()


library(VennDiagram)
library(UpSetR)
top_genes <- read.delim('analysis/summary/top_genes.csv', sep = ',', header = TRUE,check.names = FALSE)

gene_sets <- lapply(unique(top_genes$Antibiotic), function(ab) {
  top_genes$aeature[top_genes$Antibiotic == ab]
})
names(gene_sets) <- unique(top_genes$Antibiotic)

# Create a Venn diagram
venn.diagram(
  x = gene_sets,
  category.names = names(gene_sets),
  filename = "venn_diagram.png"
)



# Create a binary membership matrix
gene_list <- unique(unlist(significant_genes))  # All unique significant genes
binary_matrix <- sapply(significant_genes, function(gene_set) {
  gene_list %in% gene_set
})
rownames(binary_matrix) <- gene_list
colnames(binary_matrix) <- antibiotics

# Convert binary matrix to a data frame for ComplexUpset
binary_df <- as.data.frame(binary_matrix)
binary_df$Gene <- rownames(binary_matrix)

# Melt the data frame for compatibility with ComplexUpset
binary_long <- binary_df %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Antibiotic",
    values_to = "Membership"
  ) %>%
  filter(Membership == TRUE) %>%
  pivot_wider(
    names_from = Antibiotic,
    values_from = Membership,
    values_fill = FALSE
  )

library(ComplexUpset)
# Plot the UpSet plot
upset_plot <- upset(
  binary_long,
  antibiotics,
  annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE
    )
  ),
  min_size = 2,  # Minimum size of intersections to display
  width_ratio = 0.2,  # Adjust ratios if needed
  height_ratio = 0.8
)

# Save the plot
ggsave(
  filename = "analysis/summary/upset_plot.png",
  plot = upset_plot,
  width = 10, height = 6, dpi = 300
)

