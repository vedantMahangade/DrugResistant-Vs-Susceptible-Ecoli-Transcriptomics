library(ggplot2)
library(pheatmap)
library(dplyr)
library(readr)
library(stringr)
library(gridExtra)
library(knitr)

output_dir <- "analysis/"
summary_dir <- "summary/"
dir.create(summary_dir, showWarnings = FALSE)

# Aggregate Results from all_results.tsv
aggregate_results <- function(output_dir) {
  result_list <- list()
  for (folder in list.dirs(output_dir, recursive = FALSE)) {
    ab_name <- unlist(str_split(basename(folder),'_'))[2]
    result_path <- file.path(folder, "tweedieverse", "all_results.tsv")
    if (file.exists(result_path)) {
      ab_results <- read.delim(result_path)
      ab_results$Antibiotic <- ab_name
      result_list[[ab_name]] <- ab_results
    }
  }
  combined_results <- bind_rows(result_list)
  write.csv(combined_results, file.path(summary_dir, "combined_results.csv"), row.names = FALSE)
  return(combined_results)
}

combined_results <- aggregate_results(output_dir)


# Extract Top Differentially Expressed Genes
top_genes <- combined_results %>%
  filter(qval < 0.05) %>%
  group_by(Antibiotic) %>%
  arrange(qval) %>%
  slice_head(n = 20) %>%
  select(Antibiotic, feature, coef, pval, qval)

write.csv(top_genes, file.path(summary_dir, "top_genes.csv"), row.names = FALSE)


# Create a Heatmap of Top Genes
heatmap_data <- combined_results %>%
  select(feature, Antibiotic, coef) %>%
  pivot_wider(names_from = Antibiotic, values_from = coef, values_fill = 0)

heatmap_matrix <- as.matrix(heatmap_data[, -1])  # Remove 'feature' column
rownames(heatmap_matrix) <- heatmap_data$feature

# Generate heatmap
png(file.path(summary_dir, "heatmap_top_genes.png"), width = 800, height = 600)
pheatmap(
  heatmap_matrix,
  # cluster_rows = TRUE,
  # cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of Coefficients for Top Differentially Expressed Genes",
  fontsize = 8
)
dev.off()


volcano_summary <- function(output_dir, summary_dir) {
  # Initialize a list to store paths to volcano plots
  volcano_images <- list()
  
  # Loop through each antibiotic folder to find the volcano plot
  for (folder in list.dirs(output_dir, recursive = FALSE)) {
    ab_name <- basename(folder)
    volcano_path <- file.path(folder, "volcano_plot.png")
    if (file.exists(volcano_path)) {
      volcano_images[[ab_name]] <- volcano_path
    }
  }
  
  # Combine images into a single PDF with one image per page
  pdf(file.path(summary_dir, "volcano_summary.pdf"), width = 8, height = 6)
  for (ab_name in names(volcano_images)) {
    # Start a new page
    grid::grid.newpage()
    
    # Read the volcano plot image
    img <- png::readPNG(volcano_images[[ab_name]])
    
    # Add the image to the current page
    grid::grid.raster(img, interpolate = FALSE)
    
    # Add a title for context
    # grid::grid.text(paste("Volcano Plot -", ab_name), y = 0.95, gp = grid::gpar(fontsize = 14, fontface = "bold"))
  }
  dev.off()
}

# Run the function
volcano_summary(output_dir, summary_dir)



# Generate Summary Report
generate_report <- function(summary_dir, top_genes, heatmap_path) {
  report_file <- file.path(summary_dir, "analysis_summary.md")
  sink(report_file)
  cat("# Analysis Summary Report\n")
  cat("\n## Top Differentially Expressed Genes\n")
  print(kable(top_genes))
  cat("\n## Heatmap of Top Genes\n")
  cat("![Heatmap](heatmap_top_genes.png)\n")
  cat("\n## Volcano Plot Summary\n")
  cat("See volcano_summary.pdf for all volcano plots.\n")
  sink()
}

generate_report(summary_dir, top_genes, "heatmap_top_genes.png")
