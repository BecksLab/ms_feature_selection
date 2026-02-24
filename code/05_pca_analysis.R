# ================================================
# 05_pca_analysis_MODIFIED.R
# ================================================

library(FactoMineR)
library(factoextra)
library(tidyverse)
library(ggrepel)

# Load datasets
predictors <- read.csv("data/cleaned/all_networks.csv") %>%
  as_tibble() %>%
  vibe_check(-c(ρ, complexity, robustness))

# 2. Run the PCA on Predictors only
pca_res <- PCA(predictors, scale.unit = TRUE, ncp = 10, graph = FALSE)

# 3. DIAGNOSTIC: How many dimensions matter?
# Look for the 'elbow' or Eigenvalues > 1
fviz_eig(pca_res, addlabels = TRUE)

# Check the eigenvalues
eig.val <- get_eigenvalue(pca_res)
print(eig.val)

# 4. Get the contribution matrix for all metrics across all 6 dimensions
contrib_matrix <- var_info$contrib[, 1:6]

# Function to find the top metric for each dimension, 
# but ensuring we don't pick the same metric twice
get_unique_representatives <- function(contrib_mat) {
  reps <- character(ncol(contrib_mat))
  used_metrics <- c()
  
  for (i in 1:ncol(contrib_mat)) {
    # Sort metrics for this dimension by contribution descending
    sorted_metrics <- sort(contrib_mat[, i], decreasing = TRUE)
    
    # Pick the top one that isn't already in our 'used' list
    candidate <- names(sorted_metrics)[!(names(sorted_metrics) %in% used_metrics)][1]
    
    reps[i] <- candidate
    used_metrics <- c(used_metrics, candidate)
  }
  return(reps)
}

# Apply the function
unique_reps <- get_unique_representatives(contrib_matrix)

# Create your final Table
representative_table <- data.frame(
  Dimension = paste("Dim", 1:6),
  Eigenvalue = eig.val[1:6, 1],
  Variance_Explained = eig.val[1:6, 2],
  Representative_Metric = unique_reps
)

print(representative_table)

# 5. VISUALIZE: The "Causal Hierarchy" Biplot
# Map your categories back onto the arrows
pca_arrows <- as.data.frame(var_stats$coord) %>%
  rownames_to_column("Property") %>%
  mutate(Category = case_when(
    Property %in% geometry_metrics ~ "Geometry",
    Property %in% path_metrics ~ "Path",
    Property %in% node_metrics ~ "Node",
    TRUE ~ "Other"
  ))

ggplot(pca_arrows, aes(x=Dim.1, y=Dim.2)) +
  geom_hline(yintercept=0, linetype=2, color="grey70") +
  geom_vline(xintercept=0, linetype=2, color="grey70") +
  geom_segment(aes(x=0, y=0, xend=Dim.1, yend=Dim.2, color=Category),
               arrow = arrow(length=unit(0.1,"inches")), linewidth = 1) +
  geom_text_repel(aes(label=Property, color=Category), max.overlaps = 15) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Network Structural Space",
       subtitle = "Dimensions colored by Energy-Flow Hierarchy",
       x = paste0("PC1 (", round(pca_res$eig[1,2], 1), "%)"),
       y = paste0("PC2 (", round(pca_res$eig[2,2], 1), "%)")) +
  theme_minimal()

ggsave("../figures/pca.png", combined_plot,
       width=4000, height=6500, units="px", dpi=600)
