library(dplyr)
library(arrow)
library(patchwork)

# -------------------------
# Cohort CSVs
# -------------------------
p1 <- "data/pseudo_bulk_consolidated_clone.csv"
p2 <- "data/pseudo_bulk_exact_clone.csv"
p3 <- "data/pseudo_bulk_consolidated_clone_statistics.csv"

if (!file.exists(p1)) p1 <- "pseudo_bulk_consolidated_clone.csv"
if (!file.exists(p2)) p2 <- "pseudo_bulk_exact_clone.csv"
if (!file.exists(p3)) p3 <- "pseudo_bulk_consolidated_clone_statistics.csv"

pseudo_bulk_consolidated_clone <- read.csv(p1, stringsAsFactors = FALSE)
pseudo_bulk_exact_clone        <- read.csv(p2, stringsAsFactors = FALSE)

# This is the PAIRWISE stats table (Protein, p_val, avg_log2FC, p_val_adj, clone1, clone2)
pseudo_bulk_consolidated_clone_statistics <- read.csv(p3, stringsAsFactors = FALSE)

# -------------------------
# Single-sample Parquet
# -------------------------
cm <- "data/cell_meta.parquet"
ex <- "data/expr_long.parquet"
if (!file.exists(cm)) cm <- "cell_meta.parquet"
if (!file.exists(ex)) ex <- "expr_long.parquet"

cell_meta <- arrow::read_parquet(cm)
expr_long <- arrow::read_parquet(ex)

# Standardize names if needed
if ("UMAP_X" %in% names(cell_meta) && !"UMAP_1" %in% names(cell_meta)) {
  cell_meta <- cell_meta %>% rename(UMAP_1 = UMAP_X)
}
if ("UMAP_Y" %in% names(cell_meta) && !"UMAP_2" %in% names(cell_meta)) {
  cell_meta <- cell_meta %>% rename(UMAP_2 = UMAP_Y)
}
if ("Stage_code" %in% names(cell_meta) && !"Stage" %in% names(cell_meta)) {
  cell_meta <- cell_meta %>% rename(Stage = Stage_code)
}

# Ensure key columns are character
cell_meta <- cell_meta %>% mutate(Cell = as.character(Cell))
expr_long <- expr_long %>% mutate(Cell = as.character(Cell),
                                  Protein = as.character(Protein))


stats_file <- "data/sample_statistics.csv"
if (!file.exists(stats_file)) stats_file <- "sample_statistics.csv"

sample_stats <- read.csv(stats_file, stringsAsFactors = FALSE) %>%
  dplyr::mutate(
    Stage = factor(Stage, levels = c("Diagnosis","CR","Relapse")),
    Epi_group = factor(Epi_group, levels = c("None","DNMT3A","TET2","IDH","Combo")),
    Sig_group = factor(Sig_group, levels = c("None","RAS","FLT3","R+F"))
  )

sample_metrics <- c(
  "Number_of_clones", "Number_of_mutations",
  "Number_of_mutations_in_dominant_clone",
  "Dominant_clone_size", "Shannon"
)

sample_groups <- c("Epi_group", "Sig_group", "Stage")

color_scheme <- list(
  Epi_group = c("None"="grey65","DNMT3A"="#FDA44B","TET2"="#85010C","IDH"="#1A62A4","Combo"="#CCA4FC"),
  Sig_group = c("None"="grey65","RAS"="#EA5B00","FLT3"="#298C1F","R+F"="#562887"),
  Stage = c("Diagnosis"="#E4C04C","CR"="#9C3130","Relapse"="#A8C2DA")
)