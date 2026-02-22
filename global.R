# global.R
library(shiny)
library(bs4Dash)
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(Matrix)
library(RColorBrewer)
library(SingleCellExperiment)
library(scDNA, include.only = c("clonograph"))


make_clone_cluster_palette <- function(clusters) {
  clusters <- as.character(clusters)
  n <- length(clusters)
  # Brewer max is 12; interpolate if needed
  base_cols <- pals::kelly(min(12, n)+2)[-c(1,2)]
  pal <- colorRampPalette(base_cols)(n)
  names(pal) <- clusters
  pal
}

make_cell_cluster_palette <- function(clusters) {
  clusters <- as.character(clusters)
  n <- length(clusters)
  # Brewer max is 12; interpolate if needed
  base_cols <- pals::tol(min(12, n)+1)
  pal <- colorRampPalette(base_cols)(n)
  names(pal) <- clusters
  pal
}

# ---- Load the cohort Seurat object ----
seurat_obj <- readRDS("./Data/global_seurat.rds")
stopifnot(inherits(seurat_obj, "Seurat"))

# ---- Standard column / reduction naming (adjustable) ----
SAMPLE_COL <- if ("Sample" %in% colnames(seurat_obj@meta.data)) "Sample" else "orig.ident"
UMAP_REDUCTION <- if ("umap.harmony" %in% names(seurat_obj@reductions)) "umap.harmony" else "umap"

# Which column to use as "cluster" when coloring by cluster
CLUSTER_COL <- if ("harmony_clusters" %in% colnames(seurat_obj@meta.data)) "harmony_clusters" else "seurat_clusters"

# ---- Helpers ----
get_samples <- function(obj, sample_col = SAMPLE_COL) {
  sort(unique(obj@meta.data[[sample_col]]))
}

get_proteins <- function(obj) {
  # Your object has Protein assay active; keep this general anyway
  assay <- if ("Protein" %in% Assays(obj)) "Protein" else Assays(obj)[1]
  rownames(GetAssayData(obj, assay = assay, slot = "data"))
}

subset_to_sample <- function(obj, sample_value, sample_col = SAMPLE_COL) {
  md <- obj@meta.data
  keep_cells <- rownames(md)[md[[sample_col]] == sample_value]
  subset(obj, cells = keep_cells)
}

safe_umap_plot <- function(obj, group_by = NULL, reduction = UMAP_REDUCTION) {
  if (!reduction %in% names(obj@reductions)) {
    return(ggplot() + theme_void() + ggtitle(paste0("Reduction not found: ", reduction)))
  }
  tryCatch({
    if (!is.null(group_by) && group_by %in% colnames(obj@meta.data)) {
      DimPlot(obj, reduction = reduction, group.by = group_by) + NoLegend()
    } else {
      DimPlot(obj, reduction = reduction) + NoLegend()
    }
  }, error = function(e) ggplot() + theme_void() + ggtitle(e$message))
}

safe_protein_umap <- function(obj, protein, reduction = UMAP_REDUCTION) {
  if (!reduction %in% names(obj@reductions)) {
    return(ggplot() + theme_void() + ggtitle(paste0("Reduction not found: ", reduction)))
  }
  if (!"Protein" %in% Assays(obj)) {
    return(ggplot() + theme_void() + ggtitle("Protein assay not found"))
  }
  tryCatch({
    DefaultAssay(obj) <- "Protein"
    FeaturePlot(obj, reduction = reduction, features = protein)
  }, error = function(e) ggplot() + theme_void() + ggtitle(e$message))
}


# ---- Load SCEs (choose ONE pattern) ----
# Pattern A: list
sce_list <- readRDS("./Data/DP45_diet_sce_list.rds")

has_sce_list <- exists("sce_list") && is.list(sce_list) && length(sce_list) > 0
has_sce_obj  <- exists("sce_obj")  && inherits(sce_obj, "SingleCellExperiment")

get_sce_for_sample <- function(sample_id) {
  if (has_sce_list) {
    if (sample_id %in% names(sce_list)) return(sce_list[[sample_id]])
    return(NULL)
  }
  if (has_sce_obj) {
    cd <- as.data.frame(SummarizedExperiment::colData(sce_obj))
    if (!("Sample" %in% colnames(cd))) return(NULL)
    keep <- rownames(cd)[cd$Sample == sample_id]
    if (length(keep) == 0) return(NULL)
    return(sce_obj[, keep, drop = FALSE])
  }
  NULL
}


summarized_out<-readRDS("./Data/DP45_summarized_out_2025_12_20.rds")

