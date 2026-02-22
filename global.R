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

clonograph<-function (sce, complete_only = FALSE, color_pal = "Reds", QC_stats = FALSE) 
{
  consolidated_clonal_abundance <- sce@metadata$Clones %>% 
    dplyr::group_by(Clone, Group) %>% dplyr::mutate(AF_med = mean(AF_med), 
                                                    DP_med = mean(DP_med), GQ_med = mean(GQ_med)) %>% dplyr::select(!c(variants, 
                                                                                                                       Count)) %>% dplyr::distinct() %>% dplyr::rowwise() %>% 
    dplyr::mutate(Count = ifelse(Group == "Other", n_Other, 
                                 n_Complete)) %>% dplyr::mutate(Total_count = n_Other + 
                                                                  n_Complete) %>% dplyr::ungroup() %>% dplyr::arrange(.data$Total_count)
  if (complete_only == TRUE) {
    consolidated_clonal_abundance <- consolidated_clonal_abundance %>% 
      dplyr::filter(.data$Group == "Complete") %>% dplyr::arrange(.data$Count) %>% 
      dplyr::select(-LCI, -UCI) %>% dplyr::rename(LCI = Complete_LCI, 
                                                  UCI = Complete_UCI)
  }
  clonal_architecture <- sce@metadata$Architecture
  mutant_order <- setdiff(colnames(sce@metadata$NGT), c("Cell", 
                                                        "Clone", "Group"))
  clonal_architecture$Clone <- factor(clonal_architecture$Clone, 
                                      levels = unique(rev(consolidated_clonal_abundance$Clone)))
  consolidated_clonal_abundance$Clone <- factor(consolidated_clonal_abundance$Clone, 
                                                levels = levels(clonal_architecture$Clone))
  consolidated_clonal_abundance$Group <- factor(consolidated_clonal_abundance$Group, 
                                                levels = c("Complete", "Other"))
  gg_clonal_barplot <- ggplot(data = consolidated_clonal_abundance, 
                              aes(x = Clone, y = Count, fill = Group)) + geom_bar(fun = "identity", 
                                                                                  stat = "summary", position = "stack") + theme_classic(base_size = 7) + 
    scale_y_continuous(expand = c(0.01, 0)) + ylab("Cell Count") + 
    geom_errorbar(aes(ymin = LCI, ymax = UCI), width = 0.2) + 
    scale_fill_manual(values = c(Other = "Grey70", Complete = RColorBrewer::brewer.pal(n = 5, 
                                                                                       name = color_pal)[5])) + theme(axis.title.x = element_blank(), 
                                                                                                                      axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                                                                                                                      axis.line.x = element_blank(), legend.position = "right", 
                                                                                                                      plot.margin = unit(c(0, 0, 0, 0), "cm"))
  gg_heatmap <- ggplot(data = clonal_architecture, aes(x = Clone, 
                                                       y = final_annot, fill = Genotype)) + geom_tile() + scale_fill_manual(values = c(WT = RColorBrewer::brewer.pal(7, 
                                                                                                                                                                     color_pal)[1], Heterozygous = RColorBrewer::brewer.pal(7, 
                                                                                                                                                                                                                            color_pal)[3], Homozygous = RColorBrewer::brewer.pal(7, 
                                                                                                                                                                                                                                                                                 color_pal)[6], Unknown = "grey50"), name = "Genotype") + 
    theme_classic(base_size = 7) + ylab("Mutation") + scale_y_discrete(limits = rev(mutant_order)) + 
    theme(legend.position = "right", legend.direction = "vertical", 
          axis.text.x = element_blank(), axis.line = element_blank(), 
          axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
          plot.margin = unit(c(0, 0, 0, 0), "cm"))
  if (QC_stats == FALSE) {
    return(cowplot::plot_grid(gg_clonal_barplot, gg_heatmap, 
                              ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 
                                                                                  0.5)))
  }
  else {
    gg_QC_heatmap_GQ <- ggplot(data = consolidated_clonal_abundance, 
                               aes(x = Clone, y = Group, fill = GQ_med)) + geom_tile() + 
      colorspace::scale_fill_continuous_divergingx(palette = "RdBu", 
                                                   mid = 30, rev = TRUE, na.value = "grey80") + 
      theme_classic(base_size = 7) + theme(legend.position = "right", 
                                           legend.direction = "horizontal", axis.text.x = element_blank(), 
                                           axis.line = element_blank(), axis.title.x = element_blank(), 
                                           axis.ticks.x = element_blank(), plot.margin = unit(c(0, 
                                                                                                0, 0, 0), "cm"))
    gg_QC_heatmap_DP <- ggplot(data = consolidated_clonal_abundance, 
                               aes(x = Clone, y = Group, fill = DP_med)) + geom_tile() + 
      colorspace::scale_fill_continuous_divergingx(palette = "RdBu", 
                                                   mid = 10, rev = TRUE, na.value = "grey80") + 
      theme_classic(base_size = 7) + theme(legend.position = "right", 
                                           legend.direction = "horizontal", axis.text.x = element_blank(), 
                                           axis.line = element_blank(), axis.title.x = element_blank(), 
                                           axis.ticks.x = element_blank(), plot.margin = unit(c(0, 
                                                                                                0, 0, 0), "cm"))
    return(cowplot::plot_grid(gg_clonal_barplot, gg_heatmap, 
                              gg_QC_heatmap_GQ, gg_QC_heatmap_DP, ncol = 1, align = "v", 
                              axis = "lr", rel_heights = c(1, 0.5, 0.25, 0.25, 
                                                           0.25)))
  }
}

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

