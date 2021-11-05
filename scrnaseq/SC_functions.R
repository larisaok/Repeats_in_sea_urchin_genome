# Slightly modified from BUSpaRse, just to avoid installing a few dependencies not used here
read_count_output <- function(filename) {
  m <- readMM(filename)
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  genes_file <- paste0(tools::file_path_sans_ext(filename), ".genes.txt")
  barcodes_file <- paste0(tools::file_path_sans_ext(filename), ".barcodes.txt")
  genes <- readLines(file(genes_file))
  barcodes <- readLines(file(barcodes_file))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}


knee_plot_1 <- function(dgCMatrix) {
  bc_rank <- barcodeRanks(dgCMatrix, lower = 10)
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]],
                     inflection = metadata(bc_rank)[["inflection"]],
                     rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]])) %>%
    distinct() %>%
    dplyr::filter(total > 0)
  return(knee_plt)
}

knee_plot <- Vectorize(knee_plot_1, SIMPLIFY = F)


pca <- function(dgCMatrix){
  tmp <- prcomp_irlba(dgCMatrix, n = 2) # scales and centers by default
  tmp_df <- as.data.frame(tmp$x)
  return(tmp_df)
}
qc_pca <- Vectorize(pca, SIMPLIFY = F)


to_tibble <- function(dgCMatrix){
  tibl <- tibble(nCount = rowSums(dgCMatrix),
                nGene = rowSums(dgCMatrix > 0),
                id = deparse(substitute(dgCMatrix)))
  return(tibl)
}
qc_gene_num <- Vectorize(to_tibble, SIMPLIFY = F)

knee_filter <- function(dgCMatrix){
  bc_rank <- barcodeRanks(dgCMatrix, lower = 10)
  tot_counts <- Matrix::colSums(dgCMatrix)
  dgCMatrix_filtered <- dgCMatrix[, tot_counts > metadata(bc_rank)$inflection]
}


make_plots <- function(seu_object, path, top10_variable, good_markers){
  #Dotplots
  dp <- DotPlot(seu_object, features = top10_variable, cols = c("blue", "red")) +
  coord_flip()
  ggsave(file = paste(path, "top10_dot.png", sep = ""), plot = dp, dpi = 320)

  dp <- DotPlot(seu_object, features = unique(good_markers), cols = c("blue", "red")) +
  coord_flip()
  ggsave(file = paste(path, "good_marks_dot.png", sep = ""), plot = dp, dpi = 320)


  #Feature plots
  fp <- FeaturePlot(seu_object, top10_variable, cols = c("blue", "red"),
                    reduction="umap", order = T, ncol = 5, slot = "scale.data",
                    keep.scale = "all", repel = T,  label = T, 
                    interactive = F, combine = F, coord.fixed = T)
  for (i in fp){
    plot_name <- paste(path, "fp_top10_", i$labels$title, ".png", sep = "")
    ggsave(filename = plot_name, plot = i, device = "png", dpi = 320)
  } 
  
  
  fp <- FeaturePlot(seu_object, good_markers, cols = c("blue", "red"),
                    reduction="umap", order = T, slot = "scale.data", 
                    keep.scale = "all", repel = T,  label = T, 
                    interactive = F, combine = F, coord.fixed = T)
  
  for (i in fp){
    plot_name <- paste(path, "fp_good_marks_", i$labels$title, ".png", sep = "")
    ggsave(filename = plot_name, plot = i, device = "png", dpi = 320)
  } 

  #Violin plots

  vp <- VlnPlot(seu_object, slot = "scale.data", good_markers, pt.size = 0.01, combine = F)
  for (i in vp){
    plot_name <- paste(path, "vp_good_markers_", i$labels$title, ".png", sep = "")
    ggsave(filename = plot_name, plot = i, device = "png", dpi = 320)
  } 
  
  vp <- VlnPlot(seu_morula, slot = "scale.data", top10_variable, pt.size = 0.01, combine = F)
  vp
  for (i in vp){
    plot_name <- paste(path, "vp_top10_", i$labels$title, ".png", sep = "")
    ggsave(filename = plot_name, plot = i, device = "png", dpi = 320)
  } 
  
  VFP <- VariableFeaturePlot(seu_object, log = FALSE)
  VFP <- LabelPoints(plot = VFP, points = top10_variable, repel = TRUE)
  ggsave(paste(path, "VFP.png", sep = ""), plot = VFP, dpi = 320)
  #return (TRUE)
}

