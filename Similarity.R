# Script for Correlation based similarity analysis

# Load packages
suppressPackageStartupMessages(library(scrabbitr))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(miloR))
suppressPackageStartupMessages(library(DelayedArray))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ggrastr))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(ggalluvial))
suppressPackageStartupMessages(library(ggrepel))

#convert to sce
dssce <- Seurat::as.SingleCellExperiment(ds)

# put dimensional reductions in place
SingleCellExperiment::reducedDim(dssce, "RPCA") <- ds@reductions$rpca@cell.embeddings
SingleCellExperiment::reducedDim(dssce, "UMAP") <- ds@reductions$umap@cell.embeddings

# compute neighborhoods
dsmilo <- miloR::Milo(dssce)
dsmilo <- miloR::buildGraph(dsmilo, k = 30, d = 30, reduced.dim = reduced.dims)
dsmilo <- miloR::makeNhoods(dsmilo, prop = 0.1, k = 30, d = 30, refined = T, reduced_dims = "RPCA")
dsmilo <- miloR::buildNhoodGraph(dsmilo)

# compute neighborhoods of reference

dssce_ref <- Seurat::as.SingleCellExperiment(ds_ref)

# put dimensional reductions in place
SingleCellExperiment::reducedDim(dssce_ref, "RPCA") <- ds_ref@reductions$rpca@cell.embeddings
SingleCellExperiment::reducedDim(dssce_ref, "UMAP") <- ds_ref@reductions$umap@cell.embeddings

dsmilo_ref <- miloR::Milo(dssce_ref)
dsmilo_ref <- miloR::buildGraph(dsmilo_ref, k = 30, d = 30, reduced.dim = reduced.dims)
dsmilo_ref <- miloR::makeNhoods(dsmilo_ref, prop = 0.1, k = 30, d = 30, refined = T, reduced_dims = "RPCA")
dsmilo_ref <- miloR::buildNhoodGraph(dsmilo_ref)

# get genes for similarity
genes <- ds_ref@assays$RNA@meta.data$var.features[!is.na(ds_ref@assays$RNA@meta.data$var.features)]
genes <- genes[-(grep("^MT", genes))]

# get maximum correlation as measure of similarity
nhsim <- scrabbitr::calcNhoodSim(dsmilo, dsmilo_ref, sim_features = sim_genes, sim_preprocessing="gene_spec", sim_measure="pearson", verbose = TRUE)
max_nhsim <- scrabbitr::getMaxMappings(nhsim$nhood_sim, 1, long_format=FALSE)

