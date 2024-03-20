# R script for analysis of bronchioid data derived from two donoers after 21days of culture

# Load data
ds1 <- Seurat::Read10x("path/to/patient1/filtered_feature_bc_matrix/")
ds2 <- Seurat::Read10x("path/to/patient2/filtered_feature_bc_matrix/")

# add patient meta.data
ds1@meta.data$Patient <- "patient1"
ds2@meta.data$Patient <- "patient2"

# merge data
ds <- merge(ds1, ds2)

# Aplly quality filtering
ds <- subset(ds, subset = nCount_RNA > 2000 & percent.mt < 10)

# Seurat analysis and data integration
ds <- Seurat::NormalizeData(ds)
ds <- Seurat::FindVariableFeatures(ds)
ds <- Seurat::ScaleData(ds)
ds <- Seurat::IntegrateLayers(
  ds,
  method = "RPCAIntegration",
  new.reduction = "rpca",
  verbose = TRUE
)
ds <- SeuratObject::JoinLayers(ds)
ds <- Seurat::FindNeighbors(ds, reduction = "rpca", dims = 1:30)
ds <- Seurat::RunUMAP(ds, reduction = "rpca", dims = 1:30)
ds <- Seurat::FindClusters(ds, resolution = 0.5)

# Annotate seurat_clusters
index <- c("Ciliated", "Basal", "Secretory", "Basal Differentiating", "Deuterosomal", "Goblet", "Proliferating")
names(index) <- c("0", "1", "2", "3", "4", "5", "6")
ds@meta.data$Manual.Annotation <- index[ds@meta.data$seurat_clusters]

# Get Ionocytes and Hillock-like cells with UMAP coordinates
plot <- Seurat::DimPlot(ds)
Ionocyte <- Seurat::CellSelector(plot) # select Ionocytes forming a separate cluster on UMAP
ds@meta.data[Ionocytes, ]$Manual.Annotation <- "Ionocyte"
Hillock_like <- Seurat::CellSelector(plot) # select Hillock-like cells forming a separate cluster on UMAP
ds@meta.data[Hillock_like, ]$Manual.Annotation <- "Hillock-like"

# CellTypist predicitions
# make input file for CellTypist
table <- Seurat::GetAssayData(ds, layer = "counts", assay = "RNA")
write.csv(table, "path/to/CellTypist/input/input.csv")

# Celltypist Python package was used as command line tool
# CellTypist models were downloaded using the command:
# celltypist --update-models

# CellTypist was then run using following command::
# celltypist --indata /path/to/input.csv --model Human_Lung_Atlas.pkl --majority-voting

# add Celltypist predictions
CT.predictions <- read.csv("path/to/predicted_labels.csv")
index <- CT.predictions$majority_voting
names(index) <- gsub("\\.", "-", CT.predictions$X)
ds@meta.data$Celltypist.prediction <- index[colnames(ds)]

# Grouping Celltypist predictions into broader categories to facilitate comparability with manual annotation
index <- c("Basal", "Basal", "Deuterosomal", "Multiciliated", "Secretory", "Secretory", "Ionocyte", "Goblet", "Hillock-like")
names(index) <- c("Suprabasal", "Basal resting", "Deuterosomal", "Multiciliated (non-nasal", "pre-TB secretory",
                  "SMG duct", "Ionocyte", "Goblet (nasal)", "Hillock-like")
ds@meta.data$Celltypist.prediction <- index[ds@meta.data$Celltypist.prediction]

# for ARI comparison make Multiciliated Ciliated, liek in teh Manual annotation
ds@meta.data$Celltypist.prediction.ari <- ds@meta.data$Celltypist.prediction
ds@meta.data$Celltypist.prediction.ari[ds@meta.data$Celltypist.prediction.ari == "Multiciliated"] <- "Ciliated"

# Annotate Basal Differentiating as Basal for ARi 
ds@meta.data$Manual.Annotation.ari <- ds@meta.data$Manual.Annotation
ds@meta.data$Manual.Annotation.ari[ds@meta.data$Manual.Annotation.ari == "Basal Differentiating"] <- "Basal"

# Compute Adjusted Rand Index between manual annotation and CellTypist prediciton
ARI <- mclust::adjustedRandIndex(ds_upload@meta.data$Celltypist.prediction.ari, ds_upload@meta.data$Manual.Annotation.ari)
