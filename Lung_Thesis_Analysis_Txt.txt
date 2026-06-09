#sCRNA Sequencing Data Analysis for CISH Lungs
#WT and KO CISH Lungs infected with 1000pfu of X31 IAV

#Upload data to seven bridges

install.packages("BiocManager")

BiocManager::install("sevenbridges")

library("sevenbridges")

#Generate Authentication token from seven bridges and create value for 'a'

a <- Auth(token = "85c8aeff9cae4859aa47941fcbf80c90", url = "https://api.sbgenomics.com/v2")

#create value for 'p' after searching for a$project()
p <- a$project(id = "rfayaz/cish-sequencing-d3")


p$upload("U:/Research/Stambas-Lab/Rifqa/CISH scRNAseq/Sequencing fastQ files", metadata = list(platform = "Illumina"))


a$project()

#Seurat old package, and make sure I am uploading the same things

#save and  load  workspace
("C:\\Users\\rfayaz\\OneDrive - Deakin University\\ScRNA CISH\\Thesis analysis\\Lung_Thesis_2")
save.image(file = "CISH_X31_RF")
load(file = "CISH_X31_RF")

#############DO NOT RUN WHATS IN BETWEEN#############################
update.packages(ask = FALSE, checkBuilt = TRUE)
old.packages()
update.packages(ask = FALSE, checkBuilt = TRUE, lib = "C:/Users/rfayaz/AppData/Local/R/win-library/4.4")
.libPaths()


packageVersion("rlang")  # should show >= 1.1.7

packageVersion("fs")     # should show >= 1.6.7

packageVersion("Rcpp")   # should show >= 1.1.0

#############DO NOT RUN WHATS IN BETWEEN#############################

#install the necessary packages 
BiocManager::install("ComplexHeatmap")
install.packages("installr")
install.packages("devtools")
install.packages('gprofiler2')
install.packages("rlang")
BiocManager::install("SingleR")
BiocManager::install("cli")
BiocManager::install("dplyr")
BiocManager::install("patchwork")
BiocManager::install("tidyverse")
BiocManager::install("data.table")
BiocManager::install("cowplot")
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("dittoSeq")
BiocManager::install("AnnotationHub")
BiocManager::install("EnhancedVolcano")
BiocManager::install("AnnotationHub")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("ensembldb")
BiocManager::install("fgsea")
BiocManager::install("GOplot")
BiocManager::install("circlize")
BiocManager::install("DESeq2")
BiocManager::install("SingleCellExperiment")
BiocManager::install("MAST")
BiocManager::install("gprofiler2")
BiocManager::install("clustree")
BiocManager::install("ggpubr")
BiocManager::install('multtest')
BiocManager::install('metap')
BiocManager::install('ggrastr')
BiocManager::install('scatterpie')
BiocManager::install('tidydr')
BiocManager::install('enrichplot')
BiocManager::install('openxlsx')
BiocManager::install('writexl')
BiocManager::install("celldex")
BiocManager::install("dittoSeq")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("clusterProfiler")
BiocManager::install("ggplot2")
BiocManager::install("ensembldb")
BiocManager::install("DESeq2")
BiocManager::install("MAST")
BiocManager::install("DT")
install.packages("ggvenn")


# load all necessary packages
library(rlang)
library(fs)
library(Rcpp)
library(devtools)
library(installr)
library(ggvenn)
library(stringr) 
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(data.table)
library(cowplot)
library(SingleR)
library(celldex)
library(dittoSeq)
library(EnhancedVolcano)
library(AnnotationHub)
library(org.Mm.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(ensembldb)
library(fgsea)
library(GOplot)
library(circlize)
library(DESeq2)
library(SingleCellExperiment)
library(MAST)
library(gprofiler2)
library(clustree)
library(ggplot2)
library(ggpubr)
library(org.Mm.eg.db)
library(ggrastr)
library(scatterpie)
library(tidydr)
library(enrichplot)
library(openxlsx)
library(writexl)
library(RColorBrewer)
library(DOSE)
library(ggrastr)
library(tidydr)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(DT)

#loadGeneCount
CISHWT_L_count <- Read10X(data.dir= "U:\\Rifqa\\PhD Documents\\ScRNA_Sequencing\\CISH_scRNAseq\\7Bridges_NewData_Aug\\LUNGS\\CISHWT_LUNG")
CISHKO_L_count <- Read10X(data.dir= "U:\\Rifqa\\PhD Documents\\ScRNA_Sequencing\\CISH_scRNAseq\\7Bridges_NewData_Aug\\LUNGS\\CISHKO_LUNG")


# create CISHWT and CISHKO as Seurat object
options(Seurat.object.assay.version = "v3")
CISHWT_Seu <- CreateSeuratObject(CISHWT_L_count, min.cells = 3, min.features = 150, project = "CISHWT")
class(CISHWT_Seu[["RNA"]])
CISHWT_Seu$Type <- "WT"
CISHWT_Seu <- RenameCells(CISHWT_Seu, add.cell.id = "WT")

CISHKO_Seu <- CreateSeuratObject(CISHKO_L_count, min.cells = 3, min.features = 150, project = "CISHKO")
class(CISHKO_Seu[["RNA"]])
CISHKO_Seu$Type <- "KO"
CISHKO_Seu <- RenameCells(CISHKO_Seu, add.cell.id = "KO")

sample_cols <- c("CISHWT" = "grey60", "CISHKO" = "lightblue", "WT" = "grey60", "KO" = "lightblue")

print(paste0("Cells in CISHWT_Seu: ", ncol(CISHWT_Seu)))
print(paste0("Cells in CISHKO_Seu: ", ncol(CISHKO_Seu)))


# Merge two data for first QC
CISH_merged <- merge(CISHWT_Seu, CISHKO_Seu) 
CISH_merged[["percent.mt"]] <- PercentageFeatureSet(CISH_merged, pattern = "^mt-")
CISH_merged[["percent.rbc"]] <- PercentageFeatureSet(CISH_merged, pattern = "^Hb[ab]-")
Fig_QC_A_before <- VlnPlot(CISH_merged, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbc"), group.by = "orig.ident", ncol = 4, pt.size = 0.1) & scale_fill_manual(values = sample_cols)
Fig_QC_A_before


# Numbers to help set up cut-off gating
table(CISH_merged$orig.ident) # Number of cells in each condition
summary(CISH_merged$nCount_RNA) # Transcript number
summary(CISH_merged$nFeature_RNA) # Number of unique genes expressed
summary(CISH_merged$percent.mt) # Percent of transcripts of mitochondrial origin
summary(CISH_merged$percent.rbc) # Percent of transcripts of RBC origin

sd(CISH_merged$nCount_RNA) 
sd(CISH_merged$nFeature_RNA)
sd(CISH_merged$percent.mt) 

p1 <- ggplot(CISH_merged@meta.data, mapping = aes(x = nFeature_RNA)) + geom_histogram(bins = 120) + geom_vline(aes(xintercept = 4000, color = 'red')) + NoLegend()
p2 <- ggplot(CISH_merged@meta.data, mapping = aes(x = percent.mt)) + geom_histogram(bins = 120) + geom_vline(aes(xintercept = 20, color = 'red')) + NoLegend()
p3 <- ggplot(CISH_merged@meta.data, mapping = aes(x = nCount_RNA)) + geom_histogram(bins = 120) + geom_vline(aes(xintercept = 15000, color = 'red')) + NoLegend()
print(p1 + p3 + p2)


#After applying QC thresholds, distribution of quality metrics is roughly normal. Another round of cluster based QC is applied downstream.
CISH_merged.f <- subset(CISH_merged, subset = percent.mt < 20 & percent.rbc < 5 & nFeature_RNA < 4000 & nFeature_RNA > 300 & nCount_RNA < 15000)
p1 <- ggplot(CISH_merged.f@meta.data, mapping = aes(x = nFeature_RNA)) + geom_histogram(bins = 120) + geom_vline(aes(xintercept = 4000), color = 'red') + NoLegend()
p2 <- ggplot(CISH_merged.f@meta.data, mapping = aes(x = percent.mt)) + geom_histogram(bins = 120) + geom_vline(aes(xintercept = 20, color = 'red')) + NoLegend()
p3 <- ggplot(CISH_merged.f@meta.data, mapping = aes(x = nCount_RNA)) + geom_histogram(bins = 120) + geom_vline(aes(xintercept = 15000, color = 'red')) + NoLegend()
print(p1 + p3 + p2) 


Fig_QC_B_after <- VlnPlot(CISH_merged.f, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rbc"), group.by = "orig.ident", ncol = 4, pt.size = 0.1) & scale_fill_manual(values = sample_cols)
Fig_QC_B_after

#inspect summary statistics


# Inspect cell numbers after creating Seurat objects
print(paste0("Cells in CISHWT_Seu: ", ncol(CISHWT_Seu)))
print(paste0("Cells in CISHKO_Seu: ", ncol(CISHKO_Seu)))
print(paste0("Total cells before merge: ", ncol(CISHWT_Seu) + ncol(CISHKO_Seu)))

# Inspect cell numbers after merging
print(paste0("Total cells in CISH_merged: ", ncol(CISH_merged)))
print(table(CISH_merged$orig.ident))
print(table(CISH_merged$Type))

# Inspect cell numbers after initial QC filtering
print(paste0("Total cells after QC filtering: ", ncol(CISH_merged.f)))
print(table(CISH_merged.f$orig.ident))
print(table(CISH_merged.f$Type))

# Inspect number and percent of cells kept after QC
cells_before_qc <- ncol(CISH_merged)
cells_after_qc <- ncol(CISH_merged.f)
cells_removed_qc <- cells_before_qc - cells_after_qc
percent_kept_qc <- round((cells_after_qc / cells_before_qc) * 100, 2)
percent_removed_qc <- round((cells_removed_qc / cells_before_qc) * 100, 2)

print(paste0("Cells before QC: ", cells_before_qc))
print(paste0("Cells after QC: ", cells_after_qc))
print(paste0("Cells removed by QC: ", cells_removed_qc))
print(paste0("Percent kept after QC: ", percent_kept_qc, "%"))
print(paste0("Percent removed by QC: ", percent_removed_qc, "%"))

# Inspect WT/KO cells kept and removed after QC
qc_cell_summary <- data.frame(
  Sample = names(table(CISH_merged$orig.ident)),
  Before_QC = as.numeric(table(CISH_merged$orig.ident)),
  After_QC = as.numeric(table(CISH_merged.f$orig.ident))
)

qc_cell_summary$Removed_QC <- qc_cell_summary$Before_QC - qc_cell_summary$After_QC
qc_cell_summary$Percent_Kept <- round((qc_cell_summary$After_QC / qc_cell_summary$Before_QC) * 100, 2)
qc_cell_summary$Percent_Removed <- round((qc_cell_summary$Removed_QC / qc_cell_summary$Before_QC) * 100, 2)

print(qc_cell_summary)
library(tidyr); library(dplyr); library(ggplot2)
#plot the summaries
library(tidyr); library(dplyr); library(ggplot2); library(patchwork)

qc_plot <- qc_cell_summary %>%
  pivot_longer(cols = c(Before_QC, After_QC, Removed_QC, Percent_Kept), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         Before_QC = "Before QC",
                         After_QC = "After QC",
                         Removed_QC = "Removed by QC",
                         Percent_Kept = "% Kept"
  ))

facet_labels_counts <- c(
  "Before QC" = "Before QC\nTotal cells before any filtering",
  "After QC" = "After QC\nCells remaining after filtering",
  "Removed by QC" = "Removed by QC\nCells excluded by QC thresholds"
)

facet_labels_pct <- c(
  "% Kept" = "% Kept\nPercentage of cells retained after QC"
)
qc_counts <- qc_plot %>% dplyr::filter(Metric %in% c("Before QC","After QC","Removed by QC"))
qc_percent <- qc_plot %>% dplyr::filter(Metric == "% Kept")

p_counts <- ggplot(qc_counts, aes(x = Sample, y = Value, fill = Sample)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Value), vjust = -0.3, size = 3.5) +
  facet_wrap(~Metric, scales = "fixed", nrow = 1, labeller = labeller(Metric = facet_labels_counts)) +
  scale_fill_manual(values = c("CISHWT" = "grey60", "CISHKO" = "lightblue")) +
  coord_cartesian(ylim = c(0, 4000)) +
  labs(title = "QC Cell Summary", x = "Sample", y = "Cell Number", fill = "Sample") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, face = "bold"))

p_percent <- ggplot(qc_percent, aes(x = Sample, y = Value, fill = Sample)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Value, "%")), vjust = -0.3, size = 3.5) +
  facet_wrap(~Metric, scales = "fixed", nrow = 1, labeller = labeller(Metric = facet_labels_pct)) +
  scale_fill_manual(values = c("CISHWT" = "grey60", "CISHKO" = "lightblue")) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Sample", y = "Percent", fill = "Sample") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, face = "bold"))

p_counts | p_percent
## end of initial QC
##Cish expression
## end of initial QC
## Cish expression
## USE RNA here because this is gene expression
DefaultAssay(CISH_merged.f) <- "RNA"
CISH_merged.f <- NormalizeData(CISH_merged.f, verbose = FALSE)

VlnPlot(CISH_merged.f, features = "Cish", group.by = "Type", pt.size = 0.1) + NoLegend()
DotPlot(CISH_merged.f, features = "Cish", group.by = "Type") + RotatedAxis()


# Pre-integration PCA / UMAP on QC-filtered merged object
# USE RNA here because this is before integration
DefaultAssay(CISH_merged.f) <- "RNA"

CISH_merged.f <- CISH_merged.f %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(verbose = FALSE) %>% 
  ScaleData(features = VariableFeatures(CISH_merged.f), verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)

## Cell-cycle scoring based on known genes
## USE RNA here because cell-cycle scoring uses gene expression
DefaultAssay(CISH_merged.f) <- "RNA"

mmus_s = gorth(cc.genes.updated.2019$s.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name
mmus_g2m = gorth(cc.genes.updated.2019$g2m.genes, source_organism = "hsapiens", target_organism = "mmusculus")$ortholog_name

CISH_merged.f <- CellCycleScoring(CISH_merged.f, s.features = mmus_s, g2m.features = mmus_g2m, set.ident = TRUE)
CISH_merged.f$CC.Difference <- CISH_merged.f$S.Score - CISH_merged.f$G2M.Score
head(CISH_merged.f)

DimPlot(CISH_merged.f, label = TRUE, split.by = "Phase") + NoLegend()

Fig_QC_CellCycle <- DimPlot(CISH_merged.f, reduction = "umap", group.by = "Phase") +
  ggtitle("UMAP Colored by Cell-Cycle Phase") +
  theme_classic()

Fig_QC_CellCycle


# Calculate % variance explained by PC1 and PC2
p1 <- DimPlot(CISH_merged.f, reduction = "pca", group.by = "Type", shuffle = TRUE)
p2 <- DimPlot(CISH_merged.f, reduction = "umap", group.by = "Type", shuffle = TRUE)
p3 <- ElbowPlot(CISH_merged.f)

percent_var = ((CISH_merged.f@reductions$pca@stdev)^2 / sum(rowVars(CISH_merged.f@assays$RNA@scale.data))) * 100
print(percent_var[1:2])
print(p1 + p2 + p3)

Fig_QC_PCA_Type <- DimPlot(CISH_merged.f, reduction = "pca", group.by = "Type", shuffle = TRUE) +
  scale_color_manual(values = c("WT" = "grey60", "KO" = "lightblue")) +
  ggtitle("PCA Colored by Genotype") +
  theme_classic()

Fig_QC_UMAP_Type <- DimPlot(CISH_merged.f, reduction = "umap", group.by = "Type", shuffle = TRUE) +
  scale_color_manual(values = c("WT" = "grey60", "KO" = "lightblue")) +
  ggtitle("UMAP Colored by Genotype") +
  theme_classic()

Fig_QC_PCA_Type + Fig_QC_UMAP_Type

print(DimHeatmap(CISH_merged.f, dims = 1:2, nfeatures = 50, balanced = TRUE, reduction = "pca"))


# Integration step
# Split merged object for integration
# Protein-coding gene list

hub <- AnnotationHub()
mouse <- query(hub, c("EnsDb", "Mus musculus"))
edb <- hub[["AH113713"]]

keys <- keys(edb, "GENENAME")
columns <- c("GENEID", "ENTREZID", "GENEBIOTYPE")

filter <- ~ gene_name %in% keys & gene_biotype == "protein_coding"
tbl <- ensembldb::select(edb, filter, columns) %>% as_tibble()

mRNA_gene_list <- unique(tbl$GENENAME)

# Integration step
CISH_merged.list <- SplitObject(CISH_merged.f, split.by = "orig.ident")

CISH_merged.list <- lapply(CISH_merged.list, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  return(x)
})

features <- SelectIntegrationFeatures(object.list = CISH_merged.list, nfeatures = 3000)

var_regex <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

features_filtered <- grep(var_regex, features, invert = TRUE, value = TRUE)
features_filtered <- intersect(features_filtered, mRNA_gene_list)

CISH_merged.list <- lapply(CISH_merged.list, function(x) {
  x <- ScaleData(x, features = features_filtered, verbose = FALSE)
  x <- RunPCA(x, features = features_filtered, verbose = FALSE)
  return(x)
})

CISH_merged.anchors <- FindIntegrationAnchors(
  object.list = CISH_merged.list,
  anchor.features = features_filtered,
  reduction = "rpca",
  dims = 1:20,
  verbose = FALSE
)

CISH_integrated <- IntegrateData(
  anchorset = CISH_merged.anchors,
  dims = 1:20,
  verbose = FALSE
)

# After integration, use integrated assay for clustering / dimensional reduction
DefaultAssay(CISH_integrated) <- "integrated"

CISH_integrated <- ScaleData(CISH_integrated, verbose = FALSE)

CISH_integrated <- RunPCA(
  CISH_integrated,
  npcs = 20,
  features = features_filtered,
  verbose = FALSE
)

CISH_integrated <- RunUMAP(CISH_integrated, dims = 1:20, verbose = FALSE)
CISH_integrated <- FindNeighbors(CISH_integrated, dims = 1:20, verbose = FALSE)
CISH_integrated <- FindClusters(CISH_integrated, resolution = 1.2, verbose = FALSE)


# Resolution testing
# USE integrated here because this tests clustering resolution
DefaultAssay(CISH_integrated) <- "integrated"

object <- Seurat::FindClusters(object = CISH_integrated, resolution = c(0.1, 0.5, 0.8, 1.0, 1.2))


head(object)
clustree(object)

colnames(object@meta.data)
grep("snn_res", colnames(object@meta.data), value = TRUE)

DimPlot(object, reduction = "umap", group.by = c("integrated_snn_res.0.1", "integrated_snn_res.0.5", "integrated_snn_res.0.8", "integrated_snn_res.1", "integrated_snn_res.1.2"))

Idents(CISH_integrated) <- "integrated_snn_res.1.2"

# Check final object
p4 <- ElbowPlot(CISH_integrated, n = 20)
CISH_integrated


# Final integrated UMAP plots
# USE integrated here because these are UMAP / clustering plots
DefaultAssay(CISH_integrated) <- "integrated"

p1 <- DimPlot(CISH_integrated, reduction = "umap", group.by = "orig.ident", shuffle = TRUE) + ggtitle("Mouse type")
p2 <- DimPlot(CISH_integrated, reduction = "umap", split.by = "orig.ident", repel = TRUE, shuffle = TRUE)
p3 <- DimPlot(CISH_integrated, reduction = "umap", group.by = "Phase", repel = TRUE, shuffle = TRUE)

print(p1 + p2 + p3)


# Contaminant marker checking
# USE RNA here because this checks real gene expression
DefaultAssay(CISH_integrated) <- "RNA"

FeaturePlot(CISH_integrated, features = c("Ptprc", "Pecam1", "Cdh5", "Col1a1", "Pdgfra", "Epcam"), ncol = 3)

library(patchwork)

p_feat <- FeaturePlot(CISH_integrated, features = c("Ptprc", "Pecam1", "Cdh5", "Col1a1", "Pdgfra", "Epcam","Rhd", "Hemgn", "Stmn1", "Hbb-bt"), ncol = 3)
p_feat + plot_annotation(
  caption = paste(
    "Marker guide:",
    "Ptprc = pan-immune cells (CD45)",
    "Pecam1 = endothelial cells",
    "Cdh5 = endothelial cells",
    "Col1a1 = fibroblast/stromal cells",
    "Pdgfra = fibroblast/stromal cells",
    "Epcam = epithelial cells",
    "Rhd = erythroid / red blood cell lineage",
    "Hemgn = erythroid / hematopoietic progenitor marker",
    "Stmn1 = proliferating / cycling cells",
    "Hbb-bt = red blood cells / hemoglobin contamination",
    sep = "\n"
  )
) & theme(plot.caption = element_text(size = 10, hjust = 0))


# Show integrated UMAP after marker checking
# USE integrated here because this is cluster visualization

DefaultAssay(CISH_integrated) <- "integrated"

DimPlot(CISH_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)

# Marker identification
# USE RNA here because FindAllMarkers should use real gene expression
DefaultAssay(CISH_integrated) <- "RNA"

# Marker identification using filtered protein-coding genes only
Idents(CISH_integrated) <- "integrated_snn_res.1.2"

marker_features <- rownames(CISH_integrated)
marker_features <- intersect(marker_features, mRNA_gene_list)

remove_marker_regex <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"
marker_features <- grep(remove_marker_regex, marker_features, invert = TRUE, value = TRUE)

Cluster.markers <- FindAllMarkers(
  CISH_integrated,
  features = marker_features,
  only.pos = TRUE,
  logfc.threshold = 0.1,
  test.use = "roc"
)

head(Cluster.markers)
table(Cluster.markers$cluster)

CISH_Lung_gene_Cluster1 <- Cluster.markers %>% 
  group_by(cluster) %>% 
  slice_max(order_by = avg_log2FC, n = 20)

View(CISH_Lung_gene_Cluster1)
write.csv((CISH_Lung_gene_Cluster1),"CISH_Lung_gene_Cluster1.csv")

FeaturePlot(CISH_integrated, feature = c ("Cd3e", "Nkg7", "Klrb1c"))


# RBC / erythroid and cell-cycle marker checking
# USE RNA here because this is gene expression
DefaultAssay(CISH_integrated) <- "RNA"

FeaturePlot(CISH_integrated, features = c("Rhd", "Hemgn", "Stmn1", "Hbb-bt"), ncol = 2, order = TRUE)

# UMAP / clustering display
# USE integrated here because this is cluster visualization
DefaultAssay(CISH_integrated) <- "integrated"

DimPlot(CISH_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)
DimPlot(CISH_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)

######
#removeclusters
# Edit this only after confirming the clusters are unwanted.
Idents(CISH_integrated) <- "integrated_snn_res.1.2"
table(Idents(CISH_integrated))
clusters_to_remove <- c(3,12,13)

CISH_integrated <- subset(CISH_integrated, idents = clusters_to_remove, invert = TRUE)

DimPlot(CISH_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)

# 1. Integrated assay → dimensional reduction + clustering ONLY
DefaultAssay(CISH_integrated) <- "integrated"

CISH_integrated <- ScaleData(CISH_integrated, verbose = FALSE)
CISH_integrated <- RunPCA(CISH_integrated, npcs = 20, verbose = FALSE)
CISH_integrated <- RunUMAP(CISH_integrated, dims = 1:20, verbose = FALSE)
CISH_integrated <- FindNeighbors(CISH_integrated, dims = 1:20, verbose = FALSE)
CISH_integrated <- FindClusters(CISH_integrated, resolution = 1.2, verbose = FALSE)

Idents(CISH_integrated) <- "seurat_clusters"


# 2. RNA assay → normalization + gene expression work
DefaultAssay(CISH_integrated) <- "RNA"

CISH_integrated <- NormalizeData(CISH_integrated, verbose = FALSE)
CISH_integrated <- FindVariableFeatures(CISH_integrated, nfeatures = 3000, verbose = FALSE)

colnames(object@meta.data)
colnames(CISH_integrated@meta.data)

#View both KO and WT combined     
DimPlot(CISH_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)
DimPlot(CISH_integrated, reduction = "umap", split.by =  "orig.ident", label = TRUE, repel = TRUE, shuffle = TRUE)

FeaturePlot(CISH_integrated, feature = c ("Cd3e", "Nkg7", "Klrb1c"))

#Run clusterfinder again 
#Run clusterfinder again 
# Run cluster marker finder again using filtered protein-coding genes only
DefaultAssay(CISH_integrated) <- "RNA"
Idents(CISH_integrated) <- "seurat_clusters"

marker_features <- rownames(CISH_integrated)
marker_features <- intersect(marker_features, mRNA_gene_list)

remove_marker_regex <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"
marker_features <- grep(remove_marker_regex, marker_features, invert = TRUE, value = TRUE)

Cluster.markers <- FindAllMarkers(
  CISH_integrated,
  features = marker_features,
  only.pos = TRUE,
  logfc.threshold = 0.1,
  test.use = "roc"
)

CISH_Lung_gene_Cluster <- Cluster.markers %>% 
  group_by(cluster) %>% 
  top_n(n = 20, wt = avg_log2FC)

View(CISH_Lung_gene_Cluster)
DimPlot(CISH_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)
FeaturePlot(CISH_integrated, feature = c ("Cd3e", "Nkg7", "Klrb1c"))
setwd("C:\\Users\\rfayaz\\OneDrive - Deakin University\\ScRNA CISH\\Thesis analysis\\Lung_Thesis_2")
write.csv((CISH_Lung_gene_Cluster),"CISH_Lung_gene_Cluster.csv")

######
DefaultAssay(CISH_integrated) <- "RNA"

#Cell Cluster Annotation
DimPlot(CISH_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)

FeaturePlot(CISH_integrated, feature = c ("Cish"))

FeaturePlot(CISH_integrated, feature = c ("Cd3e", "Cd3g","Cd28","Cd5")) #Tcells

FeaturePlot(CISH_integrated, feature = c ("Cd4", "Cd8b1", "Cd8a", "Cd3e", "Il2rb","Il2ra")) #Cd4CD8

FeaturePlot(CISH_integrated, feature = c ("Cd3e")) 

FeaturePlot(CISH_integrated, features = c("Ncr1", "Gzma", "Klrb1c", "Nkg7")) #NK Cells

FeaturePlot(CISH_integrated, feature = c ("Cd4", "Cd3e", "Cd8a"))

FeaturePlot(CISH_integrated, feature = c ("Csf1r","C1qb", "Ms4a6c", "Fcgr1")) #MacsandMonos

FeaturePlot(CISH_integrated, features = c("Ly6c2", "Ccr2", "Cx3cr1", "Ms4a7", "C1qb", "C1qc"))
FeaturePlot(CISH_integrated, feature = c ("C1qc", "C1qb", "Lgmn", "Ms4a7")) #IMs

FeaturePlot(CISH_integrated, feature = c ("Cst3", "Ifitm6","Clec4a3", "Treml4")) #iMons

FeaturePlot(CISH_integrated, feature = c ("Ms4a1","Ighm", "Iglc2", "H2-Eb1")) #Bcells

FeaturePlot(CISH_integrated, feature = c ("S100a9", "Csf3r", "Mmp9", "Retnlg")) #Neutrophils

FeaturePlot(CISH_integrated, feature = c ("Ccr2", "Ly6c2")) #Classical mons

FeaturePlot(CISH_integrated, feature = c ("Cx3cr1", "Cd43", "Spn")) #Non-Classical mons

FeaturePlot(CISH_integrated, feature = c ("Siglech", "Flt3", "Zbtb46", "Rogdi", "Cd209a", "Ccr7")) #Dcs

FeaturePlot(CISH_integrated, feature = c ("Siglecf")) #AMs

FeaturePlot(CISH_integrated, feature = c ("Sftpc", "Lamp3", "Hopx", "Ager"))

FeaturePlot(CISH_integrated, feature = c ("Col1a1", "Col1a2", "Sod3", "Col6a2")) #AFs

FeaturePlot(CISH_integrated, feature = c ("Ccr2"))


####################### SINGLER WORKFLOW OPTIONAL#######################
library(SingleCellExperiment)

# Convert your integrated Seurat object to SCE, without touching the original
CISH_sce_for_SingleR <- as.SingleCellExperiment(CISH_integrated)


#Load reference dataset
library(celldex)

# Prebuilt ImmGen reference for mouse
immgen_ref <- ImmGenData()

library(SingleR)

singleR_results <- SingleR(
  test = CISH_sce_for_SingleR,  # separate copy
  ref = immgen_ref,
  labels = immgen_ref$label.main
)

colnames(CISH_integrated@meta.data)

# See predicted cell types
head(singleR_results$labels)
table(singleR_results$labels)

# Copy labels to Seurat object without modifying existing structure
CISH_integrated$SingleR_label <- singleR_results$labels

table(CISH_integrated$seurat_clusters, CISH_integrated$SingleR_label)

#you can create a heatmap to visualise better
library(pheatmap)

# Create a cross-tab of clusters vs SingleR labels
cluster_label_table <- table(CISH_integrated$seurat_clusters, CISH_integrated$SingleR_label)

# Optional: convert counts to proportions per cluster
cluster_label_prop <- prop.table(cluster_label_table, margin = 1)

my_colors <- colorRampPalette(c("white", "blue", "red"))(100)
pheatmap(cluster_label_prop,
         color = my_colors,
         cluster_rows = FALSE,      # no clustering for clusters
         cluster_cols = FALSE,      # no clustering for cell types
         show_rownames = TRUE,      # clusters on y-axis
         show_colnames = TRUE,      # cell types on x-axis
         angle_col = 45,            # rotate x-axis labels for readability
         border_color = NA,         # remove grid lines
         fontsize_row = 12,         # cluster labels
         fontsize_col = 12,         # cell type labels
         cellwidth = 25,            # adjust cell width
         cellheight = 15,           # adjust row height
         main = "Cluster vs Cell Type Proportions",
         legend = TRUE)

##do deeper analysis 
# ---------------------------
# Load libraries
# ---------------------------
library(Seurat)
library(SingleCellExperiment)
library(celldex)
library(SingleR)
library(pheatmap)
library(dplyr)
library(tidyr)

# ---------------------------
# Convert Seurat -> SCE
# ---------------------------
CISH_sce_for_SingleR_lung <- as.SingleCellExperiment(CISH_integrated)

# ---------------------------
# Load ImmGen reference
# ---------------------------
immgen_ref <- ImmGenData()

# ---------------------------
# Run SingleR (fine labels for CD4/CD8 resolution)
# ---------------------------
singleR_results <- SingleR(
  test = CISH_sce_for_SingleR_lung,
  ref = immgen_ref,
  labels = immgen_ref$label.fine
)

# Inspect predictions
head(singleR_results$labels)
table(singleR_results$labels)

# ---------------------------
# Add labels to Seurat object
# ---------------------------
CISH_integrated$SingleR_label <- singleR_results$labels

# ---------------------------
# Cluster vs Cell Type counts (long format)
# ---------------------------
cluster_type_counts <- as.data.frame(table(
  Cluster = CISH_integrated$seurat_clusters,
  CellType = CISH_integrated$SingleR_label
))

# Remove zero entries and sort
cluster_type_counts <- cluster_type_counts[cluster_type_counts$Freq > 0, ]
cluster_type_counts <- cluster_type_counts[order(
  cluster_type_counts$Cluster,
  -cluster_type_counts$Freq
), ]

# View table
cluster_type_counts

# ---------------------------
# Wide format (matrix-style table)
# ---------------------------
cluster_type_matrix <- xtabs(
  Freq ~ Cluster + CellType,
  data = cluster_type_counts
)

cluster_type_matrix

# Add total cells per cluster
cluster_totals <- rowSums(cluster_type_matrix)
cluster_summary <- cbind(cluster_type_matrix, TotalCells = cluster_totals)

view(cluster_summary)

write.csv((cluster_summary),"cluster_summary.csv")


###################################################################################

DefaultAssay(CISH_integrated) <- "integrated"
Idents(CISH_integrated) <- "seurat_clusters"

# cell annotation and rename each cluster 

celltype=data.frame(ClusterID=0:16,
                    celltype='NA')
celltype[celltype$ClusterID %in% c( 0),2]='Neutrophils'
celltype[celltype$ClusterID %in% c( 1),2]='B-cells'
celltype[celltype$ClusterID %in% c( 2),2]='Macrophages'
celltype[celltype$ClusterID %in% c( 3),2]='NK cells'
celltype[celltype$ClusterID %in% c( 4),2]='Macrophages'
celltype[celltype$ClusterID %in% c( 5),2]='T-cells'
celltype[celltype$ClusterID %in% c( 6),2]='Neutrophils'
celltype[celltype$ClusterID %in% c( 7),2]='T-cells'
celltype[celltype$ClusterID %in% c( 8),2]='Classical monocytes'
celltype[celltype$ClusterID %in% c( 9),2]='Classical monocytes'
celltype[celltype$ClusterID %in% c( 10),2]='Non-classical monocytes'
celltype[celltype$ClusterID %in% c( 11),2]='DCs'
celltype[celltype$ClusterID %in% c( 12),2]='Macrophages'
celltype[celltype$ClusterID %in% c( 13),2]='Macrophages'
celltype[celltype$ClusterID %in% c( 14),2]='Neutrophils'
celltype[celltype$ClusterID %in% c( 15),2]='DCs'
celltype[celltype$ClusterID %in% c( 16),2]='Other granulocytes'

table(CISH_integrated$seurat_clusters)
print(view(celltype))

table(celltype$celltype)
CISH_integrated@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  CISH_integrated@meta.data[which(CISH_integrated@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(CISH_integrated@meta.data$celltype)
head(CISH_integrated)

DefaultAssay(CISH_integrated) <- "integrated"
DimPlot(CISH_integrated, reduction = "umap", group.by = "celltype", 
        label = TRUE, repel = TRUE) + 
  ggtitle("Cell annotations")

pcluster1 <- DimPlot(CISH_integrated, 
                     reduction = "umap", 
                     label = TRUE, 
                     group.by = "celltype",
                     repel = TRUE, 
                     pt.size = 1.3,  # Adjust point size
                     alpha = 100, # Adjust transparency
                     label.size = 4.5) +
  ggtitle("Lung Immune Cell Types")

CISH_integrated$celltype_wrapped <- stringr::str_wrap(CISH_integrated$celltype, width = 14)

pcluster <- DimPlot(CISH_integrated, 
                    reduction = "umap", 
                    label = TRUE, 
                    group.by = "celltype_wrapped",
                    repel = TRUE, 
                    pt.size = 1.3,
                    alpha = 1,
                    label.size = 3.8) +
  ggtitle("Lung Immune Cell Types")

pcluster2 <- pcluster +
  theme_dr() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(hjust = 0, vjust = -1),  
    axis.title.y = element_text(hjust = 0, vjust = 1, angle = 90),  
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),  
    axis.text = element_blank(),
    axis.title = element_text(size = 15)
  ) +
  labs(title = "Lung Immune Cell Types") +
  guides(fill = "none")

pcluster2

pcluster_legend_only <- DimPlot(
  CISH_integrated,
  reduction = "umap",
  group.by = "celltype_wrapped",
  label = FALSE,
  repel = FALSE,
  pt.size = 1.3,
  alpha = 1
) +
  theme_dr() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.title = element_text(size = 15),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    legend.position = "right"
  ) +
  labs(title = "Lung Immune Cell Types")

pcluster_legend_only

## END OF CELL CLUSTERING ##

##FEATURES BY CELL TYPE
DefaultAssay(CISH_integrated) <- "RNA"

DotPlot(CISH_integrated, 
        features = c("Nkg7", "Cd3e", "Ms4a1", "S100a9", "Il1rn", "S100a4", 
                     "Adgre4", "Csf3r", "Ms4a7", "C1qc",
                     "Syngr2", "Rogdi", "Chil3", "Ccr2", "Ly6c2", "Cx3cr1", "Mcpt8", "Siglecf", "Prg2", 
                     "Rab44", "Hdc"), 
        group.by = "celltype") +
  scale_color_gradient(low = "lightblue", high = "darkblue") +  # Adjust color gradient
  labs(title = "Selected Features by Cell Type",  # Add title
       x = "Features", 
       y = "Cell Type") +
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(angle = 360, hjust = 1, size = 15),
        axis.title = element_text(size = 12),  # Adjust axis title size
        plot.title = element_text(hjust = 0.5, size = 14),  # Center and size the title
        legend.title = element_text(size = 12),  # Adjust legend title size
        legend.text = element_text(size = 12))  # Adjust legend text size


VlnPlot(
  CISH_integrated,
  features = c("Nkg7", "Cd3e", "Ms4a1", "S100a9", "Il1rn", "S100a4", 
               "Adgre4", "Csf3r", "Ms4a7", "C1qc",
               "Syngr2", "Rogdi", "Chil3", "Ccr2", "Ly6c2", "Cx3cr1", "Mcpt8", "Siglecf", "Prg2", 
               "Rab44", "Hdc"),
  group.by = "celltype",
  pt.size = 0,
  stack = TRUE,
  flip = TRUE,
  fill.by = "ident"
) +
  ggtitle("Selected Feature Expression by Cell Type") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 18),
    strip.text.y.right = element_text(angle = 20, hjust = 0, size = 18, face = "bold"),
    strip.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none"
  )

########################## END OF OPTIONAL ADDITIONAL UMAP STYLES###########################################


#cell percentage 
celltype_ratio <- CISH_integrated@meta.data %>%
  group_by(Type, celltype) %>%
  summarise(n = n()) %>%
  mutate(relative_freq = 100*n/sum(n))
celltype_ratio$celltype <- factor(celltype_ratio$celltype)


# 2. Define custom colors using hex codes (replace with any you like)
custom_colors <- c(
  "WT" = "#D3D3D3",      # light grey
  "KO" = "#8dc7cb"       # turquoise/teal
)

# 3. Create frequency plot
freq_plot <- ggplot(celltype_ratio, aes(x = celltype, y = relative_freq)) + 
  geom_col(aes(fill = Type), color = "black", position = "dodge", width = 0.8) +
  ylab("Relative Frequency (%)") +
  scale_fill_manual(values = custom_colors) +  # use your custom colors here
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~ celltype, ncol = 3, scales = "free") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 11),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_line(color = "black")  # makes x-axis visible
  )

# 4. Combine into a single row if desired
Frequencyplot <- plot_grid(freq_plot, align = "h", axis = "b", nrow = 1)

Frequencyplot

####
# Total cells per cluster
# Total cells per annotated cluster/cell type
# Total cells per cell type
celltype_counts <- CISH_integrated@meta.data %>%
  group_by(celltype) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(desc(n_cells))

celltype_counts
celltype_count_plot <- ggplot(celltype_counts, aes(x = reorder(celltype, -n_cells), y = n_cells, fill = celltype)) +
  geom_col(color = "black", width = 0.8) +
  geom_text(aes(label = n_cells), vjust = -0.3, size = 3.5) +
  theme_classic(base_size = 12) +
  labs(
    title = "Total Cells per Cell Type",
    x = "Cell Type",
    y = "Number of Cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

celltype_count_plot

# combined some clusters so explain in thesis

###DIFFERENTIAL GENE EXPRESSION ANALYSIS######
# differential gene  expression analysis (WT v.s KO)
#You need to modify test.use  with  multiple variables from "wilcox", "LR", “MAST" and ‘t’ and figure out what stat fit your sample.
setwd("C:\\Users\\rfayaz\\OneDrive - Deakin University\\ScRNA CISH\\Thesis analysis\\Lung_Thesis_2\\DEG")


DefaultAssay(CISH_integrated) <- "RNA"
CISH_integrated <- FindVariableFeatures(CISH_integrated, nfeatures = 3000, verbose = FALSE)
length(VariableFeatures(CISH_integrated))

Idents(CISH_integrated) <- "Type"


Deq_all_wilcox <- FindMarkers(CISH_integrated, ident.1 = "KO", ident.2 = "WT", logfc.threshold = 0, test.use = "wilcox")
Deq_all_MAST <- FindMarkers(CISH_integrated, ident.1 = "KO", ident.2 = "WT", logfc.threshold = 0, test.use = "MAST")
Deq_all_LR <- FindMarkers(CISH_integrated, ident.1 = "KO", ident.2 = "WT", logfc.threshold = 0, test.use = "LR")

print(view(Deq_all_wilcox))
print(view(Deq_all_LR))
print(view(Deq_all_MAST))

Deq_all_MAST$gene <- rownames(Deq_all_MAST)
view(Deq_all_MAST)

sig_genes_MAST <- subset(Deq_all_MAST, p_val_adj < 0.05)
print(view(sig_genes_MAST))

write.csv(Deq_all_MAST , "MAST_LUNG.csv")
write.csv(Deq_all_LR , "LR_LUNG.csv")

#now remove the non coding genes

pattern <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

Deq_all_wilcox$gene <- rownames(Deq_all_wilcox)
Deq_all_wilcox <- Deq_all_wilcox[!grepl(pattern, Deq_all_wilcox$gene), ]
View(Deq_all_wilcox)

sig_genes_wilcox <- subset(Deq_all_wilcox, p_val_adj < 0.05)
print(view(sig_genes_wilcox))

write.csv(Deq_all_wilcox , "Wilcox_LUNG.csv")
write.csv(sig_genes_wilcox , "sig_genes_wilcox_LUNG.csv")

#LR
Deq_all_LR$gene <- rownames(Deq_all_LR)
Deq_all_LR <- Deq_all_LR[!grepl(pattern, Deq_all_LR$gene), ]
View(Deq_all_LR)

sig_genes_LR <- subset(Deq_all_LR, p_val_adj < 0.05)
View(sig_genes_LR)

write.csv(Deq_all_LR, "LR_LUNG.csv")
write.csv(sig_genes_LR, "sig_genes_LR_LUNG.csv")

#Mast
Deq_all_MAST$gene <- rownames(Deq_all_MAST)
Deq_all_MAST <- Deq_all_MAST[!grepl(pattern, Deq_all_MAST$gene), ]
View(Deq_all_MAST)

sig_genes_MAST <- subset(Deq_all_MAST, p_val_adj < 0.05)
View(sig_genes_MAST)

write.csv(Deq_all_MAST, "MAST_LUNG.csv")
write.csv(sig_genes_MAST, "sig_genes_MAST_LUNG.csv")

# check how many genes are <0.05 p adjust in mast and wilcox
sum(sig_genes_wilcox$p_val_adj < 0.05)
sum(sig_genes_MAST$p_val_adj < 0.05)

#then check if that number is the same if the >0.25 in avglog2fcis added
sum(abs(sig_genes_wilcox$avg_log2FC) > 0.25) #only one gene dropped out
sum(abs(sig_genes_MAST$avg_log2FC) > 0.25)
sum(abs(sig_genes_LR$avg_log2FC) > 0.25)

#check how many genes are in common - are most genes in common?
common_sig_genes <- Reduce(intersect, list(sig_genes_wilcox$gene,sig_genes_MAST$gene,sig_genes_LR$gene))
common_sig_genes
length(common_sig_genes)

# Pairwise overlaps- check the amount of overlap
length(intersect(sig_genes_wilcox$gene, sig_genes_MAST$gene))
length(intersect(sig_genes_wilcox$gene, sig_genes_LR$gene))
length(intersect(sig_genes_MAST$gene, sig_genes_LR$gene))

# Percentage of each test that overlaps all 3 tests
length(common_sig_genes) / nrow(sig_genes_wilcox) * 100
length(common_sig_genes) / nrow(sig_genes_MAST) * 100
length(common_sig_genes) / nrow(sig_genes_LR) * 100


#do a summary table for overlap
overlap_summary <- data.frame(
  Test = c("Wilcox", "MAST", "LR"),
  Total_Significant_Genes = c(
    nrow(sig_genes_wilcox),
    nrow(sig_genes_MAST),
    nrow(sig_genes_LR)
  ),
  Common_With_All_3 = c(
    length(common_sig_genes),
    length(common_sig_genes),
    length(common_sig_genes)
  )
)

overlap_summary$Percent_Common <- round(
  overlap_summary$Common_With_All_3 / overlap_summary$Total_Significant_Genes * 100, 1
)

overlap_summary

library(gridExtra)
overlap_table_plot <- tableGrob(
  overlap_summary,
  rows = NULL,
  theme = ttheme_minimal(
    base_size = 12,
    core = list(fg_params = list(fontface = "plain")),
    colhead = list(fg_params = list(fontface = "bold"))
  )
)

grid::grid.newpage()
grid::grid.draw(overlap_table_plot)

# If Percent_Common is high, for example >70%, most significant genes overlap across tests
#this will be your high-confidence genes

#now make a table and graph of the high confidence genes
common_DEG_table <- Deq_all_MAST[Deq_all_MAST$gene %in% common_sig_genes, ]
common_DEG_table <- common_DEG_table[order(common_DEG_table$p_val_adj), ]
common_DEG_table[, c("gene", "avg_log2FC", "p_val_adj")]
View(common_DEG_table)
write.csv(common_DEG_table, "common_DEG_wilcox_MAST_LR_LUNG.csv")

library(Seurat)
DefaultAssay(CISH_integrated) <- "RNA"
CISH_integrated <- ScaleData(
  CISH_integrated,
  features = common_sig_genes,
  verbose = FALSE
)
DoHeatmap(
  CISH_integrated,
  features = common_sig_genes,
  group.by = "Type"
) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  ggtitle("Common DE Genes Across All Statistical Tests")

DefaultAssay(CISH_integrated) <- "RNA"
library(ggplot2)
common_DEG_table$Direction <- ifelse(common_DEG_table$avg_log2FC > 0, "Up in CISH KO", "Down CISH KO")

ggplot(common_DEG_table, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_col(color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("Up in CISH KO" = "#8dc7cb", "Down CISH KO" = "#D3D3D3")) +
  theme_classic(base_size = 12) +
  labs(
    title = "Core DEG Signature Genes Common Across All Statistical Tests",
    x = "Gene",
    y = "avg_log2FC",
    fill = "Direction"
  )


#seeing two cells with high CISH expression in KO, it might be ambient or background

Mast_sig_final <- subset(
  Deq_all_MAST,
  p_val_adj < 0.05 & 
    abs(avg_log2FC) > 0.25 &
    (pct.1 > 0.10 | pct.2 > 0.10)
)

View(Mast_sig_final)
write.csv(Mast_sig_final, "Mast_sig_final_LUNG.csv")

#using mast from here onwards
##############################################################################################

# Summary statistics for all filtered MAST genes

# set up threshold for log2foldchange
summary_positive <- summary(Deq_all_MAST$avg_log2FC[Deq_all_MAST$avg_log2FC > 0])

# Calculate the summary statistics for FC < 0
summary_negative <- summary(Deq_all_MAST$avg_log2FC[Deq_all_MAST$avg_log2FC < 0])
## > 0 means upregulated vice versa 

# Extract the 1st, mean, and 3rd quartile from each summary
# Summary stats for genes upregulated in KO

first_positive <- summary_positive[2]
mean_positive <- mean(Deq_all_MAST$avg_log2FC[Deq_all_MAST$avg_log2FC > 0])
third_positive <- summary_positive[5]

# Summary stats for genes downregulated in KO (i.e. up in WT)

first_negative <- summary_negative[2]
mean_negative <- mean(Deq_all_MAST$avg_log2FC[Deq_all_MAST$avg_log2FC < 0])
third_negative <- summary_negative[5]


# Print the summaries
cat("Summary for FC > 0:\n") #FC > 0 (genes upregulated in KO)
cat("1st Quartile:", first_positive, "\n")
cat("Mean:", mean_positive, "\n")
cat("3rd Quartile:", third_positive, "\n\n")

cat("Summary for FC < 0:\n") #FC < 0 (genes downregulated in KO, i.e., up in WT)
cat("1st Quartile:", first_negative, "\n")
cat("Mean:", mean_negative, "\n")
cat("3rd Quartile:", third_negative, "\n")

#3rd quartile for FC>0 which is upregulated and 1st for FC<0 which is downregulated
#Q3 for upregulated → captures the strongest upregulated genes (top 25%)
#Q1 for downregulated → captures the strongest downregulated genes (bottom 25%)
#“Take the most strongly upregulated genes as those in the 3rd quartile of FC>0, and the most strongly downregulated genes as those in the 1st quartile of FC<0.”
# “1st quartile” for negative values actually means the most negative fold changes (largest magnitude of downregulation), which is what you want.

# Assign group based on your quartiles

Deq_all_MAST$group <- case_when(
  Deq_all_MAST$p_val_adj < 0.05 & Deq_all_MAST$avg_log2FC > third_positive  ~ "up",
  Deq_all_MAST$p_val_adj < 0.05 & Deq_all_MAST$avg_log2FC < first_negative ~ "down",
  TRUE ~ "none"
)

# Assign significance based on adjusted p-value
Deq_all_MAST$Sig <- case_when(
  Deq_all_MAST$p_val_adj < 0.05 ~ "Sig",
  TRUE ~ "NS"
)

table(Deq_all_MAST$group)


##########
#Create a mapping table of gene symbols to Entrez IDs
mapped_genes <- data.frame(GeneName = rownames(Deq_all_MAST),
                           
                           ensemblID = mapIds(org.Mm.eg.db, keys =Deq_all_MAST$gene, keytype = "SYMBOL", column="ENTREZID"))

#Add Entrez IDs to your DE table
Deq_all_MAST$ensemblID <- mapped_genes$ensemblID

print(view(Deq_all_MAST))
write.csv((Deq_all_MAST),"DGE_MAST_CISHLUNG_EntrezID.csv")

# Replace zero p-values to avoid plotting errors

Deq_all_MAST$p_val_adj <- ifelse(
  Deq_all_MAST$p_val_adj == 0, 
  1e-304, 
  Deq_all_MAST$p_val_adj
)

Deq_all_MAST_top <- Deq_all_MAST

Deq_all_MAST_top <- Deq_all_MAST %>%
  dplyr::filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.1,
    pct.1 > 0.10 | pct.2 > 0.10
  )

print(view(Deq_all_MAST_top))

######################
##########################
# Function to get top N DEGs
##########################
get_top_DEGs <- function(de_table, n = 20, direction = c("up", "down")) {
  direction <- match.arg(direction)
  
  if (direction == "up") {
    genes <- rownames(de_table)[de_table$group == "up"]
    genes <- genes[order(-de_table[genes, "avg_log2FC"])]
  } else {
    genes <- rownames(de_table)[de_table$group == "down"]
    genes <- genes[order(de_table[genes, "avg_log2FC"])]
  }
  
  head(genes, n)
}
# Top upregulated (KO vs WT)
top_10_DEG_up <- get_top_DEGs(Deq_all_MAST_top, 10, "up")
top_20_DEG_up <- get_top_DEGs(Deq_all_MAST_top, 20, "up")
top_30_DEG_up <- get_top_DEGs(Deq_all_MAST_top, 30, "up")

# Top downregulated (up in WT)
top_10_DEG_down <- get_top_DEGs(Deq_all_MAST_top, 10, "down")
top_20_DEG_down <- get_top_DEGs(Deq_all_MAST_top, 20, "down")
top_30_DEG_down <- get_top_DEGs(Deq_all_MAST_top, 30, "down")

# Combine for volcano or table
top_10_DEG <- c(top_10_DEG_up, top_10_DEG_down)
top_20_DEG <- c(top_20_DEG_up, top_20_DEG_down)
top_30_DEG <- c(top_30_DEG_up, top_30_DEG_down)

##########################
# Create table with log2FC and direction
##########################
create_DEG_table <- function(genes, de_table) {
  data.frame(
    Gene = genes,
    log2FC = de_table[genes, "avg_log2FC"],
    p_val_adj = de_table[genes, "p_val_adj"],
    Direction = ifelse(de_table[genes, "avg_log2FC"] > 0, "Up in KO", "Up in WT")
  )
}

top_20_table <- create_DEG_table(top_20_DEG, Deq_all_MAST_top)
View(top_20_table)

##########################
# Volcano Plot
##########################

# Define custom colors based on quartiles
# Define custom colors using final DEG cutoff
keyvals <- ifelse(
  Deq_all_MAST_top$group == "up", "red",
  ifelse(Deq_all_MAST_top$group == "down", "royalblue", "grey80")
)

names(keyvals)[keyvals == "red"] <- "Up in KO"
names(keyvals)[keyvals == "royalblue"] <- "Up in WT"
names(keyvals)[keyvals == "grey80"] <- "NS"

Deq_plot <- EnhancedVolcano(
  Deq_all_MAST_top,
  x = "avg_log2FC",
  y = "p_val_adj",
  colCustom = keyvals,
  lab = rownames(Deq_all_MAST_top),
  selectLab = top_30_DEG,
  drawConnectors = TRUE,
  pCutoff = 0.05,
  FCcutoff = 0.25,
  max.overlaps = 60,
  title = "Lung KO vs WT"
)

Deq_plot

##########################
# Heatmap
##########################

# Build the heatmap matrix
heatmap_matrix <- data.frame(
  log2FC = Deq_all_MAST_top[top_20_DEG, "avg_log2FC"]
)
rownames(heatmap_matrix) <- top_20_DEG
heatmap_matrix <- as.matrix(heatmap_matrix)

breaks <- seq(min(heatmap_matrix), max(heatmap_matrix), length.out = 40)
colors <- colorRampPalette(c("blue", "white", "red"))(length(breaks)-1)

pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colors,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top Differentially Expressed Genes",
  legend = TRUE
)
table(Deq_all_MAST_top$group)

####


## Gene ontology and GSEA analysis is next
###############################
# 1. Check and set active assay
###############################

nrow(Deq_all_MAST_top)

sum(Deq_all_MAST_top$p_val_adj < 0.05)
sum(Deq_all_MAST_top$p_val_adj < 0.05 & abs(Deq_all_MAST_top$avg_log2FC) > 0.25)
sum(Deq_all_MAST_top$p_val_adj < 0.05 & abs(Deq_all_MAST_top$avg_log2FC) > 0.25 & Deq_all_MAST_top$pct.1 > 0.10)

# See all assays
names(CISH_integrated@assays)

# Check current active assay
DefaultAssay(CISH_integrated)

# For gene expression / DE / plotting expression, use RNA
DefaultAssay(CISH_integrated) <- "RNA"

Idents(CISH_integrated) <- "Type"

###############################
# 2. Prepare DEG table
###############################

# Use your existing MAST DEG table
# If Deq_all_MAST_top does not exist, create it from Deq_all_MAST
if (!exists("Deq_all_MAST_top")) {
  Deq_all_MAST_top <- Deq_all_MAST
}

# Add gene column if missing
if (!"gene" %in% colnames(Deq_all_MAST_top)) {
  Deq_all_MAST_top$gene <- rownames(Deq_all_MAST_top)
}

# Remove unwanted genes
remove_gene_pattern <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

Deq_all_MAST_top <- Deq_all_MAST_top %>%
  dplyr::filter(!grepl(remove_gene_pattern, gene))

# Replace zero adjusted p-values for plotting safety
Deq_all_MAST_top$p_val_adj <- ifelse(
  Deq_all_MAST_top$p_val_adj == 0,
  1e-304,
  Deq_all_MAST_top$p_val_adj
)

###############################
# 3. Define significant genes
###############################

sig_genes <- Deq_all_MAST_top %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  pull(gene)
length(sig_genes)

###############################
# 4. Map SYMBOL to ENSEMBL
###############################

mapped <- bitr(
  sig_genes,
  fromType = "SYMBOL",
  toType = "ENSEMBL",
  OrgDb = org.Mm.eg.db
)

mapped <- mapped %>%
  distinct(SYMBOL, .keep_all = TRUE)

head(mapped)
## blended featureplots

FeaturePlot(
  CISH_integrated,
  features = c("Lrrc8b", "Shb"),
  blend = TRUE,
  order = TRUE
)

###############################
# 5. GO enrichment: BP, MF, CC
###############################

GOBP_all <- enrichGO(
  gene          = mapped$ENSEMBL,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "BP",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  readable      = TRUE
)

GOMF_all <- enrichGO(
  gene          = mapped$ENSEMBL,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "MF",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  readable      = TRUE
)

GOCC_all <- enrichGO(
  gene          = mapped$ENSEMBL,
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENSEMBL",
  ont           = "CC",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  readable      = TRUE
)
###None found so dont run the below is nothing is found
###############################
# 6. GO dotplots
###############################

dotplot(GOBP_all, showCategory = 10) +
  ggtitle("GO Biological Process Enrichment in Lungs")

dotplot(GOMF_all, showCategory = 10) +
  ggtitle("GO Molecular Function Enrichment in Lungs")

dotplot(GOCC_all, showCategory = 10) +
  ggtitle("GO Cellular Component Enrichment in Lungs")

###############################
# 7. Prepare log2FC table for GO heatmap
###############################

log2FC_df <- Deq_all_MAST_top %>%
  dplyr::filter(gene %in% sig_genes) %>%
  dplyr::select(gene, avg_log2FC) %>%
  left_join(mapped, by = c("gene" = "SYMBOL"))

head(log2FC_df)

###############################
# 8. Function to calculate mean log2FC per GO term
###############################

make_GO_heatmap_matrix <- function(go_object, log2FC_df, top_n = 5) {
  
  go_df <- as.data.frame(go_object)
  
  if (nrow(go_df) == 0) {
    return(NULL)
  }
  
  go_df <- go_df %>%
    rowwise() %>%
    mutate(
      mean_log2FC = mean(
        log2FC_df$avg_log2FC[
          log2FC_df$gene %in% strsplit(geneID, "/")[[1]]
        ],
        na.rm = TRUE
      )
    ) %>%
    ungroup() %>%
    filter(!is.na(mean_log2FC))
  
  combined <- bind_rows(
    go_df %>% arrange(desc(mean_log2FC)) %>% slice_head(n = top_n),
    go_df %>% arrange(mean_log2FC) %>% slice_head(n = top_n)
  ) %>%
    distinct(Description, .keep_all = TRUE)
  
  heatmap_matrix <- matrix(
    combined$mean_log2FC,
    nrow = nrow(combined),
    ncol = 1,
    dimnames = list(combined$Description, "")
  )
  
  return(heatmap_matrix)
}

###############################
# 9. Create GO heatmap matrices
###############################

heatmap_matrix_BP <- make_GO_heatmap_matrix(GOBP_all, log2FC_df, top_n = 5)
heatmap_matrix_MF <- make_GO_heatmap_matrix(GOMF_all, log2FC_df, top_n = 5)
heatmap_matrix_CC <- make_GO_heatmap_matrix(GOCC_all, log2FC_df, top_n = 5)

###############################
# 10. Draw GO heatmaps
###############################

wrap_row_names <- function(names, width = 40) {
  sapply(names, function(x) paste(strwrap(x, width = width), collapse = "\n"))
}

all_values <- c(
  as.vector(heatmap_matrix_BP),
  as.vector(heatmap_matrix_MF),
  as.vector(heatmap_matrix_CC)
)

col_fun <- colorRamp2(
  c(min(all_values, na.rm = TRUE), 0, max(all_values, na.rm = TRUE)),
  c("blue", "white", "red")
)

rownames(heatmap_matrix_BP) <- wrap_row_names(rownames(heatmap_matrix_BP), width = 40)
rownames(heatmap_matrix_MF) <- wrap_row_names(rownames(heatmap_matrix_MF), width = 40)
rownames(heatmap_matrix_CC) <- wrap_row_names(rownames(heatmap_matrix_CC), width = 40)

p_BP <- Heatmap(
  heatmap_matrix_BP,
  name = "mean_log2FC",
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  border = TRUE,
  row_names_gp = gpar(fontsize = 8),
  row_names_side = "right",
  row_names_max_width = unit(15, "cm"),
  row_title = "Biological Process",
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 10, fontface = "bold")
)

p_MF <- Heatmap(
  heatmap_matrix_MF,
  name = "mean_log2FC",
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  border = TRUE,
  row_names_gp = gpar(fontsize = 8),
  row_names_side = "right",
  row_names_max_width = unit(15, "cm"),
  row_title = "Molecular Function",
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 10, fontface = "bold")
)

p_CC <- Heatmap(
  heatmap_matrix_CC,
  name = "mean_log2FC",
  col = col_fun,
  show_row_names = TRUE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  border = TRUE,
  row_names_gp = gpar(fontsize = 8),
  row_names_side = "right",
  row_names_max_width = unit(15, "cm"),
  row_title = "Cellular Component",
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 10, fontface = "bold")
)

combined_GO_heatmap <- p_BP %v% p_MF %v% p_CC

draw(combined_GO_heatmap, heatmap_legend_side = "left")

#######
## GSEA using MAST DEG results

DefaultAssay(CISH_integrated) <- "RNA"

Deq_all_MAST_pct_filtered <- Deq_all_MAST %>%
  dplyr::filter(
    pct.1 >= 0.10 | pct.2 >= 0.10
  )

View(Deq_all_MAST_pct_filtered)

# Use your filtered MAST DEG table
gsea_input <- Deq_all_MAST

# Make sure gene column exists
gsea_input$gene <- rownames(gsea_input)

# Remove unwanted genes
remove_gene_pattern <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

gsea_input <- gsea_input %>%
  dplyr::filter(!grepl(remove_gene_pattern, gene))

# Create ranked gene list using avg_log2FC
gene_list_GSEA <- gsea_input$avg_log2FC
names(gene_list_GSEA) <- gsea_input$gene

gene_list_GSEA <- gene_list_GSEA[!is.na(gene_list_GSEA)]
gene_list_GSEA <- gene_list_GSEA[!is.na(names(gene_list_GSEA))]
gene_list_GSEA <- gene_list_GSEA[!duplicated(names(gene_list_GSEA))]
gene_list_GSEA <- sort(gene_list_GSEA, decreasing = TRUE)

# Convert SYMBOL to ENSEMBL
mapping_GSEA <- bitr(
  names(gene_list_GSEA),
  fromType = "SYMBOL",
  toType = "ENSEMBL",
  OrgDb = org.Mm.eg.db
) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

gene_list_ens <- gene_list_GSEA[mapping_GSEA$SYMBOL]
names(gene_list_ens) <- mapping_GSEA$ENSEMBL

gene_list_ens <- gene_list_ens[!is.na(names(gene_list_ens))]
gene_list_ens <- gene_list_ens[!duplicated(names(gene_list_ens))]
gene_list_ens <- sort(gene_list_ens, decreasing = TRUE)

# Run GSEA GO Biological Process
gsea_BP <- gseGO(
  geneList      = gene_list_ens,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  keyType       = "ENSEMBL",
  minGSSize     = 15,
  maxGSSize     = 300,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = FALSE
)

# Significant results
gsea_BP_sig <- as.data.frame(gsea_BP) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust)

View(gsea_BP_sig)
dotplot(gsea_BP, showCategory = gsea_BP_sig$Description) +
  ggtitle("GSEA GO BP: KO vs WT")
##none 

##molecular function
gsea_MF <- gseGO(
  geneList      = gene_list_ens,
  OrgDb         = org.Mm.eg.db,
  ont           = "MF",
  keyType       = "ENSEMBL",
  minGSSize     = 15,
  maxGSSize     = 300,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = FALSE
)

gsea_MF_sig <- as.data.frame(gsea_MF) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust)

View(gsea_MF_sig)
dotplot(gsea_MF, showCategory = gsea_MF_sig$Description) +
  ggtitle("GSEA GO MF: KO vs WT")
##CC

gsea_CC <- gseGO(
  geneList      = gene_list_ens,
  OrgDb         = org.Mm.eg.db,
  ont           = "CC",
  keyType       = "ENSEMBL",
  minGSSize     = 15,
  maxGSSize     = 300,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = FALSE
)

gsea_CC_sig <- as.data.frame(gsea_CC) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust)

View(gsea_CC_sig)

dotplot(gsea_CC, showCategory = gsea_CC_sig$Description) +
  ggtitle("GSEA GO CC: KO vs WT")

# Dotplots only if you find any significant ones 

#next, subsetting

DefaultAssay(CISH_integrated) <- "RNA"

NK_obj <- subset(CISH_integrated, subset = celltype == "NK cells")
Idents(NK_obj) <- "Type"

NK_DEG <- FindMarkers(
  NK_obj,
  ident.1 = "KO",
  ident.2 = "WT",
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "MAST"
)

NK_DEG$gene <- rownames(NK_DEG)

# Remove unwanted genes
pattern <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

NK_DEG_filtered <- NK_DEG %>%
  dplyr::filter(!grepl(pattern, gene))

# Significant NK DEGs
NK_DEG_sig <- NK_DEG_filtered %>%
  dplyr::filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25,
    pct.1 > 0.10 | pct.2 > 0.10
  ) %>%
  dplyr::arrange(p_val_adj)

View(NK_DEG_filtered)
View(NK_DEG_sig)
write.csv((NK_DEG_filtered),"NK_DEG_all.csv")

###None found, try macrophages

DefaultAssay(CISH_integrated) <- "RNA"

# Subset macrophages
Mac_obj <- subset(CISH_integrated, subset = celltype == "Macrophages")
Idents(Mac_obj) <- "Type"

table(Mac_obj$Type)

# General macrophage DEG
Mac_DEG <- FindMarkers(
  Mac_obj,
  ident.1 = "KO",
  ident.2 = "WT",
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "MAST"
)

Mac_DEG$gene <- rownames(Mac_DEG)

pattern <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

Mac_DEG_filtered <- Mac_DEG %>%
  dplyr::filter(!grepl(pattern, gene))
View(Mac_DEG_filtered)
write.csv(Mac_DEG_filtered, "Mac_DEG_filtered_MAST.csv", row.names = FALSE)


Mac_DEG_sig <- Mac_DEG_filtered %>%
  dplyr::filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25,
    (pct.1 >= 0.10 | pct.2 >= 0.10)
  ) %>%
  dplyr::arrange(p_val_adj)

View(Mac_DEG_sig)

write.csv(Mac_DEG_sig, "Mac_DEG_sig_MAST.csv", row.names = FALSE)

## monocytes 

# Subset monocytes
Mono_obj <- subset(
  CISH_integrated,
  subset = celltype %in% c("Classical monocytes", "Non-classical monocytes")
)

Idents(Mono_obj) <- "Type"

table(Mono_obj$celltype, Mono_obj$Type)

# General monocyte DEG
Mono_DEG <- FindMarkers(
  Mono_obj,
  ident.1 = "KO",
  ident.2 = "WT",
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "MAST"
)

Mono_DEG$gene <- rownames(Mono_DEG)

pattern <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

Mono_DEG_filtered <- Mono_DEG %>%
  dplyr::filter(!grepl(pattern, gene))

View(Mono_DEG_filtered)

Mono_DEG_sig <- Mono_DEG_filtered %>%
  dplyr::filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25,
    (pct.1 >= 0.10 | pct.2 >= 0.10)
  ) %>%
  dplyr::arrange(p_val_adj)

View(Mono_DEG_sig)

write.csv(Mono_DEG, "Mono_DEG_all_MAST.csv", row.names = FALSE)
write.csv(Mono_DEG_filtered, "Mono_DEG_filtered_MAST.csv", row.names = FALSE)
# barplot to visualise all the significant mono and mac genes 

# Subset neutrophils
Neut_obj <- subset(
  CISH_integrated,
  subset = celltype == "Neutrophils"
)

Idents(Neut_obj) <- "Type"

table(Neut_obj$celltype, Neut_obj$Type)

# General neutrophil DEG
Neut_DEG <- FindMarkers(
  Neut_obj,
  ident.1 = "KO",
  ident.2 = "WT",
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "MAST"
)

Neut_DEG$gene <- rownames(Neut_DEG)

pattern <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

Neut_DEG_filtered <- Neut_DEG %>%
  dplyr::filter(!grepl(pattern, gene))

View(Neut_DEG_filtered)
write.csv(Neut_DEG_filtered, "Neut_DEG_filtered_MAST.csv", row.names = FALSE)

Neut_DEG_sig <- Neut_DEG_filtered %>%
  dplyr::filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25,
    (pct.1 >= 0.10 | pct.2 >= 0.10)
  ) %>%
  dplyr::arrange(p_val_adj)

View(Neut_DEG_sig)

##### graph
DefaultAssay(CISH_integrated) <- "RNA"

# Subset macrophages
Mac_obj <- subset(CISH_integrated, subset = celltype == "Macrophages")
Idents(Mac_obj) <- "Type"

table(Mac_obj$Type)

# General macrophage DEG
Mac_DEG <- FindMarkers(
  Mac_obj,
  ident.1 = "KO",
  ident.2 = "WT",
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "MAST"
)

Mac_DEG$gene <- rownames(Mac_DEG)

pattern <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

Mac_DEG_filtered <- Mac_DEG %>%
  dplyr::filter(!grepl(pattern, gene))
View(Mac_DEG_filtered)
write.csv(Mac_DEG_filtered, "Mac_DEG_filtered_MAST.csv", row.names = FALSE)


Mac_DEG_sig <- Mac_DEG_filtered %>%
  dplyr::filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25,
    (pct.1 >= 0.10 | pct.2 >= 0.10)
  ) %>%
  dplyr::arrange(p_val_adj)

View(Mac_DEG_sig)

write.csv(Mac_DEG_sig, "Mac_DEG_sig_MAST.csv", row.names = FALSE)

## monocytes 

# Subset monocytes
Mono_obj <- subset(
  CISH_integrated,
  subset = celltype %in% c("Classical monocytes", "Non-classical monocytes")
)

Idents(Mono_obj) <- "Type"

table(Mono_obj$celltype, Mono_obj$Type)

# General monocyte DEG
Mono_DEG <- FindMarkers(
  Mono_obj,
  ident.1 = "KO",
  ident.2 = "WT",
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "MAST"
)

Mono_DEG$gene <- rownames(Mono_DEG)

pattern <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

Mono_DEG_filtered <- Mono_DEG %>%
  dplyr::filter(!grepl(pattern, gene))

View(Mono_DEG_filtered)

Mono_DEG_sig <- Mono_DEG_filtered %>%
  dplyr::filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25,
    (pct.1 >= 0.10 | pct.2 >= 0.10)
  ) %>%
  dplyr::arrange(p_val_adj)

View(Mono_DEG_sig)

write.csv(Mono_DEG, "Mono_DEG_all_MAST.csv", row.names = FALSE)
write.csv(Mono_DEG_filtered, "Mono_DEG_filtered_MAST.csv", row.names = FALSE)
# barplot to visualise all the significant mono and mac genes 

# Subset neutrophils
Neut_obj <- subset(
  CISH_integrated,
  subset = celltype == "Neutrophils"
)

Idents(Neut_obj) <- "Type"

table(Neut_obj$celltype, Neut_obj$Type)

# General neutrophil DEG
Neut_DEG <- FindMarkers(
  Neut_obj,
  ident.1 = "KO",
  ident.2 = "WT",
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "MAST"
)

Neut_DEG$gene <- rownames(Neut_DEG)

pattern <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

Neut_DEG_filtered <- Neut_DEG %>%
  dplyr::filter(!grepl(pattern, gene))

View(Neut_DEG_filtered)

Neut_DEG_sig <- Neut_DEG_filtered %>%
  dplyr::filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.25,
    (pct.1 >= 0.10 | pct.2 >= 0.10)
  ) %>%
  dplyr::arrange(p_val_adj)

View(Neut_DEG_sig)


##

# Combine significant DEGs from 3 populations
sig_DEG_combined <- dplyr::bind_rows(
  Mac_DEG_sig %>% dplyr::mutate(Celltype = "Macrophages"),
  Mono_DEG_sig %>% dplyr::mutate(Celltype = "Monocytes"),
  Neut_DEG_sig %>% dplyr::mutate(Celltype = "Neutrophils")
) %>%
  dplyr::mutate(Direction = ifelse(avg_log2FC > 0, "Up in KO", "Up in WT")) %>%
  dplyr::arrange(Celltype, avg_log2FC)

# Barplot
sig_DEG_barplot <- ggplot(
  sig_DEG_combined,
  aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = Direction)
) +
  geom_col(color = "black") +
  coord_flip() +
  facet_wrap(~ Celltype, scales = "free_y") +
  scale_fill_manual(values = c("Up in KO" = "lightblue", "Up in WT" = "grey70")) +
  theme_classic(base_size = 12) +
  labs(
    title = "Significant DEGs in Macrophages, Monocytes, and Neutrophils",
    x = "Gene",
    y = "avg_log2FC",
    fill = "Direction"
  ) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

sig_DEG_barplot
############



