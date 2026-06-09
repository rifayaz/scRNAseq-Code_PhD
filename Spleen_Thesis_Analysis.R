

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



CISHWT_SP_count <- Read10X(data.dir= "U:\\Rifqa\\PhD Documents\\ScRNA_Sequencing\\CISH_scRNAseq\\7Bridges_NewData_Aug\\SPLEEN\\CISHWT_SPLEEN")
CISHKO_SP_count <- Read10X(data.dir= "U:\\Rifqa\\PhD Documents\\ScRNA_Sequencing\\CISH_scRNAseq\\7Bridges_NewData_Aug\\SPLEEN\\CISHKO_SPLEEN")

# create CISHWT_SP and CISHKO_SP as Seurat objects
options(Seurat.object.assay.version = "v3")

CISHWT_SP_Seu <- CreateSeuratObject(CISHWT_SP_count, min.cells = 3, min.features = 150, project = "CISHWT_SP")
class(CISHWT_SP_Seu[["RNA"]])
CISHWT_SP_Seu$Type <- "WT"
CISHWT_SP_Seu <- RenameCells(CISHWT_SP_Seu, add.cell.id = "WT_SP")

CISHKO_SP_Seu <- CreateSeuratObject(CISHKO_SP_count, min.cells = 3, min.features = 150, project = "CISHKO_SP")
class(CISHKO_SP_Seu[["RNA"]])
CISHKO_SP_Seu$Type <- "KO"
CISHKO_SP_Seu <- RenameCells(CISHKO_SP_Seu, add.cell.id = "KO_SP")

sample_cols_SP <- c(
  "CISHWT_SP" = "grey60",
  "CISHKO_SP" = "lightblue",
  "WT" = "grey60",
  "KO" = "lightblue"
)

print(paste0("Cells in CISHWT_SP_Seu: ", ncol(CISHWT_SP_Seu)))
print(paste0("Cells in CISHKO_SP_Seu: ", ncol(CISHKO_SP_Seu)))


# Merge two SP data for first QC
CISH_SP_merged <- merge(CISHWT_SP_Seu, CISHKO_SP_Seu)

CISH_SP_merged[["percent.mt"]] <- PercentageFeatureSet(CISH_SP_merged, pattern = "^mt-")
CISH_SP_merged[["percent.rbc"]] <- PercentageFeatureSet(CISH_SP_merged, pattern = "^Hb[ab]-")

Fig_QC_A_before_SP <- VlnPlot(
  CISH_SP_merged,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbc"),
  group.by = "orig.ident",
  ncol = 4,
  pt.size = 0.1
) & scale_fill_manual(values = sample_cols_SP)

Fig_QC_A_before_SP


# Numbers to help set up cut-off gating
table(CISH_SP_merged$orig.ident)
summary(CISH_SP_merged$nCount_RNA)
summary(CISH_SP_merged$nFeature_RNA)
summary(CISH_SP_merged$percent.mt)
summary(CISH_SP_merged$percent.rbc)

sd(CISH_SP_merged$nCount_RNA)
sd(CISH_SP_merged$nFeature_RNA)
sd(CISH_SP_merged$percent.mt)

p1_SP <- ggplot(CISH_SP_merged@meta.data, mapping = aes(x = nFeature_RNA)) +
  geom_histogram(bins = 120) +
  geom_vline(aes(xintercept = 4000, color = "red")) +
  NoLegend()

p2_SP <- ggplot(CISH_SP_merged@meta.data, mapping = aes(x = percent.mt)) +
  geom_histogram(bins = 120) +
  geom_vline(aes(xintercept = 20, color = "red")) +
  NoLegend()

p3_SP <- ggplot(CISH_SP_merged@meta.data, mapping = aes(x = nCount_RNA)) +
  geom_histogram(bins = 120) +
  geom_vline(aes(xintercept = 15000, color = "red")) +
  NoLegend()

print(p1_SP + p3_SP + p2_SP)


# After applying QC thresholds
CISH_SP_merged.f <- subset(
  CISH_SP_merged,
  subset = percent.mt < 20 &
    percent.rbc < 5 &
    nFeature_RNA < 4000 &
    nFeature_RNA > 300 &
    nCount_RNA < 15000
)

p1_SP <- ggplot(CISH_SP_merged.f@meta.data, mapping = aes(x = nFeature_RNA)) +
  geom_histogram(bins = 120) +
  geom_vline(aes(xintercept = 4000), color = "red") +
  NoLegend()

p2_SP <- ggplot(CISH_SP_merged.f@meta.data, mapping = aes(x = percent.mt)) +
  geom_histogram(bins = 120) +
  geom_vline(aes(xintercept = 20, color = "red")) +
  NoLegend()

p3_SP <- ggplot(CISH_SP_merged.f@meta.data, mapping = aes(x = nCount_RNA)) +
  geom_histogram(bins = 120) +
  geom_vline(aes(xintercept = 15000, color = "red")) +
  NoLegend()

print(p1_SP + p3_SP + p2_SP)


Fig_QC_B_after_SP <- VlnPlot(
  CISH_SP_merged.f,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rbc"),
  group.by = "orig.ident",
  ncol = 4,
  pt.size = 0.1
) & scale_fill_manual(values = sample_cols_SP)

Fig_QC_B_after_SP


# Inspect cell numbers after creating Seurat objects
print(paste0("Cells in CISHWT_SP_Seu: ", ncol(CISHWT_SP_Seu)))
print(paste0("Cells in CISHKO_SP_Seu: ", ncol(CISHKO_SP_Seu)))
print(paste0("Total cells before merge: ", ncol(CISHWT_SP_Seu) + ncol(CISHKO_SP_Seu)))

# Inspect cell numbers after merging
print(paste0("Total cells in CISH_SP_merged: ", ncol(CISH_SP_merged)))
print(table(CISH_SP_merged$orig.ident))
print(table(CISH_SP_merged$Type))

# Inspect cell numbers after initial QC filtering
print(paste0("Total cells after QC filtering: ", ncol(CISH_SP_merged.f)))
print(table(CISH_SP_merged.f$orig.ident))
print(table(CISH_SP_merged.f$Type))

# Inspect number and percent of cells kept after QC
cells_before_qc_SP <- ncol(CISH_SP_merged)
cells_after_qc_SP <- ncol(CISH_SP_merged.f)
cells_removed_qc_SP <- cells_before_qc_SP - cells_after_qc_SP
percent_kept_qc_SP <- round((cells_after_qc_SP / cells_before_qc_SP) * 100, 2)
percent_removed_qc_SP <- round((cells_removed_qc_SP / cells_before_qc_SP) * 100, 2)

print(paste0("Cells before QC: ", cells_before_qc_SP))
print(paste0("Cells after QC: ", cells_after_qc_SP))
print(paste0("Cells removed by QC: ", cells_removed_qc_SP))
print(paste0("Percent kept after QC: ", percent_kept_qc_SP, "%"))
print(paste0("Percent removed after QC: ", percent_removed_qc_SP, "%"))

# Inspect WT/KO cells kept and removed after QC
qc_cell_summary_SP <- data.frame(
  Sample = names(table(CISH_SP_merged$orig.ident)),
  Before_QC = as.numeric(table(CISH_SP_merged$orig.ident)),
  After_QC = as.numeric(table(CISH_SP_merged.f$orig.ident))
)

qc_cell_summary_SP$Removed_QC <- qc_cell_summary_SP$Before_QC - qc_cell_summary_SP$After_QC
qc_cell_summary_SP$Percent_Kept <- round((qc_cell_summary_SP$After_QC / qc_cell_summary_SP$Before_QC) * 100, 2)
qc_cell_summary_SP$Percent_Removed <- round((qc_cell_summary_SP$Removed_QC / qc_cell_summary_SP$Before_QC) * 100, 2)

print(qc_cell_summary_SP)


library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)

qc_plot_SP <- qc_cell_summary_SP %>%
  pivot_longer(
    cols = c(Before_QC, After_QC, Removed_QC, Percent_Kept),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      Before_QC = "Before QC",
      After_QC = "After QC",
      Removed_QC = "Removed by QC",
      Percent_Kept = "% Kept"
    )
  )

facet_labels_counts_SP <- c(
  "Before QC" = "Before QC\nTotal cells before any filtering",
  "After QC" = "After QC\nCells remaining after filtering",
  "Removed by QC" = "Removed by QC\nCells excluded by QC thresholds"
)

facet_labels_pct_SP <- c(
  "% Kept" = "% Kept\nPercentage of cells retained after QC"
)

qc_counts_SP <- qc_plot_SP %>% dplyr::filter(Metric %in% c("Before QC", "After QC", "Removed by QC"))
qc_percent_SP <- qc_plot_SP %>% dplyr::filter(Metric == "% Kept")

p_counts_SP <- ggplot(qc_counts_SP, aes(x = Sample, y = Value, fill = Sample)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Value), vjust = -0.3, size = 3.5) +
  facet_wrap(~Metric, scales = "fixed", nrow = 1, labeller = labeller(Metric = facet_labels_counts_SP)) +
  scale_fill_manual(values = c("CISHWT_SP" = "grey60", "CISHKO_SP" = "lightblue")) +
  coord_cartesian(ylim = c(0, 8000)) +
  labs(title = "SP QC Cell Summary", x = "Sample", y = "Cell Number", fill = "Sample") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, face = "bold"))

p_percent_SP <- ggplot(qc_percent_SP, aes(x = Sample, y = Value, fill = Sample)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(Value, "%")), vjust = -0.3, size = 3.5) +
  facet_wrap(~Metric, scales = "fixed", nrow = 1, labeller = labeller(Metric = facet_labels_pct_SP)) +
  scale_fill_manual(values = c("CISHWT_SP" = "grey60", "CISHKO_SP" = "lightblue")) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Sample", y = "Percent", fill = "Sample") +
  theme_classic() +
  theme(strip.text = element_text(size = 10, face = "bold"))

p_counts_SP | p_percent_SP


## end of initial QC
## Cish expression
## USE RNA here because this is gene expression

DefaultAssay(CISH_SP_merged.f) <- "RNA"
CISH_SP_merged.f <- NormalizeData(CISH_SP_merged.f, verbose = FALSE)

VlnPlot(CISH_SP_merged.f, features = "Cish", group.by = "Type", pt.size = 0.1) + NoLegend()
DotPlot(CISH_SP_merged.f, features = "Cish", group.by = "Type") + RotatedAxis()


# Pre-integration PCA / UMAP on QC-filtered merged SP object
# USE RNA here because this is before integration

DefaultAssay(CISH_SP_merged.f) <- "RNA"

CISH_SP_merged.f <- CISH_SP_merged.f %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(features = VariableFeatures(CISH_SP_merged.f), verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  FindNeighbors(dims = 1:20, verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE)


## Cell-cycle scoring based on known genes
## USE RNA here because cell-cycle scoring uses gene expression

DefaultAssay(CISH_SP_merged.f) <- "RNA"

mmus_s_SP <- gorth(
  cc.genes.updated.2019$s.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name

mmus_g2m_SP <- gorth(
  cc.genes.updated.2019$g2m.genes,
  source_organism = "hsapiens",
  target_organism = "mmusculus"
)$ortholog_name

CISH_SP_merged.f <- CellCycleScoring(
  CISH_SP_merged.f,
  s.features = mmus_s_SP,
  g2m.features = mmus_g2m_SP,
  set.ident = TRUE
)

CISH_SP_merged.f$CC.Difference <- CISH_SP_merged.f$S.Score - CISH_SP_merged.f$G2M.Score

head(CISH_SP_merged.f)

DimPlot(CISH_SP_merged.f, label = TRUE, split.by = "Phase") + NoLegend()

Fig_QC_CellCycle_SP <- DimPlot(CISH_SP_merged.f, reduction = "umap", group.by = "Phase") +
  ggtitle("SP UMAP Colored by Cell-Cycle Phase") +
  theme_classic()

Fig_QC_CellCycle_SP


# Calculate % variance explained by PC1 and PC2

p1_SP <- DimPlot(CISH_SP_merged.f, reduction = "pca", group.by = "Type", shuffle = TRUE)

p2_SP <- DimPlot(CISH_SP_merged.f, reduction = "umap", group.by = "Type", shuffle = TRUE)

p3_SP <- ElbowPlot(CISH_SP_merged.f)

percent_var_SP <- ((CISH_SP_merged.f@reductions$pca@stdev)^2 /
                     sum(rowVars(CISH_SP_merged.f@assays$RNA@scale.data))) * 100

print(percent_var_SP[1:2])
print(p1_SP + p2_SP + p3_SP)

Fig_QC_PCA_Type_SP <- DimPlot(CISH_SP_merged.f, reduction = "pca", group.by = "Type", shuffle = TRUE) +
  scale_color_manual(values = c("WT" = "grey60", "KO" = "lightblue")) +
  ggtitle("SP PCA Colored by Genotype") +
  theme_classic()

Fig_QC_UMAP_Type_SP <- DimPlot(CISH_SP_merged.f, reduction = "umap", group.by = "Type", shuffle = TRUE) +
  scale_color_manual(values = c("WT" = "grey60", "KO" = "lightblue")) +
  ggtitle("SP UMAP Colored by Genotype") +
  theme_classic()

Fig_QC_PCA_Type_SP + Fig_QC_UMAP_Type_SP

print(DimHeatmap(CISH_SP_merged.f, dims = 1:2, nfeatures = 50, balanced = TRUE, reduction = "pca"))

#######if pca plot is good, then proceed

# Integration step for SP
# Split merged SP object for integration
# Protein-coding gene list

hub_SP <- AnnotationHub()
mouse_SP <- query(hub_SP, c("EnsDb", "Mus musculus"))
edb_SP <- hub_SP[["AH113713"]]

keys_SP <- keys(edb_SP, "GENENAME")
columns_SP <- c("GENEID", "ENTREZID", "GENEBIOTYPE")

filter_SP <- ~ gene_name %in% keys_SP & gene_biotype == "protein_coding"
tbl_SP <- ensembldb::select(edb_SP, filter_SP, columns_SP) %>% as_tibble()

mRNA_gene_list_SP <- unique(tbl_SP$GENENAME)


# Integration step
CISH_SP_merged.list <- SplitObject(CISH_SP_merged.f, split.by = "orig.ident")

CISH_SP_merged.list <- lapply(CISH_SP_merged.list, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  return(x)
})

features_SP <- SelectIntegrationFeatures(object.list = CISH_SP_merged.list, nfeatures = 3000)

var_regex_SP <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

features_filtered_SP <- grep(var_regex_SP, features_SP, invert = TRUE, value = TRUE)
features_filtered_SP <- intersect(features_filtered_SP, mRNA_gene_list_SP)

CISH_SP_merged.list <- lapply(CISH_SP_merged.list, function(x) {
  x <- ScaleData(x, features = features_filtered_SP, verbose = FALSE)
  x <- RunPCA(x, features = features_filtered_SP, verbose = FALSE)
  return(x)
})

library(future)

options(future.globals.maxSize = 3000 * 1024^2)  # 3 GB
plan(sequential)

CISH_SP_merged.anchors <- FindIntegrationAnchors(
  object.list = CISH_SP_merged.list,
  anchor.features = features_filtered_SP,
  reduction = "rpca",
  dims = 1:20,
  verbose = FALSE
)

CISH_SP_integrated <- IntegrateData(
  anchorset = CISH_SP_merged.anchors,
  dims = 1:20,
  verbose = FALSE
)
CISH_SP_integrated
Assays(CISH_SP_integrated)
DefaultAssay(CISH_SP_integrated)

# After integration, use integrated assay for clustering / dimensional reduction
DefaultAssay(CISH_SP_integrated) <- "integrated"

CISH_SP_integrated <- ScaleData(CISH_SP_integrated, verbose = FALSE)

CISH_SP_integrated <- RunPCA(
  CISH_SP_integrated,
  npcs = 20,
  features = features_filtered_SP,
  verbose = FALSE
)

CISH_SP_integrated <- RunUMAP(CISH_SP_integrated, dims = 1:20, verbose = FALSE)
CISH_SP_integrated <- FindNeighbors(CISH_SP_integrated, dims = 1:20, verbose = FALSE)
CISH_SP_integrated <- FindClusters(CISH_SP_integrated, resolution = 0.8, verbose = FALSE)


# Resolution testing
# USE integrated here because this tests clustering resolution
DefaultAssay(CISH_SP_integrated) <- "integrated"

object_SP <- Seurat::FindClusters(
  object = CISH_SP_integrated,
  resolution = c(0.1, 0.5, 0.8, 1.0, 1.2)
)

head(object_SP)

clustree(object_SP)

colnames(object_SP@meta.data)
grep("snn_res", colnames(object_SP@meta.data), value = TRUE)

DimPlot(
  object_SP,
  reduction = "umap",
  group.by = c(
    "integrated_snn_res.0.1",
    "integrated_snn_res.0.5",
    "integrated_snn_res.0.8",
    "integrated_snn_res.1",
    "integrated_snn_res.1.2"
  )
)

Idents(CISH_SP_integrated) <- "integrated_snn_res.0.8"


# Check final object
p4_SP <- ElbowPlot(CISH_SP_integrated, n = 20)

CISH_SP_integrated


# Final integrated UMAP plots
# USE integrated here because these are UMAP / clustering plots
DefaultAssay(CISH_SP_integrated) <- "integrated"

p1_SP <- DimPlot(
  CISH_SP_integrated,
  reduction = "umap",
  group.by = "orig.ident",
  shuffle = TRUE
) + ggtitle("SP Mouse type")

p2_SP <- DimPlot(
  CISH_SP_integrated,
  reduction = "umap",
  repel = TRUE,
  shuffle = TRUE
)

p3_SP <- DimPlot(
  CISH_SP_integrated,
  reduction = "umap",
  group.by = "Phase",
  repel = TRUE,
  shuffle = TRUE
)

print(p1_SP + p2_SP + p3_SP)


##############

# Contaminant marker checking for SP
DefaultAssay(CISH_SP_integrated) <- "RNA"

FeaturePlot(CISH_SP_integrated, features = c("Ptprc", "Pecam1", "Cdh5", "Col1a1", "Pdgfra", "Epcam"), ncol = 3)

library(patchwork)

p_feat_SP <- FeaturePlot(CISH_SP_integrated, features = c("Ptprc", "Pecam1", "Cdh5", "Col1a1", 
                                                          "Pdgfra", "Epcam", "Rhd", "Hemgn", 
                                                          "Stmn1", "Hbb-bt"), ncol = 3)

p_feat_SP + plot_annotation(caption = paste("Marker guide:", "Ptprc = pan-immune cells (CD45)", 
                                            "Pecam1 = endothelial cells", 
                                            "Cdh5 = endothelial cells", 
                                            "Col1a1 = fibroblast/stromal cells", 
                                            "Pdgfra = fibroblast/stromal cells", 
                                            "Epcam = epithelial cells", 
                                            "Rhd = erythroid / red blood cell lineage", 
                                            "Hemgn = erythroid / hematopoietic progenitor marker", 
                                            "Stmn1 = proliferating / cycling cells", 
                                            "Hbb-bt = red blood cells / hemoglobin contamination", sep = "\n")) & theme(plot.caption = element_text(size = 10, hjust = 0))
FeaturePlot(CISH_SP_integrated, features = c("Pecam1","Cdh5"))
DefaultAssay(CISH_SP_integrated) <- "integrated"

DimPlot(CISH_SP_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)

#####next

# Marker identification for SP (resolution 0.8)
DefaultAssay(CISH_SP_integrated) <- "RNA"
Idents(CISH_SP_integrated) <- "integrated_snn_res.0.8"

marker_features_SP <- rownames(CISH_SP_integrated)
marker_features_SP <- intersect(marker_features_SP, mRNA_gene_list_SP)

remove_marker_regex_SP <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"
marker_features_SP <- grep(remove_marker_regex_SP, marker_features_SP, invert = TRUE, value = TRUE)

Cluster.markers_SP <- FindAllMarkers(CISH_SP_integrated, features = marker_features_SP, only.pos = TRUE, logfc.threshold = 0.1, test.use = "roc")

head(Cluster.markers_SP)
table(Cluster.markers_SP$cluster)

CISH_SP_gene_Cluster1 <- Cluster.markers_SP %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 20)

View(CISH_SP_gene_Cluster1)
write.csv(CISH_SP_gene_Cluster1, "CISH_SP_gene_Cluster1_res0.8.csv")

FeaturePlot(CISH_SP_integrated, features = c("Cd3e","Nkg7","Klrb1c","Ms4a1","Cd79a","Lyz2"), ncol = 3)


# Extended RBC / erythroid / contamination / cell-cycle checking
DefaultAssay(CISH_SP_integrated) <- "RNA"

FeaturePlot(CISH_SP_integrated,
            features = c(
              "Rhd","Hemgn","Hbb-bt","Hba-a1","Hba-a2",      # RBC
              "Pecam1","Cdh5","Col1a1","Pdgfra","Epcam",     # non-immune contamination
              "Mki67","Top2a","Stmn1","Pcna"                 # cell cycle
            ),
            ncol = 4, order = TRUE
)


# UMAP / clustering display
DefaultAssay(CISH_SP_integrated) <- "integrated"

DimPlot(CISH_SP_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)

#####

######
# remove clusters for SP
# Edit this only after confirming the clusters are unwanted.
Idents(CISH_SP_integrated) <- "integrated_snn_res.0.8"
table(Idents(CISH_SP_integrated))
clusters_to_remove_SP <- c(20)

CISH_SP_integrated <- subset(CISH_SP_integrated, idents = clusters_to_remove_SP, invert = TRUE)

DefaultAssay(CISH_SP_integrated) <- "integrated"
DimPlot(CISH_SP_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)

CISH_SP_integrated <- ScaleData(CISH_SP_integrated, verbose = FALSE)
CISH_SP_integrated <- RunPCA(CISH_SP_integrated, npcs = 20, verbose = FALSE)
CISH_SP_integrated <- RunUMAP(CISH_SP_integrated, dims = 1:20, verbose = FALSE)
CISH_SP_integrated <- FindNeighbors(CISH_SP_integrated, dims = 1:20, verbose = FALSE)
CISH_SP_integrated <- FindClusters(CISH_SP_integrated, resolution = 0.8, verbose = FALSE)

Idents(CISH_SP_integrated) <- "seurat_clusters"

DefaultAssay(CISH_SP_integrated) <- "RNA"
CISH_SP_integrated <- NormalizeData(CISH_SP_integrated, verbose = FALSE)
CISH_SP_integrated <- FindVariableFeatures(CISH_SP_integrated, nfeatures = 3000, verbose = FALSE)

colnames(object_SP@meta.data)
colnames(CISH_SP_integrated@meta.data)

DimPlot(CISH_SP_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)
DimPlot(CISH_SP_integrated, reduction = "umap", split.by = "orig.ident", label = TRUE, repel = TRUE, shuffle = TRUE)

FeaturePlot(CISH_SP_integrated, features = c("Cd3e","Nkg7","Klrb1c"))

####

# Run cluster marker finder again for SP using filtered protein-coding genes only
DefaultAssay(CISH_SP_integrated) <- "RNA"
Idents(CISH_SP_integrated) <- "seurat_clusters"

marker_features_SP <- rownames(CISH_SP_integrated)
marker_features_SP <- intersect(marker_features_SP, mRNA_gene_list_SP)

remove_marker_regex_SP <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"
marker_features_SP <- grep(remove_marker_regex_SP, marker_features_SP, invert = TRUE, value = TRUE)

Cluster.markers_SP <- FindAllMarkers(CISH_SP_integrated, features = marker_features_SP, only.pos = TRUE, logfc.threshold = 0.1, test.use = "roc")

CISH_SP_gene_Cluster <- Cluster.markers_SP %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

View(CISH_SP_gene_Cluster)
DimPlot(CISH_SP_integrated, reduction = "umap", label = TRUE, repel = TRUE, shuffle = TRUE)
FeaturePlot(CISH_SP_integrated, features = c("Cd3e","Nkg7","Klrb1c"))

setwd("C:\\Users\\rfayaz\\OneDrive - Deakin University\\ScRNA CISH\\Thesis analysis\\Spleen_Thesis_2")
write.csv(CISH_SP_gene_Cluster, "CISH_SP_gene_Cluster.csv")

#### next is annotation

DefaultAssay(CISH_SP_integrated) <- "RNA"

#Cell Cluster Annotation
FeaturePlot(CISH_SP_integrated, feature = c ("Cish"))
FeaturePlot(CISH_SP_integrated, feature = c ("Cd3e", "Cd3g","Cd28","Cd5")) #Tcells
FeaturePlot(CISH_SP_integrated, feature = c ("Cd4", "Cd8b1", "Cd8a", "Cd3e", "Il2rb","Il2ra")) #Cd4CD8
FeaturePlot(CISH_SP_integrated, feature = c ("Cd3e")) 
FeaturePlot(CISH_SP_integrated, features = c("Ncr1", "Gzma", "Klrb1c", "Nkg7")) #NK Cells
FeaturePlot(CISH_SP_integrated, feature = c ("Cd4", "Cd3e", "Cd8a"))
FeaturePlot(CISH_SP_integrated, feature = c ("Csf1r","C1qb", "Ms4a6c", "Fcgr1")) #MacsandMonos
FeaturePlot(CISH_SP_integrated, features = c("Ly6c2", "Ccr2", "Cx3cr1", "Ms4a7", "C1qb", "C1qc"))
FeaturePlot(CISH_SP_integrated, feature = c ("C1qc", "C1qb", "Lgmn", "Ms4a7")) #IMs
FeaturePlot(CISH_SP_integrated, feature = c ("Cst3", "Ifitm6","Clec4a3", "Treml4")) #iMons
FeaturePlot(CISH_SP_integrated, feature = c ("Ms4a1","Ighm", "Iglc2", "H2-Eb1")) #Bcells
FeaturePlot(CISH_SP_integrated, feature = c ("S100a9", "Csf3r", "Mmp9", "Retnlg")) #Neutrophils
FeaturePlot(CISH_SP_integrated, feature = c ("Ccr2", "Ly6c2")) #Classical mons
FeaturePlot(CISH_SP_integrated, feature = c ("Cx3cr1", "Cd43", "Spn")) #Non-Classical mons
FeaturePlot(CISH_SP_integrated, feature = c ("Siglech", "Flt3", "Zbtb46", "Rogdi", "Cd209a", "Ccr7")) #Dcs

FeaturePlot(CISH_SP_integrated, feature = c ("Foxp3"))

### Try SingleR
####################### SINGLER WORKFLOW OPTIONAL FOR SP #######################
library(SingleCellExperiment)
CISH_sce_for_SingleR_SP <- as.SingleCellExperiment(CISH_SP_integrated)

library(celldex)
immgen_ref_SP <- ImmGenData()

library(SingleR)
singleR_results_SP <- SingleR(test = CISH_sce_for_SingleR_SP, ref = immgen_ref_SP, labels = immgen_ref_SP$label.main)

colnames(CISH_SP_integrated@meta.data)
head(singleR_results_SP$labels)
table(singleR_results_SP$labels)

CISH_SP_integrated$SingleR_label <- singleR_results_SP$labels
table(CISH_SP_integrated$seurat_clusters, CISH_SP_integrated$SingleR_label)

library(pheatmap)
cluster_label_table_SP <- table(CISH_SP_integrated$seurat_clusters, CISH_SP_integrated$SingleR_label)
cluster_label_prop_SP <- prop.table(cluster_label_table_SP, margin = 1)

my_colors_SP <- colorRampPalette(c("white", "blue", "red"))(100)
pheatmap(cluster_label_prop_SP, color = my_colors_SP, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = TRUE, angle_col = 45, border_color = NA, fontsize_row = 12, fontsize_col = 12, cellwidth = 25, cellheight = 15, main = "SP Cluster vs Cell Type Proportions", legend = TRUE)


## deeper SingleR analysis for SP
library(Seurat)
library(SingleCellExperiment)
library(celldex)
library(SingleR)
library(pheatmap)
library(dplyr)
library(tidyr)

CISH_sce_for_SingleR_SP_fine <- as.SingleCellExperiment(CISH_SP_integrated)
immgen_ref_SP <- ImmGenData()

singleR_results_SP_fine <- SingleR(test = CISH_sce_for_SingleR_SP_fine, ref = immgen_ref_SP, labels = immgen_ref_SP$label.fine)

head(singleR_results_SP_fine$labels)
table(singleR_results_SP_fine$labels)

CISH_SP_integrated$SingleR_label <- singleR_results_SP_fine$labels

cluster_type_counts_SP <- as.data.frame(table(Cluster = CISH_SP_integrated$seurat_clusters, CellType = CISH_SP_integrated$SingleR_label))
cluster_type_counts_SP <- cluster_type_counts_SP[cluster_type_counts_SP$Freq > 0, ]
cluster_type_counts_SP <- cluster_type_counts_SP[order(cluster_type_counts_SP$Cluster, -cluster_type_counts_SP$Freq), ]

cluster_type_counts_SP

cluster_type_matrix_SP <- xtabs(Freq ~ Cluster + CellType, data = cluster_type_counts_SP)
cluster_type_matrix_SP

cluster_totals_SP <- rowSums(cluster_type_matrix_SP)
cluster_summary_SP <- cbind(cluster_type_matrix_SP, TotalCells = cluster_totals_SP)

View(cluster_summary_SP)
write.csv(cluster_summary_SP, "cluster_summary_SP.csv")

### end of singleR

### now to annotate 
DefaultAssay(CISH_SP_integrated) <- "integrated"
Idents(CISH_SP_integrated) <- "seurat_clusters"

# cell annotation and rename each cluster
celltype_SP <- data.frame(ClusterID = 0:22, celltype = "NA")
celltype_SP[celltype_SP$ClusterID %in% c(0,1,2,5,7,8,14,16,19,20), 2] <- "B-cells"
celltype_SP[celltype_SP$ClusterID %in% c(3), 2] <- "CD8 T-cells"
celltype_SP[celltype_SP$ClusterID %in% c(4), 2] <- "CD4 T-cells"
celltype_SP[celltype_SP$ClusterID %in% c(6,21), 2] <- "Neutrophils"
celltype_SP[celltype_SP$ClusterID %in% c(9), 2] <- "T-regs"
celltype_SP[celltype_SP$ClusterID %in% c(10), 2] <- "NK cells"
celltype_SP[celltype_SP$ClusterID %in% c(11), 2] <- "NKT cells"
celltype_SP[celltype_SP$ClusterID %in% c(12), 2] <- "CD8 T-cells"
celltype_SP[celltype_SP$ClusterID %in% c(13), 2] <- "Monocytes"
celltype_SP[celltype_SP$ClusterID %in% c(15), 2] <- "Macrophages"
celltype_SP[celltype_SP$ClusterID %in% c(17), 2] <- "DCs"
celltype_SP[celltype_SP$ClusterID %in% c(18), 2] <- "CD4 T-cells"
celltype_SP[celltype_SP$ClusterID %in% c(22), 2] <- "DCs"

table(CISH_SP_integrated$seurat_clusters)
View(celltype_SP)

table(celltype_SP$celltype)

CISH_SP_integrated@meta.data$celltype <- "NA"
for(i in 1:nrow(celltype_SP)){
  CISH_SP_integrated@meta.data[which(CISH_SP_integrated@meta.data$seurat_clusters == celltype_SP$ClusterID[i]), "celltype"] <- celltype_SP$celltype[i]
}

table(CISH_SP_integrated@meta.data$celltype)
head(CISH_SP_integrated)

DefaultAssay(CISH_SP_integrated) <- "integrated"

DimPlot(CISH_SP_integrated, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("SP Cell annotations")

pcluster1_SP <- DimPlot(CISH_SP_integrated, reduction = "umap", label = TRUE, group.by = "celltype", repel = TRUE, pt.size = 1.3, alpha = 1, label.size = 4.5) + ggtitle("Spleen Immune Cell Types")

CISH_SP_integrated$celltype_wrapped <- stringr::str_wrap(CISH_SP_integrated$celltype, width = 14)

pcluster_SP <- DimPlot(CISH_SP_integrated, reduction = "umap", label = TRUE, group.by = "celltype_wrapped", repel = TRUE, pt.size = 1.3, alpha = 1, label.size = 3.8) + ggtitle("Spleen Immune Cell Types")

pcluster2_SP <- pcluster_SP + theme_dr() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_text(hjust = 0, vjust = -1), axis.title.y = element_text(hjust = 0, vjust = 1, angle = 90), legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold", size = 18), axis.text = element_blank(), axis.title = element_text(size = 15)) + labs(title = "Spleen Immune Cell Types") + guides(fill = "none")

pcluster2_SP

pcluster_legend_only_SP <- DimPlot(CISH_SP_integrated, reduction = "umap", group.by = "celltype_wrapped", label = FALSE, repel = FALSE, pt.size = 1.3, alpha = 1) + theme_dr() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_blank(), axis.title = element_text(size = 15), plot.title = element_text(hjust = 0.5, face = "bold", size = 18), legend.position = "right") + labs(title = "Spleen Immune Cell Types")

pcluster_legend_only_SP

### end of cluster annotation

## FEATURES BY CELL TYPE FOR SP
DefaultAssay(CISH_SP_integrated) <- "RNA"

DotPlot(CISH_SP_integrated, features = c("Foxp3", "Ncr1", "Gzma", "Klrb1c", "Nkg7","Cd3e", "Cd4", "Cd8a", "Ms4a1","S100a9","Il1rn","S100a4","Adgre4","Csf3r","Ms4a7","C1qc","Syngr2","Rogdi","Chil3","Ccr2","Ly6c2","Cx3cr1","Hdc"), group.by = "celltype") +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Selected Features by Cell Type", x = "Features", y = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(angle = 360, hjust = 1, size = 15), axis.title = element_text(size = 12), plot.title = element_text(hjust = 0.5, size = 14), legend.title = element_text(size = 12), legend.text = element_text(size = 12))

VlnPlot(CISH_SP_integrated, features = c("Foxp3", "Ncr1", "Gzma", "Klrb1c", "Nkg7","Cd3e", "Cd4", "Cd8a", "Ms4a1","S100a9","Il1rn","S100a4","Adgre4","Csf3r","Ms4a7","C1qc","Syngr2","Rogdi","Chil3","Ccr2","Ly6c2","Cx3cr1","Hdc"), group.by = "celltype", pt.size = 0, stack = TRUE, flip = TRUE, fill.by = "ident") +
  ggtitle("Selected Feature Expression by Cell Type") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 18), strip.text.y.right = element_text(angle = 20, hjust = 0, size = 18, face = "bold"), strip.background = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")

########################## END OF OPTIONAL ADDITIONAL UMAP STYLES FOR SP ###########################################

# cell percentage
celltype_ratio_SP <- CISH_SP_integrated@meta.data %>%
  group_by(Type, celltype) %>%
  summarise(n = n()) %>%
  mutate(relative_freq = 100 * n / sum(n))

celltype_ratio_SP$celltype <- factor(celltype_ratio_SP$celltype)

custom_colors_SP <- c("WT" = "#D3D3D3", "KO" = "#8dc7cb")

freq_plot_SP <- ggplot(celltype_ratio_SP, aes(x = celltype, y = relative_freq)) +
  geom_col(aes(fill = Type), color = "black", position = "dodge", width = 0.8) +
  ylab("Relative Frequency (%)") +
  scale_fill_manual(values = custom_colors_SP) +
  scale_y_continuous(expand = c(0, 0)) +
  facet_wrap(~ celltype, ncol = 3, scales = "free") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 11), legend.title = element_blank(), legend.text = element_text(size = 12), axis.text.y = element_text(size = 10, color = "black"), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_line(color = "black"))

Frequencyplot_SP <- plot_grid(freq_plot_SP, align = "h", axis = "b", nrow = 1)

Frequencyplot_SP

####
# Total cells per cell type
celltype_counts_SP <- CISH_SP_integrated@meta.data %>%
  group_by(celltype) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  arrange(desc(n_cells))

celltype_counts_SP

celltype_count_plot_SP <- ggplot(celltype_counts_SP, aes(x = reorder(celltype, -n_cells), y = n_cells, fill = celltype)) +
  geom_col(color = "black", width = 0.8) +
  geom_text(aes(label = n_cells), vjust = -0.3, size = 3.5) +
  theme_classic(base_size = 12) +
  labs(title = "Total Cells per Cell Type", x = "Cell Type", y = "Number of Cells") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

celltype_count_plot_SP

###### now for DEG

### DIFFERENTIAL GENE EXPRESSION ANALYSIS FOR SP ######
setwd("C:\\Users\\rfayaz\\OneDrive - Deakin University\\ScRNA CISH\\Thesis analysis\\Spleen_Thesis_2\\DEG")

DefaultAssay(CISH_SP_integrated) <- "RNA"
CISH_SP_integrated <- FindVariableFeatures(CISH_SP_integrated, nfeatures = 3000, verbose = FALSE)
length(VariableFeatures(CISH_SP_integrated))

Idents(CISH_SP_integrated) <- "Type"

Deq_all_wilcox_SP <- FindMarkers(CISH_SP_integrated, ident.1 = "KO", ident.2 = "WT", logfc.threshold = 0, test.use = "wilcox")
Deq_all_MAST_SP <- FindMarkers(CISH_SP_integrated, ident.1 = "KO", ident.2 = "WT", logfc.threshold = 0, test.use = "MAST")
Deq_all_LR_SP <- FindMarkers(CISH_SP_integrated, ident.1 = "KO", ident.2 = "WT", logfc.threshold = 0, test.use = "LR")

View(Deq_all_wilcox_SP)
View(Deq_all_LR_SP)
View(Deq_all_MAST_SP)

Deq_all_MAST_SP$gene <- rownames(Deq_all_MAST_SP)
sig_genes_MAST_SP <- subset(Deq_all_MAST_SP, p_val_adj < 0.05)
View(sig_genes_MAST_SP)

write.csv(Deq_all_MAST_SP, "MAST_SP.csv")
write.csv(Deq_all_LR_SP, "LR_SP.csv")

pattern_SP <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

Deq_all_wilcox_SP$gene <- rownames(Deq_all_wilcox_SP)
Deq_all_wilcox_SP <- Deq_all_wilcox_SP[!grepl(pattern_SP, Deq_all_wilcox_SP$gene), ]
View(Deq_all_wilcox_SP)

sig_genes_wilcox_SP <- subset(Deq_all_wilcox_SP, p_val_adj < 0.05)
View(sig_genes_wilcox_SP)

write.csv(Deq_all_wilcox_SP, "Wilcox_SP.csv")
write.csv(sig_genes_wilcox_SP, "sig_genes_wilcox_SP.csv")

Deq_all_LR_SP$gene <- rownames(Deq_all_LR_SP)
Deq_all_LR_SP <- Deq_all_LR_SP[!grepl(pattern_SP, Deq_all_LR_SP$gene), ]
View(Deq_all_LR_SP)

sig_genes_LR_SP <- subset(Deq_all_LR_SP, p_val_adj < 0.05)
View(sig_genes_LR_SP)

write.csv(Deq_all_LR_SP, "LR_SP.csv")
write.csv(sig_genes_LR_SP, "sig_genes_LR_SP.csv")

Deq_all_MAST_SP$gene <- rownames(Deq_all_MAST_SP)
Deq_all_MAST_SP <- Deq_all_MAST_SP[!grepl(pattern_SP, Deq_all_MAST_SP$gene), ]
View(Deq_all_MAST_SP)

sig_genes_MAST_SP <- subset(Deq_all_MAST_SP, p_val_adj < 0.05)
View(sig_genes_MAST_SP)

write.csv(Deq_all_MAST_SP, "MAST_SP.csv")
write.csv(sig_genes_MAST_SP, "sig_genes_MAST_SP.csv")

sum(sig_genes_wilcox_SP$p_val_adj < 0.05)
sum(sig_genes_MAST_SP$p_val_adj < 0.05)

sum(abs(sig_genes_wilcox_SP$avg_log2FC) > 0.25)
sum(abs(sig_genes_MAST_SP$avg_log2FC) > 0.25)
sum(abs(sig_genes_LR_SP$avg_log2FC) > 0.25)

common_sig_genes_SP <- Reduce(intersect, list(sig_genes_wilcox_SP$gene, sig_genes_MAST_SP$gene, sig_genes_LR_SP$gene))
common_sig_genes_SP
length(common_sig_genes_SP)

length(intersect(sig_genes_wilcox_SP$gene, sig_genes_MAST_SP$gene))
length(intersect(sig_genes_wilcox_SP$gene, sig_genes_LR_SP$gene))
length(intersect(sig_genes_MAST_SP$gene, sig_genes_LR_SP$gene))

length(common_sig_genes_SP) / nrow(sig_genes_wilcox_SP) * 100
length(common_sig_genes_SP) / nrow(sig_genes_MAST_SP) * 100
length(common_sig_genes_SP) / nrow(sig_genes_LR_SP) * 100

overlap_summary_SP <- data.frame(
  Test = c("Wilcox", "MAST", "LR"),
  Total_Significant_Genes = c(nrow(sig_genes_wilcox_SP), nrow(sig_genes_MAST_SP), nrow(sig_genes_LR_SP)),
  Common_With_All_3 = c(length(common_sig_genes_SP), length(common_sig_genes_SP), length(common_sig_genes_SP))
)

overlap_summary_SP$Percent_Common <- round(overlap_summary_SP$Common_With_All_3 / overlap_summary_SP$Total_Significant_Genes * 100, 1)

overlap_summary_SP

library(gridExtra)
overlap_table_plot_SP <- tableGrob(overlap_summary_SP, rows = NULL, theme = ttheme_minimal(base_size = 12, core = list(fg_params = list(fontface = "plain")), colhead = list(fg_params = list(fontface = "bold"))))

grid::grid.newpage()
grid::grid.draw(overlap_table_plot_SP)


####### high confidence gene list plotting

# table and graph of high-confidence SP genes
common_DEG_table_SP <- Deq_all_MAST_SP[Deq_all_MAST_SP$gene %in% common_sig_genes_SP, ]
common_DEG_table_SP <- common_DEG_table_SP[order(common_DEG_table_SP$p_val_adj), ]
common_DEG_table_SP[, c("gene", "avg_log2FC", "p_val_adj")]

View(common_DEG_table_SP)
write.csv(common_DEG_table_SP, "common_DEG_wilcox_MAST_LR_SP.csv")

common_DEG_table_SP_filtered <- common_DEG_table_SP %>%
  dplyr::filter(
    gene %in% common_sig_genes_SP,
    (pct.1 >= 0.1 | pct.2 >= 0.1),
    abs(avg_log2FC) >= 0.25,
    p_val_adj < 0.05
  ) %>%
  dplyr::arrange(p_val_adj, dplyr::desc(abs(avg_log2FC)))

View(common_DEG_table_SP_filtered)
write.csv(common_DEG_table_SP_filtered, "common_DEG_SP_filtered_pct0.1_logFC0.25_padj0.05.csv")

DefaultAssay(CISH_SP_integrated) <- "RNA"

CISH_SP_integrated <- ScaleData(CISH_SP_integrated, features = common_DEG_table_SP_filtered, verbose = FALSE)

#top 10 upregulated and dowreg only
top10_up_SP <- common_DEG_table_SP_filtered %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 10)

top10_up_genes_SP <- top10_up_SP$gene

top10_down_SP <- common_DEG_table_SP_filtered %>%
  dplyr::filter(avg_log2FC < 0) %>%
  dplyr::arrange(avg_log2FC) %>%
  dplyr::slice_head(n = 10)

top10_down_genes_SP <- top10_down_SP$gene

######
library(ggplot2)

top20_DEG_SP <- rbind(top10_up_SP, top10_down_SP)

top20_DEG_SP$Direction <- ifelse(top20_DEG_SP$avg_log2FC > 0, "Up in KO", "Down in KO")

ggplot(top20_DEG_SP, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = Direction)) +
  geom_col(color = "black") +
  coord_flip() +
  scale_fill_manual(values = c("Up in KO" = "#8dc7cb", "Down in KO" = "#D3D3D3")) +
  theme_classic(base_size = 12) +
  labs(title = "Core DEG signature Genes Common Across All Statistical Tests", x = "Gene", y = "avg_log2FC", fill = "Direction")

###### next, set thresholds for up and down for further analysis
# using MAST from here onwards for SP
##############################################################################################

summary_positive_SP <- summary(Deq_all_MAST_SP$avg_log2FC[Deq_all_MAST_SP$avg_log2FC > 0])
summary_negative_SP <- summary(Deq_all_MAST_SP$avg_log2FC[Deq_all_MAST_SP$avg_log2FC < 0])

first_positive_SP <- summary_positive_SP[2]
mean_positive_SP <- mean(Deq_all_MAST_SP$avg_log2FC[Deq_all_MAST_SP$avg_log2FC > 0])
third_positive_SP <- summary_positive_SP[5]

first_negative_SP <- summary_negative_SP[2]
mean_negative_SP <- mean(Deq_all_MAST_SP$avg_log2FC[Deq_all_MAST_SP$avg_log2FC < 0])
third_negative_SP <- summary_negative_SP[5]

cat("Summary for SP FC > 0:\n")
cat("1st Quartile:", first_positive_SP, "\n")
cat("Mean:", mean_positive_SP, "\n")
cat("3rd Quartile:", third_positive_SP, "\n\n")

cat("Summary for SP FC < 0:\n")
cat("1st Quartile:", first_negative_SP, "\n")
cat("Mean:", mean_negative_SP, "\n")
cat("3rd Quartile:", third_negative_SP, "\n")

Deq_all_MAST_SP$group <- case_when(
  Deq_all_MAST_SP$p_val_adj < 0.05 & Deq_all_MAST_SP$avg_log2FC > third_positive_SP ~ "up",
  Deq_all_MAST_SP$p_val_adj < 0.05 & Deq_all_MAST_SP$avg_log2FC < first_negative_SP ~ "down",
  TRUE ~ "none"
)

Deq_all_MAST_SP$Sig <- case_when(
  Deq_all_MAST_SP$p_val_adj < 0.05 ~ "Sig",
  TRUE ~ "NS"
)

table(Deq_all_MAST_SP$group)

#####

# Create a mapping table of gene symbols to Entrez IDs for SP
mapped_genes_SP <- data.frame(GeneName = rownames(Deq_all_MAST_SP),
                              ensemblID = mapIds(org.Mm.eg.db, keys = Deq_all_MAST_SP$gene, keytype = "SYMBOL", column = "ENTREZID"))

Deq_all_MAST_SP$ensemblID <- mapped_genes_SP$ensemblID

View(Deq_all_MAST_SP)
write.csv(Deq_all_MAST_SP, "DGE_MAST_CISHSP_EntrezID.csv")

Deq_all_MAST_SP$p_val_adj <- ifelse(Deq_all_MAST_SP$p_val_adj == 0, 1e-304, Deq_all_MAST_SP$p_val_adj)

Deq_all_MAST_top_SP <- Deq_all_MAST_SP %>%
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.1, pct.1 > 0.10 | pct.2 > 0.10)

View(Deq_all_MAST_top_SP)
write.csv(Deq_all_MAST_top_SP, "Deq_all_MAST_top_SP_SignificantOnly.csv", row.names = TRUE)

get_top_DEGs_SP <- function(de_table, n = 20, direction = c("up", "down")) {
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

top_10_DEG_up_SP <- get_top_DEGs_SP(Deq_all_MAST_top_SP, 10, "up")
top_20_DEG_up_SP <- get_top_DEGs_SP(Deq_all_MAST_top_SP, 20, "up")
top_30_DEG_up_SP <- get_top_DEGs_SP(Deq_all_MAST_top_SP, 30, "up")

top_10_DEG_down_SP <- get_top_DEGs_SP(Deq_all_MAST_top_SP, 10, "down")
top_20_DEG_down_SP <- get_top_DEGs_SP(Deq_all_MAST_top_SP, 20, "down")
top_30_DEG_down_SP <- get_top_DEGs_SP(Deq_all_MAST_top_SP, 30, "down")

top_10_DEG_SP <- c(top_10_DEG_up_SP, top_10_DEG_down_SP)
top_20_DEG_SP <- c(top_20_DEG_up_SP, top_20_DEG_down_SP)
top_30_DEG_SP <- c(top_30_DEG_up_SP, top_30_DEG_down_SP)

### next is volcano plot and creating DEG tables
##########################
# Create table with log2FC and direction for SP
##########################
create_DEG_table_SP <- function(genes, de_table) {
  data.frame(
    Gene = genes,
    log2FC = de_table[genes, "avg_log2FC"],
    p_val_adj = de_table[genes, "p_val_adj"],
    Direction = ifelse(de_table[genes, "avg_log2FC"] > 0, "Up in KO", "Up in WT")
  )
}

top_20_table_SP <- create_DEG_table_SP(top_20_DEG_SP, Deq_all_MAST_top_SP)
View(top_20_table_SP)

##########################
# Volcano Plot for SP
##########################
keyvals_SP <- ifelse(Deq_all_MAST_top_SP$group == "up", "red",
                     ifelse(Deq_all_MAST_top_SP$group == "down", "royalblue", "grey80"))

names(keyvals_SP)[keyvals_SP == "red"] <- "Up in KO"
names(keyvals_SP)[keyvals_SP == "royalblue"] <- "Up in WT"
names(keyvals_SP)[keyvals_SP == "grey80"] <- "NS"

Deq_plot_SP <- EnhancedVolcano(
  Deq_all_MAST_top_SP,
  x = "avg_log2FC",
  y = "p_val_adj",
  colCustom = keyvals_SP,
  lab = rownames(Deq_all_MAST_top_SP),
  selectLab = top_30_DEG_SP,
  drawConnectors = TRUE,
  pCutoff = 0.05,
  FCcutoff = 0.25,
  max.overlaps = 60,
  title = "Spleen KO vs WT"
)

Deq_plot_SP

##########################
# Heatmap for SP
##########################
heatmap_matrix_SP <- data.frame(log2FC = Deq_all_MAST_top_SP[top_20_DEG_SP, "avg_log2FC"])
rownames(heatmap_matrix_SP) <- top_20_DEG_SP
heatmap_matrix_SP <- as.matrix(heatmap_matrix_SP)

breaks_SP <- seq(min(heatmap_matrix_SP), max(heatmap_matrix_SP), length.out = 40)
colors_SP <- colorRampPalette(c("blue", "white", "red"))(length(breaks_SP) - 1)

pheatmap(
  heatmap_matrix_SP,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colors_SP,
  show_rownames = TRUE,
  show_colnames = FALSE,
  main = "Top SP Differentially Expressed Genes",
  legend = TRUE
)

table(Deq_all_MAST_top_SP$group)

#####next, we prepare for gene ontology#######
###############################
# GO + GSEA workflow for SP
###############################

nrow(Deq_all_MAST_top_SP)

sum(Deq_all_MAST_top_SP$p_val_adj < 0.05)
sum(Deq_all_MAST_top_SP$p_val_adj < 0.05 & abs(Deq_all_MAST_top_SP$avg_log2FC) > 0.25)
sum(Deq_all_MAST_top_SP$p_val_adj < 0.05 & abs(Deq_all_MAST_top_SP$avg_log2FC) > 0.25 & Deq_all_MAST_top_SP$pct.1 > 0.10)

names(CISH_SP_integrated@assays)
DefaultAssay(CISH_SP_integrated)
DefaultAssay(CISH_SP_integrated) <- "RNA"

Idents(CISH_SP_integrated) <- "Type"

###############################
# Prepare DEG table
###############################

if (!exists("Deq_all_MAST_top_SP")) {
  Deq_all_MAST_top_SP <- Deq_all_MAST_SP
}

if (!"gene" %in% colnames(Deq_all_MAST_top_SP)) {
  Deq_all_MAST_top_SP$gene <- rownames(Deq_all_MAST_top_SP)
}

remove_gene_pattern_SP <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

Deq_all_MAST_top_SP <- Deq_all_MAST_top_SP %>%
  dplyr::filter(!grepl(remove_gene_pattern_SP, gene))

Deq_all_MAST_top_SP$p_val_adj <- ifelse(Deq_all_MAST_top_SP$p_val_adj == 0, 1e-304, Deq_all_MAST_top_SP$p_val_adj)

###############################
# Define significant genes
###############################

sig_genes_SP <- Deq_all_MAST_top_SP %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  pull(gene)

length(sig_genes_SP)

###############################
# Map SYMBOL to ENSEMBL
###############################

mapped_SP <- bitr(sig_genes_SP, fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Mm.eg.db)
mapped_SP <- mapped_SP %>% distinct(SYMBOL, .keep_all = TRUE)

head(mapped_SP)

FeaturePlot(CISH_SP_integrated, features = c("Lrrc8b","Shb"), blend = TRUE, order = TRUE)

###############################
# GO enrichment
###############################

GOBP_all_SP <- enrichGO(gene = mapped_SP$ENSEMBL, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)

GOMF_all_SP <- enrichGO(gene = mapped_SP$ENSEMBL, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)

GOCC_all_SP <- enrichGO(gene = mapped_SP$ENSEMBL, OrgDb = org.Mm.eg.db, keyType = "ENSEMBL", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", readable = TRUE)

###############################
# GO dotplots
###############################

dotplot(GOBP_all_SP, showCategory = 10) + ggtitle("GO Biological Process Enrichment in Spleen")
dotplot(GOMF_all_SP, showCategory = 10) + ggtitle("GO Molecular Function Enrichment in Spleen")
dotplot(GOCC_all_SP, showCategory = 10) + ggtitle("GO Cellular Component Enrichment in Spleen")

###############################
# Prepare log2FC table
###############################

log2FC_df_SP <- Deq_all_MAST_top_SP %>%
  dplyr::filter(gene %in% sig_genes_SP) %>%
  dplyr::select(gene, avg_log2FC) %>%
  left_join(mapped_SP, by = c("gene" = "SYMBOL"))

###############################
# Function for GO heatmap
###############################

make_GO_heatmap_matrix_SP <- function(go_object, log2FC_df, top_n = 5) {
  go_df <- as.data.frame(go_object)
  if (nrow(go_df) == 0) return(NULL)
  
  go_df <- go_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(mean_log2FC = mean(log2FC_df$avg_log2FC[log2FC_df$gene %in% strsplit(geneID, "/")[[1]]], na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(mean_log2FC))
  
  combined <- dplyr::bind_rows(
    go_df %>% dplyr::arrange(dplyr::desc(mean_log2FC)) %>% dplyr::slice_head(n = top_n),
    go_df %>% dplyr::arrange(mean_log2FC) %>% dplyr::slice_head(n = top_n)
  ) %>%
    dplyr::distinct(Description, .keep_all = TRUE)
  
  matrix(combined$mean_log2FC, nrow = nrow(combined), ncol = 1, dimnames = list(combined$Description, ""))
}

###############################
# Create heatmaps
###############################

heatmap_matrix_BP_SP <- make_GO_heatmap_matrix_SP(GOBP_all_SP, log2FC_df_SP, top_n = 10)
heatmap_matrix_MF_SP <- make_GO_heatmap_matrix_SP(GOMF_all_SP, log2FC_df_SP, top_n = 10)
heatmap_matrix_CC_SP <- make_GO_heatmap_matrix_SP(GOCC_all_SP, log2FC_df_SP, top_n = 10)
wrap_row_names <- function(names, width = 40) sapply(names, function(x) paste(strwrap(x, width = width), collapse = "\n"))

all_values_SP <- c(as.vector(heatmap_matrix_BP_SP), as.vector(heatmap_matrix_MF_SP), as.vector(heatmap_matrix_CC_SP))

col_fun_SP <- colorRamp2(c(min(all_values_SP, na.rm = TRUE), 0, max(all_values_SP, na.rm = TRUE)), c("blue","white","red"))

rownames(heatmap_matrix_BP_SP) <- wrap_row_names(rownames(heatmap_matrix_BP_SP))
rownames(heatmap_matrix_MF_SP) <- wrap_row_names(rownames(heatmap_matrix_MF_SP))
rownames(heatmap_matrix_CC_SP) <- wrap_row_names(rownames(heatmap_matrix_CC_SP))

p_BP_SP <- Heatmap(heatmap_matrix_BP_SP, name = "mean_log2FC", col = col_fun_SP, show_row_names = TRUE, show_column_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, row_title = "Biological Process")

p_MF_SP <- Heatmap(heatmap_matrix_MF_SP, name = "mean_log2FC", col = col_fun_SP, show_row_names = TRUE, show_column_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, row_title = "Molecular Function")

p_CC_SP <- Heatmap(heatmap_matrix_CC_SP, name = "mean_log2FC", col = col_fun_SP, show_row_names = TRUE, show_column_names = FALSE, cluster_rows = FALSE, cluster_columns = FALSE, row_title = "Cellular Component")

draw(p_BP_SP %v% p_MF_SP %v% p_CC_SP, heatmap_legend_side = "left")

library(openxlsx)

GO_sig_BP_SP <- as.data.frame(GOBP_all_SP) %>% dplyr::filter(p.adjust < 0.05)
GO_sig_MF_SP <- as.data.frame(GOMF_all_SP) %>% dplyr::filter(p.adjust < 0.05)
GO_sig_CC_SP <- as.data.frame(GOCC_all_SP) %>% dplyr::filter(p.adjust < 0.05)

write.xlsx(
  list(
    Biological_Process = GO_sig_BP_SP,
    Molecular_Function = GO_sig_MF_SP,
    Cellular_Component = GO_sig_CC_SP
  ),
  file = "GO_significant_terms_SP.xlsx",
  rowNames = FALSE
)
getwd()
##### next is GSEA on filtered genes 
#######
## GSEA using SP MAST DEG results
## Strictly filtered by pct.1/pct.2 and abs(avg_log2FC) > 0.25

DefaultAssay(CISH_SP_integrated) <- "RNA"

Deq_all_MAST_pct_filtered_SP <- Deq_all_MAST_SP %>%
  dplyr::filter((pct.1 >= 0.10 | pct.2 >= 0.10) & abs(avg_log2FC) > 0.25)

View(Deq_all_MAST_pct_filtered_SP)

gsea_input_SP <- Deq_all_MAST_pct_filtered_SP
gsea_input_SP$gene <- rownames(gsea_input_SP)

remove_gene_pattern_SP <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

gsea_input_SP <- gsea_input_SP %>%
  dplyr::filter(!grepl(remove_gene_pattern_SP, gene))

gene_list_GSEA_SP <- gsea_input_SP$avg_log2FC
names(gene_list_GSEA_SP) <- gsea_input_SP$gene

gene_list_GSEA_SP <- gene_list_GSEA_SP[!is.na(gene_list_GSEA_SP)]
gene_list_GSEA_SP <- gene_list_GSEA_SP[!is.na(names(gene_list_GSEA_SP))]
gene_list_GSEA_SP <- gene_list_GSEA_SP[!duplicated(names(gene_list_GSEA_SP))]
gene_list_GSEA_SP <- sort(gene_list_GSEA_SP, decreasing = TRUE)

mapping_GSEA_SP <- bitr(names(gene_list_GSEA_SP), fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Mm.eg.db) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

gene_list_ens_SP <- gene_list_GSEA_SP[mapping_GSEA_SP$SYMBOL]
names(gene_list_ens_SP) <- mapping_GSEA_SP$ENSEMBL

gene_list_ens_SP <- gene_list_ens_SP[!is.na(names(gene_list_ens_SP))]
gene_list_ens_SP <- gene_list_ens_SP[!duplicated(names(gene_list_ens_SP))]
gene_list_ens_SP <- sort(gene_list_ens_SP, decreasing = TRUE)

gsea_BP_SP <- gseGO(geneList = gene_list_ens_SP, OrgDb = org.Mm.eg.db, ont = "BP", keyType = "ENSEMBL", minGSSize = 15, maxGSSize = 300, pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0, verbose = FALSE)

gsea_BP_sig_SP <- as.data.frame(gsea_BP_SP) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust)

View(gsea_BP_sig_SP)

dotplot(gsea_BP_SP, showCategory = head(gsea_BP_sig_SP$Description, 10)) +
  ggtitle("SP GSEA GO BP: KO vs WT")

gsea_MF_SP <- gseGO(geneList = gene_list_ens_SP, OrgDb = org.Mm.eg.db, ont = "MF", keyType = "ENSEMBL", minGSSize = 15, maxGSSize = 300, pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0, verbose = FALSE)

gsea_MF_sig_SP <- as.data.frame(gsea_MF_SP) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust)

View(gsea_MF_sig_SP)

dotplot(gsea_MF_SP, showCategory = head(gsea_MF_sig_SP$Description, 10)) +
  ggtitle("SP GSEA GO MF: KO vs WT")

gsea_CC_SP <- gseGO(geneList = gene_list_ens_SP, OrgDb = org.Mm.eg.db, ont = "CC", keyType = "ENSEMBL", minGSSize = 15, maxGSSize = 300, pvalueCutoff = 1, pAdjustMethod = "BH", eps = 0, verbose = FALSE)

gsea_CC_sig_SP <- as.data.frame(gsea_CC_SP) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust)

View(gsea_CC_sig_SP)

dotplot(gsea_CC_SP, showCategory = head(gsea_CC_sig_SP$Description, 10)) +
  ggtitle("SP GSEA GO CC: KO vs WT")

# Save GO GSEA significant results for SP

write.csv(gsea_BP_sig_SP, "GSEA_GO_BP_sig_SP.csv", row.names = FALSE)
write.csv(gsea_MF_sig_SP, "GSEA_GO_MF_sig_SP.csv", row.names = FALSE)
write.csv(gsea_CC_sig_SP, "GSEA_GO_CC_sig_SP.csv", row.names = FALSE)
#####
#### Hallmark GSEA for SP gave nothing 
### try kegg
#### Simple KEGG enrichment for SP

library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)

kegg_genes_SP <- Deq_all_MAST_top_SP %>%
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25, pct.1 >= 0.10 | pct.2 >= 0.10) %>%
  dplyr::pull(gene)

kegg_mapped_SP <- bitr(kegg_genes_SP, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

kegg_SP <- enrichKEGG(
  gene = kegg_mapped_SP$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

kegg_SP <- setReadable(kegg_SP, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

kegg_sig_SP <- as.data.frame(kegg_SP) %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(p.adjust)

View(kegg_sig_SP)

dotplot(kegg_SP, showCategory = 10) +
  ggtitle("Spleen KEGG Pathway Enrichment")

write.csv(kegg_sig_SP, "KEGG_sig_SP.csv", row.names = FALSE)

### next, we can subset and see whats going on
# next, subsetting for NK cells spleen SP
setwd("C:\\Users\\rfayaz\\OneDrive - Deakin University\\ScRNA CISH\\Thesis analysis\\Spleen_Thesis_2\\DEG\\Subset")
DefaultAssay(CISH_SP_integrated) <- "RNA"

NK_obj_SP <- subset(CISH_SP_integrated, subset = celltype == "NK cells")
Idents(NK_obj_SP) <- "Type"

NK_DEG_SP <- FindMarkers(NK_obj_SP, ident.1 = "KO", ident.2 = "WT", logfc.threshold = 0, min.pct = 0.05, test.use = "MAST")

NK_DEG_SP$gene <- rownames(NK_DEG_SP)

pattern_SP <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

NK_DEG_filtered_SP <- NK_DEG_SP %>%
  dplyr::filter(!grepl(pattern_SP, gene))

NK_DEG_sig_SP <- NK_DEG_filtered_SP %>%
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25, pct.1 > 0.10 | pct.2 > 0.10) %>%
  dplyr::arrange(p_val_adj)

View(NK_DEG_filtered_SP)
View(NK_DEG_sig_SP)

## no significant DEGs so lets do an exploratory box plot

FeaturePlot(NK_obj_SP, features = c("Nkg7","Fasl","Klrd1","Klrk1","Gzma","Gzmb","Prf1"), split.by = "Type", order = TRUE)
VlnPlot(NK_obj_SP, features = c("Nkg7","Fasl","Klrd1","Klrk1","Gzma","Gzmb","Prf1"), group.by = "Type", pt.size = 0)
DotPlot(NK_obj_SP, features = c("Nkg7","Fasl","Klrd1","Klrk1","Gzma","Gzmb","Prf1"), group.by = "Type") + RotatedAxis()

write.csv(NK_DEG_filtered_SP, "NK_DEG_all_SP.csv", row.names = TRUE)
write.csv(NK_DEG_sig_SP, "NK_DEG_sig_SP.csv", row.names = TRUE)

### subset DCs
# subset DCs and run DEG for SP

DefaultAssay(CISH_SP_integrated) <- "RNA"

DC_obj_SP <- subset(CISH_SP_integrated, subset = celltype == "DCs")
Idents(DC_obj_SP) <- "Type"

table(DC_obj_SP$Type)

DC_DEG_SP <- FindMarkers(
  DC_obj_SP,
  ident.1 = "KO",
  ident.2 = "WT",
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "MAST"
)

DC_DEG_SP$gene <- rownames(DC_DEG_SP)

pattern_SP <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

DC_DEG_filtered_SP <- DC_DEG_SP %>%
  dplyr::filter(!grepl(pattern_SP, gene))

DC_DEG_sig_SP <- DC_DEG_filtered_SP %>%
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25, pct.1 > 0.10 | pct.2 > 0.10) %>%
  dplyr::arrange(p_val_adj)

View(DC_DEG_filtered_SP)
View(DC_DEG_sig_SP)

write.csv(DC_DEG_filtered_SP, "DC_DEG_all_SP.csv", row.names = TRUE)
write.csv(DC_DEG_sig_SP, "DC_DEG_sig_SP.csv", row.names = TRUE)

##B-cells

# subset B-cells and run DEG for SP

DefaultAssay(CISH_SP_integrated) <- "RNA"

Bcell_obj_SP <- subset(CISH_SP_integrated, subset = celltype == "B-cells")
Idents(Bcell_obj_SP) <- "Type"

table(Bcell_obj_SP$Type)

Bcell_DEG_SP <- FindMarkers(
  Bcell_obj_SP,
  ident.1 = "KO",
  ident.2 = "WT",
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "MAST"
)

Bcell_DEG_SP$gene <- rownames(Bcell_DEG_SP)

pattern_SP <- "^mt-|^Rpl\\d+|^Rps\\d+|^Hbb|^Hba|^Gm|Rik$|^ENSMUSG"

Bcell_DEG_filtered_SP <- Bcell_DEG_SP %>%
  dplyr::filter(!grepl(pattern_SP, gene))

Bcell_DEG_sig_SP <- Bcell_DEG_filtered_SP %>%
  dplyr::filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25, pct.1 > 0.10 | pct.2 > 0.10) %>%
  dplyr::arrange(p_val_adj)

View(Bcell_DEG_filtered_SP)
View(Bcell_DEG_sig_SP)

write.csv(Bcell_DEG_filtered_SP, "Bcell_DEG_all_SP.csv", row.names = TRUE)
write.csv(Bcell_DEG_sig_SP, "Bcell_DEG_sig_SP.csv", row.names = TRUE)

#### 
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("BiocGenerics", "BiocNeighbors", "BiocParallel", "ComplexHeatmap"))
install.packages(c("NMF", "circlize", "ggalluvial", "patchwork", "igraph", "future", "future.apply"))
remotes::install_github("sqjin/CellChat")

install.packages("remotes")
remotes::install_github("sqjin/CellChat")
library(CellChat)























