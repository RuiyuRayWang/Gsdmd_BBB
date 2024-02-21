setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Seurat)
library(tidyverse)
library(tidydr)

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})
# Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

winkler_adult_brain_endo_mtx <- readr::read_tsv("../data/Winkler_2022/Adult_Brain_Endothelial/exprMatrix.tsv.gz") |> column_to_rownames("gene")
winkler_adult_brain_endo_meta <- readr::read_tsv("../data/Winkler_2022/Adult_Brain_Endothelial/meta.tsv") |> column_to_rownames("Cell")
winkler_adult_brain_endo_umap <- readr::read_tsv("../data/Winkler_2022/Adult_Brain_Endothelial/Seurat_umap.coords.tsv.gz", col_names = c("Cell","UMAP_1","UMAP_2")) |> column_to_rownames("Cell")
winkler_adult_brain_endo <- CreateSeuratObject(counts = winkler_adult_brain_endo_mtx, meta.data = winkler_adult_brain_endo_meta)
winkler_adult_brain_endo <- NormalizeData(winkler_adult_brain_endo)
winkler_adult_brain_endo[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(winkler_adult_brain_endo_umap))

Idents(winkler_adult_brain_endo) <- "clusters"
winkler_adult_brain_endo <- RenameIdents(winkler_adult_brain_endo, `Art1`="Artery1", `Art2`="Artery2", `Art3`="Artery3", `Cap`="Capillary")
d1 <- DimPlot(winkler_adult_brain_endo, label = T, label.box = T, label.size = 5) +
  ggtitle("Winkler et al. (2022):\nadult brain endothelial") + 
  theme_dr(xlength = .2, ylength = .2) +
  theme(
    plot.title = element_text(size = 15),
    legend.position = "none",
    axis.title = element_text(size = 15, hjust = 0)
  )
ggsave(
  filename = "winkler_endo_1.eps", plot = d1, device = "eps", path = "../figures/FigS16/", width = 5, height = 5, dpi = 300, family = "Arial"
)
  
f1_1 <- FeaturePlot(winkler_adult_brain_endo, features = "GSDMD") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.text = element_text(size = 13)
  )
ggsave(
  filename = "winkler_endo_2.eps", plot = f1_1, device = "eps", path = "../figures/FigS16/", width = 5.6, height = 5, dpi = 300, family = "Arial"
)

f1_2 <- FeaturePlot(winkler_adult_brain_endo, features = "CASP4") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.text = element_text(size = 13)
  )
ggsave(
  filename = "winkler_endo_3.eps", plot = f1_2, device = "eps", path = "../figures/FigS16/", width = 5.6, height = 5, dpi = 300, family = "Arial"
)

winkler_endo_arteriovenous_malform_ctrl <- readRDS("../data/Winkler_2022/Endothelial_Arteriovenous_Malformation_Control/final_merged_ec_cb.rds")
Idents(winkler_endo_arteriovenous_malform_ctrl) <- 'Sample'
winkler_endo_arteriovenous_malform_ctrl <- RenameIdents(
  winkler_endo_arteriovenous_malform_ctrl, 
  'AVM049'="AVM", 'AVM064'="AVM", 'AVM211'="AVM", 'AVM420'="AVM", 'AVM524'="AVM",
  'CTRL085'="CTRL", 'CTRL086'="CTRL", 'CTRL099_21'="CTRL", 'CTRL099_22'="CTRL", 'CTRL12'="CTRL"
  )
winkler_endo_arteriovenous_malform_ctrl[["state"]] <- Idents(winkler_endo_arteriovenous_malform_ctrl)
winkler_endo_arteriovenous_malform_ctrl$state <- factor(winkler_endo_arteriovenous_malform_ctrl$state, levels = c("CTRL","AVM"))

# winkler_endo_arteriovenous_malform_ctrl <- NormalizeData(winkler_endo_arteriovenous_malform_ctrl)

VlnPlot(winkler_endo_arteriovenous_malform_ctrl, features = "CASP4", group.by = "state", slot = "counts")
VlnPlot(winkler_endo_arteriovenous_malform_ctrl, features = "CD14", group.by = "state", slot = "counts")
VlnPlot(winkler_endo_arteriovenous_malform_ctrl, features = "GSDMD", group.by = "state", slot = "counts")

d2 <- DimPlot(winkler_endo_arteriovenous_malform_ctrl, group.by = "clusters", label = T, label.box = T) +
  ggtitle("Arteriovenous Malformation & Control") + 
  theme_dr() +
  NoLegend()
f2_1 <- FeaturePlot(winkler_endo_arteriovenous_malform_ctrl, features = "GSDMD") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  )
f2_2 <- FeaturePlot(winkler_endo_arteriovenous_malform_ctrl, features = "CASP4") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  )

d2 | f2_1 | f2_2
