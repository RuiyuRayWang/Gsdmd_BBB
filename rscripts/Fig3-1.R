setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Seurat)
library(ggplot2)
library(tidydr)
library(RColorBrewer)


suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})

if(!dir.exists("../figures/Fig3")) {dir.create("../figures/Fig3")}

cells.combined <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/cells_use_combined.h5Seurat")

# Fig3a
cells.combined <- subset(cells.combined, subset = sample %in% c("WT PBS", "WT LPS"))
cells.combined$cell_subtype <- factor(
  cells.combined$cell_subtype, 
  levels = c("Endo","Astro","OPCs","COPs","NFOLs","MOL1","MOL2","Micro",
             "Ex-L2/3 IT","Ex-L4/5 IT","Ex-L5 IT","Ex-L5/6 NP1","Ex-L5/6 NP2","Ex-L6 IT","Ex-L6 CT","Ex-L6b","Inh-Chodl","Inh-Lamp5","Inh-Lhx6","Inh-Vip")
)
Idents(cells.combined) <- "cell_subtype"
cells.combined <- RenameIdents(cells.combined, 
                               Endo="1   Endo", Astro="2   Astro", OPCs="3   OPCs", 
                               COPs="4   COPs", NFOLs="5   NFOLs", MOL1="6   MOL1", MOL2="7   MOL2", Micro="8   Micro",
                               `Ex-L2/3 IT`="9   Ex-L2/3 IT", `Ex-L4/5 IT`="10 Ex-L4/5 IT", `Ex-L5 IT`="11 Ex-L5 IT", `Ex-L5/6 NP1`="12 Ex-L5/6 NP1", 
                               `Ex-L5/6 NP2`="13 Ex-L5/6 NP2", `Ex-L6 IT`="14 Ex-L6 IT", `Ex-L6 CT`="15 Ex-L6 CT", `Ex-L6b`="16 Ex-L6b", 
                               `Inh-Chodl`="17 Inh-Chodl", `Inh-Lamp5`="18 Inh-Lamp5", `Inh-Lhx6`="19 Inh-Lhx6", `Inh-Vip`="20 Inh-Vip")
DimPlot(cells.combined, 
        # group.by = "cell_subtype",
        reduction = "umap.int") +
  ggtitle(NULL) +
  theme_dr(xlength = .2, ylength = .2) +
  theme(
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    axis.title = element_text(hjust = .02, size = 14)
  )
ggsave(filename = "umap_int_cell_subtype.eps", device = "eps", path = "../figures/Fig3/", width = 6.8, height = 5, dpi = 300, family = "Arial")
