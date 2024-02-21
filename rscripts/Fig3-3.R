setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Augur)
library(Seurat)
library(tidyverse)
library(viridis)
source("plot_umap_refactored.R")

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})

augur_cellsubtype <- readRDS("../outs/Augur/augur_cellsubtype.rds")

# ======================================================================================================
# Lollipop
# plot_lollipop(augur_cellsubtype[["WT"]])  ## hack into the function to produce colored plot

aucs = augur_cellsubtype[["WT"]]$AUC
aucs = aucs |>
  mutate(group = recode(
    cell_type, `Endo` = "Vasculature",
    `Micro` = "Glia", `Astro` = "Glia", `MOL2` = "Glia", `MOL1` = "Glia", `OPCs` = "Glia",
    `Ex-L6 IT` = "Neuron", `Ex-L2/3 IT` = "Neuron", `Ex-L5 IT` = "Neuron", `Ex-L6b` = "Neuron", `Ex-L6 CT` = "Neuron", `Ex-L4/5 IT` = "Neuron", 
    `Ex-L5/6 NP1` = "Neuron", `Ex-L5/6 NP2` = "Neuron",
    `Inh-Lamp5` = "Neuron", `Inh-Vip` = "Neuron", `Inh-Lhx6` = "Neuron"))
size_sm = 14
size_lg = 16
range = range(aucs$auc)
expand = abs(diff(range)) * 0.2
p = aucs %>% ggplot(aes(x = reorder(cell_type, auc), y = auc)) + 
  geom_hline(aes(yintercept = 0.5), linetype = "dotted", size = 0.3) + 
  geom_point(aes(color = group), size = 1.5) + 
  geom_text(aes(label = format(auc, digits = 3), y = ifelse(auc < 0.5, 0.5, auc)), size = 4, nudge_y = expand, hjust = 0.5, color = "black") + 
  geom_segment(aes(color = group, xend = cell_type, yend = 0.5)) + scale_y_continuous("Augur AUC", limits = c(min(range[1] - expand, 0.5), range[2] + expand * 1.5)) + 
  scale_color_manual(values = brewer.pal(12,"Paired")[c(5,7,3)], breaks = c("Vasculature", "Glia", "Neuron")) +
  coord_flip() + 
  theme_bw() + 
  theme(axis.title.x = element_text(size = size_lg), axis.text.x = element_text(size = size_sm, color = "black"), 
        axis.title.y = element_blank(), axis.text.y = element_text(size = size_sm-0.6, color = "black"),
        panel.grid = element_blank(), 
        strip.text = element_text(size = size_lg), strip.background = element_blank(), 
        axis.line.y = element_blank(), axis.line.x = element_blank(), 
        legend.position = "none", legend.text = element_text(size = size_sm), 
        legend.title = element_blank(), legend.key.size = unit(0.6, "lines"), 
        legend.margin = margin(rep(0, 4)), legend.background = element_blank(), 
        plot.title = element_text(size = size_lg, hjust = 0.5))
ggsave(filename = "wt_lollipop_update.eps", plot = p, device = "eps", path = "../figures/Fig3/", width = 4, height = 4, dpi = 300, family = "Arial")

# ======================================================================================================
# Augur AUC on UMAP
cells.combined <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/cells_use_combined.h5Seurat")
sc_plot = cells.combined |> subset(genotype == "WT")

u <- plot_umap_refactored(
  augur = augur_cellsubtype[["WT"]], sc = sc_plot, palette = "OrRd", reduction = "umap.int", cell_type_col = "cell_subtype", 
  size_sm = 14, size_lg = 16, lgd.pos = c(.90,.05), lgd.name = "Augur\nAUC", direction = 1
)
ggsave(filename = "wt_augurscore_umap_update.eps", plot = u, device = "eps", path = "../figures/Fig3/", width = 4, height = 4, dpi = 300, family = "Arial")
