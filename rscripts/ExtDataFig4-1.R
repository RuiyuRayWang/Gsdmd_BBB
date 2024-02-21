setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Seurat)
library(tidyverse)
library(tidydr)
library(ggpubr)

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})

if(!dir.exists("../figures/ExtDataFig4")) {dir.create("../figures/ExtDataFig4")}

cells.clean <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/cells_clean.h5Seurat")
set.seed(42)

# nGenes violin
# ExtDataFig.4c
p1 <- cells.clean |> 
  FetchData(vars = c("nFeature_RNA","sample")) |>
  ggplot(mapping = aes(x=sample, y=nFeature_RNA)) +
  geom_violin(mapping = aes(fill = sample), color = "black") +
  geom_jitter(size = 0.3, height = 0, alpha = .3) +
  scale_x_discrete(labels = c(bquote(atop("WT","PBS")),bquote(atop("WT","LPS")),bquote(atop(italic('Gsdmd'^{"-/-"}),"PBS")),bquote(atop(italic('Gsdmd'^{"-/-"}),"LPS")))) +
  ylab("Number of genes detected") +
  theme_linedraw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 18, hjust = .5),
    axis.title.y = element_text(size = 21),
    axis.text.y = element_text(size = 18)
  )
ggexport(p1, filename = "../figures/ExtDataFig4/ExtDataFig4C.pdf", family = "Arial", width = 6, height = 5)

# Cell type composition across samples
# ExtDataFig. 4e
Idents(cells.clean) <- "cell_type"
cells.clean <- RenameIdents(cells.clean, `Proliferating OPCs` = "Prolif. OPCs")
cells.clean$cell_type_brief <- Idents(cells.clean)
cells.clean$cell_type_brief <- factor(cells.clean$cell_type_brief, 
                                      levels = c("Ex-Neu.","Inh-Neu.","Endo","Peri","SMCs","VLMCs","Astro","OPCs","Prolif. OPCs","Oligo","Micro","PVMs"))
p2 = cells.clean |> 
  FetchData(vars = c("cell_type_brief","sample")) |>
  dplyr::group_by(sample, cell_type_brief) |>
  dplyr::summarise(n = n()) |>
  mutate(prop = n/sum(n)*100) |>
  ggplot(mapping = aes(fill=cell_type_brief, y=prop, x=sample)) +
  geom_bar(position="fill", stat = "identity", color = "white", size = .2, width = .7) +
  ylab("Proportion (%)") +
  scale_fill_manual(name = "Cell types", values = scales::hue_pal()(12)) +
  scale_x_discrete(labels = c(bquote(atop("WT","PBS")),bquote(atop("WT","LPS")),bquote(atop(italic('Gsdmd'^{"-/-"}),"PBS")),bquote(atop(italic('Gsdmd'^{"-/-"}),"LPS")))) +
  scale_y_continuous(labels = scales::percent) +
  theme_linedraw() +
  theme(
    legend.position = "left",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 19, angle = 90, hjust = .5, vjust = .5),
    axis.title.y = element_text(size = 24),
    axis.text.y = element_text(size = 18)
  )
ggexport(p2, filename = "../figures/ExtDataFig4/ExtDataFig4e.pdf", family = "Arial", width = 6, height = 6)

# Preprocessed UMAP
# ExtDataFig. 4f
p3 = DimPlot(
  cells.clean, 
  group.by = "cell_type_brief",
  label = T,
  reduction = "umap",
  repel = T,
  label.size = 8
) +
  NoAxes() +
  theme_dr(xlength = .2, ylength = .2) +
  ggtitle("Preprocessed dataset") +
  theme(
    legend.position = "none", 
    panel.grid = element_blank(),
    axis.title = element_text(size = 17, hjust = 0),
    plot.title = element_text(size = 24, face = "bold"),
    legend.text = element_text(size = 18)
  )

ggexport(p3, filename = "../figures/ExtDataFig4/ExtDataFig4f.pdf", family = "Arial")

ggarrange(p2, NULL, p3, widths = c(1.4,0.1,1.6), ncol = 3, common.legend = TRUE, legend = "right") |>
  ggexport(filename = "../figures/ExtDataFig4/cell_composition_umap.pdf", width = 13, height = 6, family = "Arial")