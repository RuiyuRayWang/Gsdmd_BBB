setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(tidyverse)
library(Seurat)
library(RColorBrewer)
library(ggpubr)
source("avgLogExp.R")

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})

suppressMessages({
  # ex_neu <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/ex_neu.h5Seurat", verbose = F)
  # inh_neu <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/inh_neu.h5Seurat", verbose = F)
  astro <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/astro.h5Seurat", verbose = F)
  oligo <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/oligo.h5Seurat", verbose = F)
  # opcs <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/opcs.h5Seurat", verbose = F)
  micro <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/micro.h5Seurat", verbose = F)
  endo <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/endo.h5Seurat", verbose = F)
})

degs.cell_subtype.lps.wt_gdko.list <- readRDS("../outs/de_analysis/degs_cell_subtype_lps_wt_gdko_list.rds")
csvd.lps.mks.cell_subtype <- readRDS("../outs/de_analysis/csvd_lps_mks_cell_subtype.rds")

## Endothelial Cells
x = degs.cell_subtype.lps.wt_gdko.list[["Endo"]] |>
  slice_max(order_by = avg_log2FC, n = 200) |>
  pull("gene")
y = degs.cell_subtype.lps.wt_gdko.list[["Endo"]] |>
  slice_min(order_by = avg_log2FC, n = 200) |>
  pull("gene")

z = csvd.lps.mks.cell_subtype[["Endo"]] |>
  dplyr::filter(abs(delta_log2FC) < .1) |>
  dplyr::filter(sum_log2FC < 0) |>
  slice_min(order_by = sum_log2FC, n = 200) |>
  pull("gene")
w = csvd.lps.mks.cell_subtype[["Endo"]] |>
  dplyr::filter(abs(delta_log2FC) < .1) |>
  dplyr::filter(sum_log2FC > 0) |>
  slice_max(order_by = sum_log2FC, n = 200) |>
  pull("gene")

length(intersect(c(x,y), c(z,w))) == 0

endo$sample = factor(endo$sample, levels = c("WT PBS", "WT LPS", "GDKO PBS", "GDKO LPS"))
endo <- ScaleData(endo, features = c(x,y,z,w))
endo_h1 = DoHeatmap(
  object = endo, 
  features = c(x,y),
  group.by = "sample", 
  group.colors = brewer.pal(8,"Set1")[c(2,1,3,4)],
  group.bar = T, label = F
) +
  # scale_fill_distiller(palette = "PiYG", direction = -1) +
  # scale_fill_viridis_c(option = "cividis", direction = 1) +
  guides(color = "none") +
  theme(
    axis.text.y = element_blank(), 
    plot.margin = margin(t=0,r=0,b=0,l=0)
  )
endo_h2 = DoHeatmap(
  object = endo, 
  features = c(z,w),
  group.by = "sample",
  group.colors = brewer.pal(8,"Set1")[c(2,1,3,4)],
  group.bar = F, label = F
) +
  # scale_fill_distiller(palette = "PiYG", direction = -1) +
  # scale_fill_viridis_c(option = "cividis", direction = 1) +
  guides(color = "none") +
  theme(
    axis.text.y = element_blank(), 
    plot.margin = margin(t=0,r=0,b=0,l=0)
  )
# pdf(
#   file = "../figures/FigSZ/FigSZ_endo.pdf", 
# )
# ggarrange(endo_h1, endo_h2, common.legend = TRUE, ncol = 1, legend = "bottom")
# dev.off()
ggarrange(endo_h1, endo_h2, common.legend = TRUE, ncol = 1, legend = "bottom") |>
  ggexport(filename = "../figures/ExtDataFig5/ExtDataFig_endo.pdf", family = "Arial")


## Microglia
x = degs.cell_subtype.lps.wt_gdko.list[["Micro"]] |>
  slice_max(order_by = avg_log2FC, n = 200) |>
  pull("gene")
y = degs.cell_subtype.lps.wt_gdko.list[["Micro"]] |>
  slice_min(order_by = avg_log2FC, n = 200) |>
  pull("gene")

z = csvd.lps.mks.cell_subtype[["Micro"]] |>
  dplyr::filter(abs(delta_log2FC) < .1) |>
  dplyr::filter(sum_log2FC < 0) |>
  slice_min(order_by = sum_log2FC, n = 200) |>
  pull("gene")
w = csvd.lps.mks.cell_subtype[["Micro"]] |>
  dplyr::filter(abs(delta_log2FC) < .1) |>
  dplyr::filter(sum_log2FC > 0) |>
  slice_max(order_by = sum_log2FC, n = 200) |>
  pull("gene")

length(intersect(c(x,y), c(z,w))) == 0

micro$sample = factor(micro$sample, levels = c("WT PBS", "WT LPS", "GDKO PBS", "GDKO LPS"))
micro <- ScaleData(micro, features = c(x,y,z,w))
micro_h1 = DoHeatmap(
  object = micro, 
  features = c(x,y),
  group.by = "sample", 
  group.colors = brewer.pal(8,"Set1")[c(2,1,3,4)],
  group.bar = T, label = F
) +
  # scale_fill_viridis_c(option = "cividis", direction = 1) +
  guides(color = "none") +
  theme(
    axis.text.y = element_blank(), 
    plot.margin = margin(t=0,r=0,b=0,l=0)
  )
micro_h2 = DoHeatmap(
  object = micro, 
  features = c(z,w),
  group.by = "sample",
  group.colors = brewer.pal(8,"Set1")[c(2,1,3,4)],
  group.bar = F, label = F
) +
  # scale_fill_viridis_c(option = "cividis", direction = 1) +
  guides(color = "none") +
  theme(
    axis.text.y = element_blank(), 
    plot.margin = margin(t=0,r=0,b=0,l=0)
  )
# pdf(
#   file = "../figures/FigSZ/FigSZ_micro.pdf", 
# )
# ggpubr::ggarrange(micro_h1, micro_h2, common.legend = TRUE, ncol = 1, legend = "bottom")
# dev.off()
ggarrange(micro_h1, micro_h2, common.legend = TRUE, ncol = 1, legend = "bottom") |>
  ggexport(filename = "../figures/ExtDataFig5/ExtDataFig_micro.pdf", family = "Arial")


## Astrocytes
x = degs.cell_subtype.lps.wt_gdko.list[["Astro"]] |>
  slice_max(order_by = avg_log2FC, n = 200) |>
  pull("gene")
y = degs.cell_subtype.lps.wt_gdko.list[["Astro"]] |>
  slice_min(order_by = avg_log2FC, n = 200) |>
  pull("gene")

z = csvd.lps.mks.cell_subtype[["Astro"]] |>
  dplyr::filter(abs(delta_log2FC) < .1) |>
  dplyr::filter(sum_log2FC < 0) |>
  slice_min(order_by = sum_log2FC, n = 200) |>
  pull("gene")
w = csvd.lps.mks.cell_subtype[["Astro"]] |>
  dplyr::filter(abs(delta_log2FC) < .1) |>
  dplyr::filter(sum_log2FC > 0) |>
  slice_max(order_by = sum_log2FC, n = 200) |>
  pull("gene")

astro$sample = factor(astro$sample, levels = c("WT PBS", "WT LPS", "GDKO PBS", "GDKO LPS"))
astro <- ScaleData(astro, features = c(x,y,z,w))
astro_h1 = DoHeatmap(
  object = astro, 
  features = c(x,y),
  group.by = "sample",
  group.colors = brewer.pal(8,"Set1")[c(2,1,3,4)],
  group.bar = T, label = F
) +
  # scale_fill_viridis_c(option = "cividis", direction = 1) +
  guides(color = "none") +
  theme(
    axis.text.y = element_blank(), 
    plot.margin = margin(t=0,r=0,b=0,l=0)
  )
astro_h2 = DoHeatmap(
  object = astro, 
  features = c(z,w),
  group.by = "sample",
  group.colors = brewer.pal(8,"Set1")[c(2,1,3,4)],
  group.bar = F, label = F
) +
  # scale_fill_viridis_c(option = "cividis", direction = 1) +
  guides(color = "none") +
  theme(
    axis.text.y = element_blank(), 
    plot.margin = margin(t=0,r=0,b=0,l=0)
  )

# pdf(
#   file = "../figures/FigSZ/FigSZ_astro.pdf", 
# )
# ggpubr::ggarrange(astro_h1, astro_h2, common.legend = TRUE, ncol = 1, legend = "bottom")
# dev.off()
ggarrange(astro_h1, astro_h2, common.legend = TRUE, ncol = 1, legend = "bottom") |>
  ggexport(filename = "../figures/ExtDataFig5/ExtDataFig_astro.pdf", family = "Arial")


## MOLs
mols <- subset(oligo, subset = cell_subtype %in% c("MOL1","MOL2"))

### MOL1
x = degs.cell_subtype.lps.wt_gdko.list[["MOLs"]] |>
  slice_max(order_by = avg_log2FC, n = 200) |>
  pull("gene")
y = degs.cell_subtype.lps.wt_gdko.list[["MOLs"]] |>
  slice_min(order_by = avg_log2FC, n = 200) |>
  pull("gene")

z = csvd.lps.mks.cell_subtype[["MOLs"]] |>
  dplyr::filter(abs(delta_log2FC) < .1) |>
  dplyr::filter(sum_log2FC < 0) |>
  slice_min(order_by = sum_log2FC, n = 200) |>
  pull("gene")
w = csvd.lps.mks.cell_subtype[["MOLs"]] |>
  dplyr::filter(abs(delta_log2FC) < .1) |>
  dplyr::filter(sum_log2FC > 0) |>
  slice_max(order_by = sum_log2FC, n = 200) |>
  pull("gene")

length(intersect(c(x,y), c(z,w))) == 0

mols$sample = factor(mols$sample, levels = c("WT PBS", "WT LPS", "GDKO PBS", "GDKO LPS"))
mols <- ScaleData(mols, features = c(x,y,z,w))
mols_h1 = DoHeatmap(
  object = mols, 
  features = c(x,y),
  group.by = "sample", 
  group.colors = brewer.pal(8,"Set1")[c(2,1,3,4)],
  group.bar = T, label = F
) +
  # scale_fill_viridis_c(option = "cividis", direction = 1) +
  guides(color = "none") +
  theme(
    axis.text.y = element_blank(), 
    plot.margin = margin(t=0,r=0,b=0,l=0)
  )
mols_h2 = DoHeatmap(
  object = mols, 
  features = c(z,w),
  group.by = "sample",
  group.colors = brewer.pal(8,"Set1")[c(2,1,3,4)],
  group.bar = F, label = F
) +
  # scale_fill_viridis_c(option = "cividis", direction = 1) +
  guides(color = "none") +
  theme(
    axis.text.y = element_blank(), 
    plot.margin = margin(t=0,r=0,b=0,l=0)
  )
ggarrange(mols_h1, mols_h2, common.legend = TRUE, ncol = 1, legend = "bottom") |>
  ggexport(filename = "../figures/ExtDataFig5/ExtDataFig_mols.pdf", family = "Arial")

