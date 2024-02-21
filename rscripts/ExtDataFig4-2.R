setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(tidyverse)
library(Seurat)
library(RColorBrewer)

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})

suppressMessages({
  ex_neu <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/ex_neu.h5Seurat", verbose = F)
  inh_neu <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/inh_neu.h5Seurat", verbose = F)
  astro <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/astro.h5Seurat", verbose = F)
  oligo <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/oligo.h5Seurat", verbose = F)
  micro <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/micro.h5Seurat", verbose = F)
  endo <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/endo.h5Seurat", verbose = F)
})

degs.cell_type.list <- readRDS("../outs/de_analysis/degs_cell_type_list.rds")
degs.cell_subtype.list <- readRDS("../outs/de_analysis/degs_cell_subtype_list.rds")

# =======================================================================================================================
## Microglia
d <- DimPlot(micro, group.by = "sample", pt.size = 1.5) +
  NoAxes() +
  ggtitle(NULL) +
  scale_color_manual(
    values = brewer.pal(8,"Set1")[c(2,1,3,4)],
    labels = c("WT PBS", "WT LPS", bquote(paste(italic('Gsdmd'^{"-/-"}),' PBS')), bquote(paste(italic('Gsdmd'^{"-/-"}),' LPS')))
    ) +
  guides(color=guide_legend(nrow=2, byrow=FALSE)) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    legend.justification = "center"
  )

ABS_LFC_cutoff = 1
P_ADJ_cutoff = 0.01

v_wt <- degs.cell_subtype.list[["Micro"]][["WT"]] |>
  ggplot(mapping=aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(size = .8) +
  geom_hline(yintercept = -log10(P_ADJ_cutoff), linetype="dashed") +
  geom_vline(xintercept = ABS_LFC_cutoff, linetype="dashed") +
  geom_vline(xintercept = -ABS_LFC_cutoff, linetype="dashed") +
  xlab(NULL) + ylab(NULL) +
  # xlab(label = expression(paste('Average Log'[2]*'Fold-Change (LPS vs. Ctrl)'))) +
  # ylab(label = expression(paste('-Log'[10]*'(p-adjust)'))) +
  ggtitle('WT') +
  annotate("text", max(degs.cell_subtype.list[["Micro"]][["WT"]]$avg_log2FC)-3, Inf, 
           label = paste("up:",degs.cell_subtype.list[["Micro"]][["WT"]] |> dplyr::filter(avg_log2FC >= ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow(),
                         "\ndown:",degs.cell_subtype.list[["Micro"]][["WT"]] |> dplyr::filter(avg_log2FC <= -ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow()), 
           hjust = 0, vjust = 1, size = 5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, family = "Arial"), axis.text.x = element_text(size = 13, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"), axis.text.y = element_text(size = 13, family = "Arial"),
    plot.title = element_text(size = 18), 
    plot.margin = margin(t=0,r=10,b=0,l=0)
  )

v_gdko <- degs.cell_subtype.list[["Micro"]][["GDKO"]] |>
  ggplot(mapping=aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(size = .8) +
  geom_hline(yintercept = -log10(P_ADJ_cutoff), linetype="dashed") +
  geom_vline(xintercept = ABS_LFC_cutoff, linetype="dashed") +
  geom_vline(xintercept = -ABS_LFC_cutoff, linetype="dashed") +
  xlab(NULL) + ylab(NULL) +
  # xlab(label = expression(paste('Average Log'[2]*'Fold-Change (LPS vs. Ctrl)'))) +
  # ylab(label = expression(paste('-Log'[10]*'(p-adjust)'))) +
  ggtitle(bquote(paste(italic('Gsdmd'^{"-/-"})))) +
  annotate("text", max(degs.cell_subtype.list[["Micro"]][["GDKO"]]$avg_log2FC)-3, Inf, 
           label = paste("up:",degs.cell_subtype.list[["Micro"]][["GDKO"]] |> dplyr::filter(avg_log2FC >= ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow(),
                         "\ndown:",degs.cell_subtype.list[["Micro"]][["GDKO"]] |> dplyr::filter(avg_log2FC <= -ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow()), 
           hjust = 0, vjust = 1, size = 5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, family = "Arial"), axis.text.x = element_text(size = 13, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"), axis.text.y = element_text(size = 13, family = "Arial"),
    plot.title = element_text(size = 18), 
    plot.margin = margin(t=0,r=10,b=0,l=0)
  )

d / v_wt / v_gdko
ggsave(
  filename = "ExtDataFig4_micro_umapvolcano.eps",
  device = "eps", plot = last_plot(), path = "../figures/ExtDataFig4/", 
  width = 4, height = 9, dpi = 300, 
  family = "Arial"
)

# =======================================================================================================================
## Astrocytes
d <- DimPlot(astro, group.by = "sample", pt.size = 1.5) +
  NoAxes() +
  ggtitle(NULL) +
  scale_color_manual(
    values = brewer.pal(8,"Set1")[c(2,1,3,4)],
    labels = c("WT PBS", "WT LPS", bquote(paste(italic('Gsdmd'^{"-/-"}),' PBS')), bquote(paste(italic('Gsdmd'^{"-/-"}),' LPS')))
  ) +
  guides(color=guide_legend(nrow=2, byrow=FALSE)) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    legend.justification = "center"
  )

ABS_LFC_cutoff = 1
P_ADJ_cutoff = 0.01

v_wt <- degs.cell_subtype.list[["Astro"]][["WT"]] |>
  ggplot(mapping=aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(size = .8) +
  geom_hline(yintercept = -log10(P_ADJ_cutoff), linetype="dashed") +
  geom_vline(xintercept = ABS_LFC_cutoff, linetype="dashed") +
  geom_vline(xintercept = -ABS_LFC_cutoff, linetype="dashed") +
  xlab(NULL) + ylab(NULL) +
  # xlab(label = expression(paste('Average Log'[2]*'Fold-Change (LPS vs. Ctrl)'))) +
  # ylab(label = expression(paste('-Log'[10]*'(p-adjust)'))) +
  ggtitle('WT') +
  annotate("text", max(degs.cell_subtype.list[["Astro"]][["WT"]]$avg_log2FC)-2, Inf, 
           label = paste("up:",degs.cell_subtype.list[["Astro"]][["WT"]] |> dplyr::filter(avg_log2FC >= ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow(),
                         "\ndown:",degs.cell_subtype.list[["Astro"]][["WT"]] |> dplyr::filter(avg_log2FC <= -ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow()), 
           hjust = 0, vjust = 1, size = 5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, family = "Arial"), axis.text.x = element_text(size = 13, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"), axis.text.y = element_text(size = 13, family = "Arial"),
    plot.title = element_text(size = 18), 
    plot.margin = margin(t=0,r=10,b=0,l=0)
  )

v_gdko <- degs.cell_subtype.list[["Astro"]][["GDKO"]] |>
  ggplot(mapping=aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(size = .8) +
  geom_hline(yintercept = -log10(P_ADJ_cutoff), linetype="dashed") +
  geom_vline(xintercept = ABS_LFC_cutoff, linetype="dashed") +
  geom_vline(xintercept = -ABS_LFC_cutoff, linetype="dashed") +
  xlab(NULL) + ylab(NULL) +
  # xlab(label = expression(paste('Average Log'[2]*'Fold-Change (LPS vs. Ctrl)'))) +
  # ylab(label = expression(paste('-Log'[10]*'(p-adjust)'))) +
  ggtitle(bquote(paste(italic('Gsdmd'^{"-/-"})))) +
  annotate("text", max(degs.cell_subtype.list[["Astro"]][["GDKO"]]$avg_log2FC)-2.8, Inf, 
           label = paste("up:",degs.cell_subtype.list[["Astro"]][["GDKO"]] |> dplyr::filter(avg_log2FC >= ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow(),
                         "\ndown:",degs.cell_subtype.list[["Astro"]][["GDKO"]] |> dplyr::filter(avg_log2FC <= -ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow()), 
           hjust = 0, vjust = 1, size = 5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, family = "Arial"), axis.text.x = element_text(size = 13, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"), axis.text.y = element_text(size = 13, family = "Arial"),
    plot.title = element_text(size = 18), 
    plot.margin = margin(t=0,r=10,b=0,l=0)
  )

d / v_wt / v_gdko
ggsave(
  filename = "ExtDataFig4_astro_umapvolcano.eps",
  device = "eps", plot = last_plot(), path = "../figures/ExtDataFig4/", 
  width = 4, height = 9, dpi = 300, 
  family = "Arial"
)

# =======================================================================================================================
## Oligodendrocytes
d1 <- DimPlot(oligo, group.by = "cell_subtype", pt.size = 1.5) +
  NoAxes() +
  ggtitle(NULL) +
  scale_color_manual(
    values = scales::hue_pal()(19)[c(4,5,6,7)]
  ) +
  guides(color=guide_legend(nrow=2, byrow=FALSE)) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    legend.justification = "center"
  )

d2 <- DimPlot(oligo, group.by = "sample", pt.size = 1.5) +
  NoAxes() +
  ggtitle(NULL) +
  scale_color_manual(
    values = brewer.pal(8,"Set1")[c(2,1,3,4)],
    labels = c("WT PBS", "WT LPS", bquote(paste(italic('Gsdmd'^{"-/-"}),' PBS')), bquote(paste(italic('Gsdmd'^{"-/-"}),' LPS')))
  ) +
  guides(color=guide_legend(nrow=2, byrow=FALSE)) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 15),
    legend.justification = "center"
  )

ABS_LFC_cutoff = 1
P_ADJ_cutoff = 0.01

v_wt <- degs.cell_subtype.list[["MOLs"]][["WT"]] |>
  ggplot(mapping=aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(size = .8) +
  geom_hline(yintercept = -log10(P_ADJ_cutoff), linetype="dashed") +
  geom_vline(xintercept = ABS_LFC_cutoff, linetype="dashed") +
  geom_vline(xintercept = -ABS_LFC_cutoff, linetype="dashed") +
  xlab(NULL) + ylab(NULL) +
  # xlab(label = expression(paste('Average Log'[2]*'Fold-Change (LPS vs. Ctrl)'))) +
  # ylab(label = expression(paste('-Log'[10]*'(p-adjust)'))) +
  ggtitle('WT MOLs') +
  annotate("text", max(degs.cell_subtype.list[["MOLs"]][["WT"]]$avg_log2FC)-2, Inf, 
           label = paste("up:",degs.cell_subtype.list[["MOLs"]][["WT"]] |> dplyr::filter(avg_log2FC >= ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow(),
                         "\ndown:",degs.cell_subtype.list[["MOLs"]][["WT"]] |> dplyr::filter(avg_log2FC <= -ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow()), 
           hjust = 0, vjust = 1, size = 5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, family = "Arial"), axis.text.x = element_text(size = 13, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"), axis.text.y = element_text(size = 13, family = "Arial"),
    plot.title = element_text(size = 18), 
    plot.margin = margin(t=0,r=10,b=0,l=0)
  )

v_gdko <- degs.cell_subtype.list[["MOLs"]][["GDKO"]] |>
  ggplot(mapping=aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(size = .8) +
  geom_hline(yintercept = -log10(P_ADJ_cutoff), linetype="dashed") +
  geom_vline(xintercept = ABS_LFC_cutoff, linetype="dashed") +
  geom_vline(xintercept = -ABS_LFC_cutoff, linetype="dashed") +
  xlab(NULL) + ylab(NULL) +
  # xlab(label = expression(paste('Average Log'[2]*'Fold-Change (LPS vs. Ctrl)'))) +
  # ylab(label = expression(paste('-Log'[10]*'(p-adjust)'))) +
  ggtitle(bquote(paste(italic('Gsdmd'^{"-/-"}),' MOLs'))) +
  annotate("text", max(degs.cell_subtype.list[["MOLs"]][["GDKO"]]$avg_log2FC)-2.2, Inf, 
           label = paste("up:",degs.cell_subtype.list[["MOLs"]][["GDKO"]] |> dplyr::filter(avg_log2FC >= ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow(),
                         "\ndown:",degs.cell_subtype.list[["MOLs"]][["GDKO"]] |> dplyr::filter(avg_log2FC <= -ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow()), 
           hjust = 0, vjust = 1, size = 5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, family = "Arial"), axis.text.x = element_text(size = 13, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"), axis.text.y = element_text(size = 13, family = "Arial"),
    plot.title = element_text(size = 18), 
    plot.margin = margin(t=0,r=10,b=0,l=0)
  )

(d1 / d2)
ggsave(
  filename = "umap_mols.eps",
  device = "eps", plot = last_plot(), path = "../figures/ExtDataFig4/", 
  width = 4, height = 7, dpi = 300, 
  family = "Arial"
)
(v_wt / v_gdko)
ggsave(
  filename = "volcano_mols.eps",
  device = "eps", plot = last_plot(), path = "../figures/ExtDataFig4/", 
  width = 4, height = 7, dpi = 300, 
  family = "Arial"
)

# =======================================================================================================================
## Neurons

ABS_LFC_cutoff = 1
P_ADJ_cutoff = 0.01

v_wt <- degs.cell_type.list[["Ex-Neu."]][["WT"]] |>
  ggplot(mapping=aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(size = .8) +
  geom_hline(yintercept = -log10(P_ADJ_cutoff), linetype="dashed") +
  geom_vline(xintercept = ABS_LFC_cutoff, linetype="dashed") +
  geom_vline(xintercept = -ABS_LFC_cutoff, linetype="dashed") +
  xlab(NULL) + ylab(NULL) +
  # xlab(label = expression(paste('Average Log'[2]*'Fold-Change (LPS vs. Ctrl)'))) +
  # ylab(label = expression(paste('-Log'[10]*'(p-adjust)'))) +
  ggtitle('WT') +
  annotate("text", max(degs.cell_type.list[["Ex-Neu."]][["WT"]]$avg_log2FC)-.6, Inf, 
           label = paste("up:",degs.cell_type.list[["Ex-Neu."]][["WT"]] |> dplyr::filter(avg_log2FC >= ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow(),
                         "\ndown:",degs.cell_type.list[["Ex-Neu."]][["WT"]] |> dplyr::filter(avg_log2FC <= -ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow()), 
           hjust = 0, vjust = 1, size = 5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, family = "Arial"), axis.text.x = element_text(size = 13, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"), axis.text.y = element_text(size = 13, family = "Arial"),
    plot.title = element_text(size = 18), 
    plot.margin = margin(t=0,r=10,b=0,l=0)
  )

v_gdko <- degs.cell_type.list[["Ex-Neu."]][["GDKO"]] |>
  ggplot(mapping=aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(size = .8) +
  geom_hline(yintercept = -log10(P_ADJ_cutoff), linetype="dashed") +
  geom_vline(xintercept = ABS_LFC_cutoff, linetype="dashed") +
  geom_vline(xintercept = -ABS_LFC_cutoff, linetype="dashed") +
  xlab(NULL) + ylab(NULL) +
  # xlab(label = expression(paste('Average Log'[2]*'Fold-Change (LPS vs. Ctrl)'))) +
  # ylab(label = expression(paste('-Log'[10]*'(p-adjust)'))) +
  ggtitle(bquote(paste(italic('Gsdmd'^{"-/-"})))) +
  annotate("text", max(degs.cell_type.list[["Ex-Neu."]][["GDKO"]]$avg_log2FC)-.6, Inf, 
           label = paste("up:",degs.cell_type.list[["Ex-Neu."]][["GDKO"]] |> dplyr::filter(avg_log2FC >= ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow(),
                         "\ndown:",degs.cell_type.list[["Ex-Neu."]][["GDKO"]] |> dplyr::filter(avg_log2FC <= -ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow()), 
           hjust = 0, vjust = 1, size = 5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, family = "Arial"), axis.text.x = element_text(size = 13, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"), axis.text.y = element_text(size = 13, family = "Arial"),
    plot.title = element_text(size = 18), 
    plot.margin = margin(t=0,r=10,b=0,l=0)
  )

v_wt / v_gdko
ggsave(
  filename = "ExtDataFig4_ex_neu_volcano.eps",
  device = "eps", plot = last_plot(), path = "../figures/ExtDataFig4/", 
  width = 4, height = 7, dpi = 300, 
  family = "Arial"
)

ABS_LFC_cutoff = 1
P_ADJ_cutoff = 0.01

v_wt <- degs.cell_type.list[["Inh-Neu."]][["WT"]] |>
  ggplot(mapping=aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(size = .8) +
  geom_hline(yintercept = -log10(P_ADJ_cutoff), linetype="dashed") +
  geom_vline(xintercept = ABS_LFC_cutoff, linetype="dashed") +
  geom_vline(xintercept = -ABS_LFC_cutoff, linetype="dashed") +
  xlab(NULL) + ylab(NULL) +
  # xlab(label = expression(paste('Average Log'[2]*'Fold-Change (LPS vs. Ctrl)'))) +
  # ylab(label = expression(paste('-Log'[10]*'(p-adjust)'))) +
  ggtitle('WT') +
  annotate("text", max(degs.cell_type.list[["Inh-Neu."]][["WT"]]$avg_log2FC)-.6, Inf, 
           label = paste("up:",degs.cell_type.list[["Inh-Neu."]][["WT"]] |> dplyr::filter(avg_log2FC >= ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow(),
                         "\ndown:",degs.cell_type.list[["Inh-Neu."]][["WT"]] |> dplyr::filter(avg_log2FC <= -ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow()), 
           hjust = 0, vjust = 1, size = 5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, family = "Arial"), axis.text.x = element_text(size = 13, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"), axis.text.y = element_text(size = 13, family = "Arial"),
    plot.title = element_text(size = 18), 
    plot.margin = margin(t=0,r=10,b=0,l=0)
  )

v_gdko <- degs.cell_type.list[["Inh-Neu."]][["GDKO"]] |>
  ggplot(mapping=aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  geom_point(size = .8) +
  geom_hline(yintercept = -log10(P_ADJ_cutoff), linetype="dashed") +
  geom_vline(xintercept = ABS_LFC_cutoff, linetype="dashed") +
  geom_vline(xintercept = -ABS_LFC_cutoff, linetype="dashed") +
  xlab(NULL) + ylab(NULL) +
  # xlab(label = expression(paste('Average Log'[2]*'Fold-Change (LPS vs. Ctrl)'))) +
  # ylab(label = expression(paste('-Log'[10]*'(p-adjust)'))) +
  ggtitle(bquote(paste(italic('Gsdmd'^{"-/-"})))) +
  annotate("text", max(degs.cell_type.list[["Inh-Neu."]][["GDKO"]]$avg_log2FC)-.6, Inf, 
           label = paste("up:",degs.cell_type.list[["Inh-Neu."]][["GDKO"]] |> dplyr::filter(avg_log2FC >= ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow(),
                         "\ndown:",degs.cell_type.list[["Inh-Neu."]][["GDKO"]] |> dplyr::filter(avg_log2FC <= -ABS_LFC_cutoff & p_val_adj <= P_ADJ_cutoff) |> nrow()), 
           hjust = 0, vjust = 1, size = 5) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 15, family = "Arial"), axis.text.x = element_text(size = 13, family = "Arial"),
    axis.title.y = element_text(size = 15, family = "Arial"), axis.text.y = element_text(size = 13, family = "Arial"),
    plot.title = element_text(size = 18), 
    plot.margin = margin(t=0,r=10,b=0,l=0)
  )

v_wt / v_gdko
ggsave(
  filename = "ExtDataFig4_inh_neu_volcano.eps",
  device = "eps", plot = last_plot(), path = "../figures/ExtDataFig4/", 
  width = 4, height = 7, dpi = 300, 
  family = "Arial"
)
