setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(tidydr)
library(RColorBrewer)

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})

if(!dir.exists("../figures/Fig3")) {dir.create("../figures/Fig3")}

endo <- SeuratDisk::LoadH5Seurat("../data/Wei_2024/endo.h5Seurat")

endo_wt <- subset(endo, subset = genotype == "WT") |>
  NormalizeData() |>
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA() |>
  FindNeighbors(dims = 1:20) |>
  FindClusters() |>
  RunUMAP(dims = 1:20)
Idents(endo_wt) <- "treatment"
p1 <- DimPlot(endo_wt, 
        group.by = "treatment",
        cols = brewer.pal(12,"Paired")[c(8,10)],
        pt.size = 1.5) +
  ggtitle(NULL) +
  NoAxes() +
  theme_dr() +
  theme(legend.position = "none",
        axis.title = element_text(size = 28, hjust = .05))
ggsave(filename = "endo_wt_umap.eps", plot = p1, device = "eps", path = "../figures/Fig3/", width = 6, height = 5, dpi = 300, family = "Arial")

# Violin Plots

## Casp4
gene = "Casp4"
p2 <- VlnPlot(endo_wt,
              features = gene,
              cols = brewer.pal(12,"Paired")[c(8,10)],
              pt.size = 1.5) +
  ylim(c(-0.01,3.2)) +
  scale_x_discrete(labels = c("","")) +
  geom_boxplot(width = .1, alpha = .3, color = "white") +
  stat_compare_means(method = "wilcox", method.args = list(alternative = "two.sided"), 
                     comparisons = list(c("LPS","PBS")), label = "p.signif", label.y = 2.5, bracket.size = 1, size = 10) +
  ylab(NULL) +
  ggtitle("Casp11") +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 28, hjust = .05, face = "bold.italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28),
        axis.text.y = element_text(size = 25))

## Cd14
gene = "Cd14"
p3 <- VlnPlot(endo_wt,
              features = gene,
              cols = brewer.pal(12,"Paired")[c(8,10)],
              pt.size = 1.5) +
  ylim(c(-0.01,3.2)) +
  geom_boxplot(width = .1, alpha = .3, color = "white") +
  stat_compare_means(method = "wilcox", method.args = list(alternative = "two.sided"), 
                     comparisons = list(c("LPS","PBS")), label = "p.signif", label.y = 2.5, bracket.size = 1, size = 10) +
  scale_x_discrete(labels = c("Ctrl", "LPS")) +
  ylab(NULL) +
  ggtitle(gene) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size = 28, hjust = .05, face = "bold.italic"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 28, color = "black"),
        axis.text.y = element_text(size = 25))

# FeaturePlots
## Casp4
gene = "Casp4"
f2 <- FeaturePlot(endo_wt,
                  features = gene,
                  pt.size = 1.2) +
  NoAxes() +
  ggtitle(NULL) +
  scale_color_distiller(palette = "RdPu", direction = 1, breaks = c(0,1,2), labels = scales::label_number(accuracy = 0.1)) +
  theme(legend.position = c(.05, .27),
        legend.text = element_text(size = 20),
        panel.border = element_rect(size = .6, color = "black"))

## Cd14
gene = "Cd14"
f3 <- FeaturePlot(endo_wt,
                  features = gene,
                  pt.size = 1.2) +
  NoAxes() +
  ggtitle(NULL) +
  scale_color_distiller(palette = "RdPu", direction = 1, breaks = c(0,1,2), labels = scales::label_number(accuracy = 0.1)) +
  theme(legend.position = c(.05, .27),
        legend.text = element_text(size = 20),
        panel.border = element_rect(size = .6, color = "black"))

pdf(file = "../figures/Fig3/endo_wt_Casp4.pdf", width = 8, height = 4, family = "ArialMT")
p2 | f2
dev.off()

pdf(file = "../figures/Fig3/endo_wt_Cd14.pdf", width = 8, height = 4, family = "ArialMT")
p3 | f3
dev.off()