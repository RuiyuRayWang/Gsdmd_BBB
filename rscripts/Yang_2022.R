setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(Seurat)
library(tidyverse)
library(tidydr)

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})

mat <- Read10X("../data/Yang_2022/")
meta <- read_tsv("../data/Yang_2022/meta.tsv") |> column_to_rownames("Cell")
umap_coords <- read_tsv("../data/Yang_2022/UMAP.coords.tsv.gz", col_names = c("Cell","UMAP_1","UMAP_2")) |> column_to_rownames("Cell")

cells <- CreateSeuratObject(mat, meta.data = meta)

endo <- subset(cells, subset = Cell_Type %in% c("Arterial","Veinous","Capillary"))
endo <- endo |>
  NormalizeData() |> 
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA() |>
  FindNeighbors(dims = 1:30) |>
  RunUMAP(dims = 1:30)
d1 <- DimPlot(endo, label = T, label.box = T, repel = T, label.size = 5)
spirious_labels <- CellSelector(d1)

endo <- subset(cells, cells = Cells(endo)[!Cells(endo) %in% spirious_labels])
endo <- endo |>
  NormalizeData() |> 
  FindVariableFeatures() |>
  ScaleData() |>
  RunPCA() |>
  FindNeighbors(dims = 1:30) |>
  RunUMAP(dims = 1:30)

Idents(endo) <- "Cell_Type"
d1 <- DimPlot(endo, cols = scales::hue_pal()(6)[c(5,3,1)], label = T, label.box = T, repel = T, label.size = 5) +
  ggtitle("Yang et al. (2022):\nhuman brain vasculature atlas") + 
  theme_dr(xlength = .2, ylength = .2) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 15),
    axis.title = element_text(size = 15, hjust = 0)
  )
ggsave(
  filename = "yang_vas_1.eps", plot = d1, device = "eps", path = "../figures/FigS16/", width = 5, height = 5, dpi = 300, family = "Arial"
)

f1 <- FeaturePlot(endo, features = "GSDMD", raster = F) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.text = element_text(size = 13)
  )
ggsave(
  filename = "yang_vas_2.eps", plot = f1, device = "eps", path = "../figures/FigS16/", width = 5.6, height = 5, dpi = 300, family = "Arial"
)

f2 <- FeaturePlot(endo, features = "CASP4", raster = F) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 15),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.text = element_text(size = 13)
  )
ggsave(
  filename = "yang_vas_3.eps", plot = f2, device = "eps", path = "../figures/FigS16/", width = 5.6, height = 5, dpi = 300, family = "Arial"
)

d1 | f1 | f2

vas <- subset(cells, subset = Cell_Type %in% c("Arterial","Veinous","Capillary","SMC","Pericyte","P. Fibro","M. Fibro","Ependymal"))
vas$Cell_Type <- factor(vas$Cell_Type, levels = c("Arterial","Veinous","Capillary","SMC","Pericyte","P. Fibro","M. Fibro","Ependymal"))

vas <- vas |> NormalizeData()

VlnPlot(vas, features = "GSDMD", group.by = "Cell_Type") +
  NoLegend() +
  ggtitle("GSDMD", subtitle = "Yang et al. (2022)") +
  ylab("Normalized Expression") +
  theme(
    plot.title = element_text(hjust = 0),
    axis.title.x = element_blank()
  )
ggsave("Yang_GSDMD_expr.eps", path = "../figures/FigS16/", width = 6, height = 4, dpi = 300)
