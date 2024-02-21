setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(tidyverse)

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})

ego <- readRDS("../outs/functional_analysis/ora/ego_bp_signature_genes_cell_subtype_dfs.rds")

# ============================================================================
# Micro
col_pal = RColorBrewer::brewer.pal(9,"Reds")
go_use = c("mononuclear cell differentiation","interferon-gamma production")
ego[["Micro"]][["x"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,5), expand = c(0,0)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "micro_ego_x.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig5/", width = 4, height = 1.6, dpi = 300, family = "Arial")


col_pal = RColorBrewer::brewer.pal(9,"Purples")
go_use = c("ATP metabolic process","phagocytosis")
ego[["Micro"]][["y"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,7), expand = c(0,0)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "micro_ego_y.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig5/", width = 4, height = 1.6, dpi = 300, family = "Arial")


col_pal = RColorBrewer::brewer.pal(9,"Blues")
go_use = c("maintenance of cell number")
ego[["Micro"]][["z"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,2), expand = c(0,0), breaks = c(0,1,2)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "micro_ego_z.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig5/", width = 4, height = 0.8, dpi = 300, family = "Arial")


col_pal = RColorBrewer::brewer.pal(9,"Greens")
go_use = c("positive regulation of cytokine production","positive regulation of innate immune response")
ego[["Micro"]][["w"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,5), expand = c(0,0), breaks = c(0,2,4,6)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "micro_ego_w.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig5/", width = 4, height = 1.6, dpi = 300, family = "Arial")

# ============================================================================
# Astro
col_pal = RColorBrewer::brewer.pal(9,"Reds")
go_use = c("synapse organization","glial cell differentiation")
ego[["Astro"]][["x"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,15), expand = c(0,0)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "astro_ego_x.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig14/ego/", width = 4, height = 1.6, dpi = 300, family = "Arial")


col_pal = RColorBrewer::brewer.pal(9,"Purples")
go_use = c("response to interferon-gamma","cytokine-mediated signaling pathway")
ego[["Astro"]][["y"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,12), expand = c(0,0)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "astro_ego_y.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig14/ego/", width = 4, height = 1.6, dpi = 300, family = "Arial")


col_pal = RColorBrewer::brewer.pal(9,"Blues")
go_use = c("Wnt signaling pathway","regulation of Notch signaling pathway")
ego[["Astro"]][["z"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,2), expand = c(0,0), breaks = c(0,1,2)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "astro_ego_z.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig14/ego/", width = 4, height = 1.6, dpi = 300, family = "Arial")


col_pal = RColorBrewer::brewer.pal(9,"Greens")
go_use = c("immune response-regulating signaling pathway","glycerolipid metabolic process")
ego[["Astro"]][["w"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,3), expand = c(0,0), breaks = c(0,2)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "astro_ego_w.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig14/ego/", width = 4, height = 1.6, dpi = 300, family = "Arial")

# ============================================================================
# MOLs
col_pal = RColorBrewer::brewer.pal(9,"Reds")
go_use = c("cell junction assembly","regulation of membrane potential")
ego[["MOLs"]][["x"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,10), expand = c(0,0), breaks = c(0,2,4,6,8,10)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "mols_ego_x.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig5/", width = 4, height = 1.6, dpi = 300, family = "Arial")


col_pal = RColorBrewer::brewer.pal(9,"Purples")
go_use = c("negative regulation of phosphorylation","positive regulation of lipid metabolic process")
ego[["MOLs"]][["y"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,2), expand = c(0,0), breaks = c(0,1,2)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "mols_ego_y.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig5/", width = 4, height = 1.6, dpi = 300, family = "Arial")


col_pal = RColorBrewer::brewer.pal(9,"Blues")
go_use = c("regulation of cell projection assembly","oligodendrocyte differentiation")
ego[["MOLs"]][["z"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,2), expand = c(0,0), breaks = c(0,1,2)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "mols_ego_z.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig5/", width = 4, height = 1.6, dpi = 300, family = "Arial")


col_pal = RColorBrewer::brewer.pal(9,"Greens")
go_use = c("cellular response to hormone stimulus","regulation of endocytosis")
ego[["MOLs"]][["w"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(fill = col_pal[6], stat = "identity", color = "black", width = .6) +
  # scale_fill_gradient(low = col_pal[5], high = col_pal[7]) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 25), position = "right") +
  scale_x_continuous(limits = c(0,3), expand = c(0,0)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "none",
    # legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(size = 1),
    axis.text.x = element_text(size = 16, hjust = 0),
    axis.title.x = element_blank(),
    # axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
    plot.margin = margin(t=0,r=0,b=0,l=10)
  )
ggsave(filename = "mols_ego_w.eps", plot = last_plot(), device = "eps", path = "../figures/ExtDataFig5/", width = 4, height = 1.6, dpi = 300, family = "Arial")
