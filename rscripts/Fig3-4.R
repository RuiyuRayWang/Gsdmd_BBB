setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(tidyverse)

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})

# Fig3e: ORA - GO BP
ego.cell_subtype.up.dfs <- readRDS("../outs/functional_analysis/ora/ego_cell_subtype_up_dfs_LFC_1_PADJ_0.01.rds")

col_pal = RColorBrewer::brewer.pal(9,"Oranges")
go_use = c("response to molecule of bacterial origin","response to lipopolysaccharide","pattern recognition receptor signaling pathway",
           "regulation of innate immune response","cytokine-mediated signaling pathway","regulation of inflammatory response")
ego.cell_subtype.up.dfs[["BP"]][["Endo"]][["WT"]] |>
  dplyr::filter(Description %in% go_use) |>
  dplyr::mutate(Count = as.numeric(Count)) |>
  ggplot(aes(x = -log10(p.adjust), y = reorder(Description, p.adjust, decreasing = TRUE))) +
  geom_bar(aes(fill = Count), stat = "identity", color = "black", width = .75) +
  scale_fill_gradient(low = col_pal[5], high = col_pal[8], limits = c(20,40)) +
  scale_y_discrete(name = "", labels = function(x) str_wrap(x, width = 28)) +
  scale_x_continuous(limits = c(0,17), expand = c(0,0)) +
  ggtitle(NULL) +
  xlab(bquote(-log[10](p.adjust))) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 17), legend.text = element_text(size = 16),
    panel.background = element_rect(size = .8, color = "black"),
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 17, hjust = 1, margin = margin(t=0,r=2,b=0,l=0)),
  )
ggsave(filename = "endo_wt_ego_up.eps", plot = last_plot(), device = "eps", path = "../figures/Fig3/", width = 6, height = 4.5, dpi = 300, family = "Arial")
