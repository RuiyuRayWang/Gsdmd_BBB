library(tidyverse)
library(ggrepel)

suppressMessages({
  # extrafont::font_import()
  extrafont::loadfonts(device="postscript")
})
tab <- read.csv("../miscs/saturation_maxc_3000_s_10_r_50_20230429.csv")
tab$Group <- factor(tab$Group, levels = c("Ex-Neu.","Inh-Neu.","Astro","Oligo","OPCs","Micro","Endo"))


# ExtFig4 d
# pdf(
#   file = "../figures/ExtDataFig13/ExtDataFig13D.pdf",
#   width = 4, height = 3
# )
p = tab |>
  group_by(Cells, Group) |>
  summarise(mu = mean(Gene)) |>
  group_by(Group) |>
  mutate(label = if_else(Cells == max(Cells), as.character(Group), NA_character_)) |>
  ggplot(aes(x = Cells, y = mu)) +
  geom_line(aes(color = Group)) +
  geom_label_repel(aes(label = label), nudge_x = 1, na.rm = TRUE, force = 10) +
  scale_color_discrete(name = "Cell type") +
  scale_x_continuous(breaks = seq(0,3000,600), name = "Number of cells considered") +
  scale_y_continuous(name = "Cumulative number of\ngenes recovered") +
  theme_linedraw() +
  theme(
    legend.position = "none",
    # axis.text = element_text(family = "Arial"),
    # axis.title = element_text(family = "Arial"),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .2, color = "gray"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 10),
    axis.text.y = element_text(angle = 90, size = 10, hjust = .5)
  )
ggexport(p, filename = "../figures/ExtendedFig/ExtDataFig4D.pdf", family = "Arial", width = 4, height = 3)
