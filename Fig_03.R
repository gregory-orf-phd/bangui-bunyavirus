# Gregory S. Orf, Ph.D.
# Virus Discovery Group, Abbott Diagnostics Division, Abbott Laboratories, U.S.A.
# Purifying selection decreases the potential for Bangui orthobunyavirus outbreaks in humans
# DOI: https://doi.org/10.1093/ve/vead018

# load libraries ----------------------------------------------------------

library(tidyverse)
library(treeio)
library(ggtree)
library(tidytree)
library(lemon)
library(hexbin)
library(cowplot)
library(svglite)

# import data -------------------------------------------------------------

L_MCC_tree <- 
  read.beast(
    "data/urlc.yule.clustalo.combined_L_cds.treefile"
  )

M_MCC_tree <- 
  read.beast(
    "data/urlc.yule.clustalo.combined_M_cds.treefile"
  )

S_MCC_tree <- 
  read.beast(
    "data/urlc.yule.clustalo.combined_S_cds.treefile"
  )

L_MCC_tree@phylo$tip.label <-
  str_replace(L_MCC_tree@phylo$tip.label, "\\|.*?\\|", "|")

L_MCC_tree@phylo$tip.label <-
  str_replace_all(L_MCC_tree@phylo$tip.label, "_", " ")

M_MCC_tree@phylo$tip.label <-
  str_replace(M_MCC_tree@phylo$tip.label, "\\|.*?\\|", "|")

M_MCC_tree@phylo$tip.label <-
  str_replace_all(M_MCC_tree@phylo$tip.label, "_", " ")

S_MCC_tree@phylo$tip.label <-
  str_replace(S_MCC_tree@phylo$tip.label, "\\|.*?\\|", "|")

S_MCC_tree@phylo$tip.label <-
  str_replace_all(S_MCC_tree@phylo$tip.label, "_", " ")

# tree options ------------------------------------------------------------

tree_labels <-
  list(
    theme_tree2(),
    vexpand(
      ratio = 0.025,
      direction = -1
    ),
    expand_limits(
      y = 24
    ),
    # hexpand(
    #   ratio = 0.1,
    #   direction = 1
    # ),
    geom_rootedge(
      rootedge = 2000,
      size = 1,
      color = "grey70"
    ),
    geom_tiplab(
      color = "black",
      size = 2
    ),
    geom_range(
      range = "height_0.95_HPD",
      color = "black",
      alpha = 0.5,
      size = 2
    ), 
    scale_color_gradient(
      low = "red",
      high = "blue",
      name = "Posterior",
      limits = c(0, 1)
    ),
    coord_capped_cart(
      bottom = "right"
    ),
    guides(
      color = guide_colorbar(
        frame.colour = "black",
        ticks.colour = "black",
        ticks = FALSE
      )
    )
  )

tree_theme <-
  theme(
    axis.text.x = element_text(color = "black", hjust = 0.5, vjust = 1),
    axis.title.x = element_text(color = "black", hjust = 0.5, vjust = -1),
    # legend.box = "horizontal",
    legend.position = c(0.2, 0.8),
    # legend.direction = "vertical",
    legend.key.height = unit(0.75, "line"),
    # legend.title = element_text(size = 12),
    legend.title.align = 0.5,
    # legend.box.just = "center",
    # legend.text = element_text(size = 10),
    panel.background = element_rect(fill = "transparent", color = NA),
    axis.line.x = element_line(),
    axis.ticks.x = element_line(color = "black"),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

# Panel A -----------------------------------------------------------------

Panel_A <-
  revts(
      ggtree(
      L_MCC_tree,
      mapping = aes(color = posterior),
      size = 1
    )
  ) +
  scale_x_continuous(
    name = expression("Time before present (years"~"\u00D7"*10^3*")"),
    limits = c(-35000, 12000),
    breaks = seq(from = -35000, to = 0, by = 5000),
    labels = seq(from = 35, to = 0, by = -5),
    expand = c(0, 0)
  ) +
  ggtitle("L Segment") +
  tree_labels +
  tree_theme

ggdraw(Panel_A)

# Panel B -----------------------------------------------------------------

Panel_B <-
  revts(  
    ggtree(
      M_MCC_tree,
      mapping = aes(color = posterior),
      size = 1
    ) 
  ) +
  scale_x_continuous(
    name = expression("Time before present (years"~"\u00D7"*10^3*")"),
    limits = c(-40000, 15000),
    breaks = seq(from = -40000, to = 0, by = 5000),
    labels = seq(from = 40, to = 0, by = -5),
    expand = c(0, 0)
  ) +
  ggtitle("M Segment") +
  tree_labels +
  tree_theme

# ggdraw(Panel_B)

# Panel C -----------------------------------------------------------------

Panel_C <-
  revts(
    ggtree(
      S_MCC_tree, 
      mapping = aes(color = posterior), 
      size = 1
    ) 
  ) +
  scale_x_continuous(
    name = expression("Time before present (years"~"\u00D7"*10^3*")"),
    limits = c(-60000, 20000),
    breaks = seq(from = -60000, to = 0, by = 5000),
    labels = seq(from = 60, to = 0, by = -5),
    expand = c(0, 0)
  ) +
  ggtitle("S Segment") +
  tree_labels +
  tree_theme

# ggdraw(Panel_C)

# arrange and save figure -------------------------------------------------

Fig_03 <-
  cowplot::plot_grid(
    Panel_A,
    Panel_B,
    Panel_C,
    NULL,
    nrow = 2,
    align = c("hv"),
    labels = c("A", "C", "B", "D"),
    label_size = 32
  )

# ggdraw(Fig_03)

ggsave(
  filename = "figures/Fig_03.png",
    plot = Fig_03,
    device = "png",
    dpi = 1200,
    height = 7,
    width = 7,
    bg = "transparent",
    units = c("in"),
    limitsize = FALSE,
    scale = 1.5
)
