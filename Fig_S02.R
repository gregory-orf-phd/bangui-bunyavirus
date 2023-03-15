# Gregory S. Orf, Ph.D.
# Virus Discovery Group, Abbott Diagnostics Division, Abbott Laboratories, U.S.A.
# Purifying selection decreases the potential for Bangui orthobunyavirus outbreaks in humans
# DOI: https://doi.org/10.1093/ve/vead018

# load libraries ----------------------------------------------------------

library(tidyverse)
library(treeio)
library(ggtree)
library(tidytree)
library(cowplot)

# import data -------------------------------------------------------------

metadata <-
  read_delim(
    file = "data/orthobunyavirus.metadata.2.tsv",
    delim = "\t",
    na = c("NA")
  )

N23_ML_L_codon <-
  read.iqtree(
    file = "data/rooted.aligned.serogroups_of_interest_L_codon.fasta.newick"
  ) %>%
    groupOTU(., c("Bangui|1970", "Bangui|2017", "Nyangole|2013")) %>%
    ggtree(
      .,
      ladderize = TRUE,
      color = "black"
    ) %<+%
    metadata

N23_ML_M_codon <-
  read.iqtree(
    file = "data/rooted.aligned.serogroups_of_interest_M_codon.fasta.newick"
  ) %>%
    groupOTU(., c("Bangui|1970", "Bangui|2017", "Nyangole|2013")) %>%
    ggtree(
      .,
      ladderize = TRUE,
      color = "black"
    ) %<+%
    metadata

N23_ML_S_codon <-
  read.iqtree(
    file = "data/rooted.aligned.serogroups_of_interest_S_codon.fasta.newick"
  ) %>%
    groupOTU(., c("Bangui|1970", "Bangui|2017", "Nyangole|2013")) %>%
    ggtree(
      .,
      ladderize = TRUE,
      color = "black"
    ) %<+%
    metadata

# configure trees ---------------------------------------------------------

Panel_A <-
  N23_ML_L_codon +
  geom_tiplab(
    mapping = aes(label = label, subset = isTip, color = group),
    geom = "label", 
    hjust = 0, 
    vjust = 0.5,
    size = 4, 
    label.padding = unit(0.1, "lines"),
    label.size = 0
    ) +
  scale_color_manual(
    values = c("black", "red"),
    guide = NULL
  ) +
  geom_rootedge(
    0.05,
    color = "black"
  ) +
  scale_x_continuous(
    limits = c(-0.05, 2.5),
    expand = c(0, 0)
  ) +
  geom_treescale(
    fontsize = 4,
    linesize = 1,
    offset = 0.5,
    x = 0.1,
    y = 15
  ) +
  geom_nodelab(
    mapping = aes(label = UFboot),
    geom = "text",
    hjust = 1,
    vjust = 0.5,
    nudge_x = -0.05,
    nudge_y = 0.5
  ) +
  # geom_tippoint(
  #   mapping = aes(
  #     fill = serogroup
  #   ),
  #   shape = 21,
  #   size = 5,
  #   color = "black",
  #   na.rm = TRUE
  # ) +
  theme(
    legend.position = "bottom",
    plot.margin = unit(c(0, 0, 2, 0), "cm"),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, size = rel(1.2), face = "bold")
  ) +
  labs(caption = "L Polyprotein, Codon Alignment")
  # guides(
  #   fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5, title = "Serogroup")
  # )

Panel_B <-
  N23_ML_M_codon +
  geom_tiplab(
    mapping = aes(label = label, subset = isTip, color = group),
    geom = "label", 
    hjust = 0, 
    vjust = 0.5,
    size = 4, 
    label.padding = unit(0.1, "lines"),
    label.size = 0
  ) +
  scale_color_manual(
    values = c("black", "red"),
    guide = NULL
  ) +
  geom_rootedge(
    0.05,
    color = "black"
  ) +
  geom_treescale(
    fontsize = 4,
    linesize = 1,
    offset = 0.5,
    x = 0.1,
    y = 15
  ) +
  scale_x_continuous(
    limits = c(-0.05, 3.5),
    expand = c(0, 0)
  ) +
  geom_nodelab(
    mapping = aes(label = UFboot),
    geom = "text",
    hjust = 1,
    vjust = 0.5,
    nudge_x = -0.05,
    nudge_y = 0.5
  ) +
  # geom_tippoint(
  #   mapping = aes(
  #     fill = serogroup
  #   ),
  #   shape = 21,
  #   size = 5,
  #   color = "black",
  #   na.rm = TRUE
  # ) +
  theme(
    legend.position = "bottom",
    plot.margin = unit(c(0, 0, 2, 0), "cm"),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, size = rel(1.2), face = "bold")
  ) +
  labs(caption = "M Polyprotein, Codon Alignment")
  # guides(
  #   fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5, title = "Serogroup")
  # )
  
Panel_C <-
  N23_ML_S_codon +
  geom_tiplab(
    mapping = aes(label = label, subset = isTip, color = group),
    geom = "label", 
    hjust = 0, 
    vjust = 0.5,
    size = 4, 
    label.padding = unit(0.1, "lines"),
    label.size = 0
  ) +
  scale_color_manual(
    values = c("black", "red"),
    guide = NULL
  ) +
  geom_rootedge(
    0.05,
    color = "black"
  ) +
  geom_treescale(
    fontsize = 4,
    linesize = 1,
    offset = 0.5,
    x = 0.1,
    y = 15
  ) +
  scale_x_continuous(
    limits = c(-0.05, 4.75),
    expand = c(0, 0)
  ) +
  geom_nodelab(
    mapping = aes(label = UFboot),
    geom = "text",
    hjust = 1,
    vjust = 0.5,
    nudge_x = -0.05,
    nudge_y = 0.5
  ) +
  # geom_tippoint(
  #   mapping = aes(
  #     fill = serogroup
  #   ),
  #   shape = 21,
  #   size = 5,
  #   color = "black",
  #   na.rm = TRUE
  # ) +
  theme(
    legend.position = "bottom",
    plot.margin = unit(c(0, 0, 2, 0), "cm"),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5),
    plot.caption = element_text(hjust = 0.5, size = rel(1.2), face = "bold")
  ) +
  labs(caption = "S Protein, Codon Alignment")
  # guides(
  #   fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5, title = "Serogroup")
  # )

# make figure -------------------------------------------------------------

Fig_S02 <-
  cowplot::plot_grid(
    Panel_A,
    Panel_B,
    Panel_C,
    ncol = 3,
    labels = c("A", "B", "C"),
    label_size = 24
  )

ggsave(
  filename = "figures/Fig_S02.svg",
  plot = Fig_S02,
  device = "svg",
  width = 6.5,
  height = 3,
  unit = "in",
  scale = 2
)

