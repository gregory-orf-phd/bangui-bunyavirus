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

L_metadata <-
  read_delim(
    file = "data/L_metadata.tsv",
    delim = "\t",
    na = c("NA")
  ) %>%
  select(taxon, accession, organism, serogroup)

# set up panel A ----------------------------------------------------------

tree1 <-
  (read.iqtree(
    file = "data/rooted.aligned.orthobunyavirus_L_AA.fasta.treefile"
  ) %>%
  groupClade(., MRCA(., "QLA47025|Tanga_virus", "QLA46988|Boraceia_virus")) %>%
  ggtree(
    .,
    ladderize = TRUE,
    mapping = aes(
      color = group
    )
  ) %<+%
  L_metadata) +
  scale_color_manual(
    values = c("gray80", "black"),
    guide = element_blank()
  ) +
  geom_rootedge(
    0.05,
    color = "gray80"
  ) +
  vexpand(
    ratio = 0.05,
    direction = -1
  ) +
  expand_limits(
    y = 730
  ) +
  geom_treescale(
    fontsize = 6,
    linesize = 1,
    offset = 10,
    x = 0.3,
    y = 200
  ) +
  geom_tippoint(
    mapping = aes(
      fill = serogroup
    ),
    shape = 21,
    size = 3,
    color = "black",
    na.rm = TRUE
  ) +
  guides(
    color = "none",
    fill = "none"
  ) +
  theme(
    legend.position = c(0.2, 0.5)
  )

tree3 <-
  read.iqtree(
    file = "data/rooted.aligned.orthobunyavirus_L_AA.fasta.treefile"
  ) %>%
  ggtree(
    .,
    ladderize = TRUE,
    color = "black"
  ) %<+%
  L_metadata

tree4 <-
  viewClade(
    tree_view = tree3, 
    node = MRCA(tree3, "QLA47025|Tanga_virus", "QLA46988|Boraceia_virus"),
    xmax_adjust = 0.2
  ) +
  geom_tiplab(
    mapping = aes(label = str_replace(str_replace_all(label, "_", " "), "orthobunya", "")),
    offset = 0.02
  ) +
  geom_nodelab(
    mapping = aes(label = UFboot),
    geom = "text",
    nudge_x = -0.055,
    nudge_y = 0.5
  ) +
  geom_tippoint(
    mapping = aes(
      fill = serogroup
    ),
    shape = 21,
    size = 5,
    color = "black",
    na.rm = TRUE
  ) +
  theme(
    legend.position = "bottom",
    plot.margin = unit(c(0, 0, 2, 0), "cm"),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  guides(
    fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5, title = "Serogroup")
  )

# import data for tanglegrams ---------------------------------------------

N23_ML_L <-
  read.iqtree(
    file = "data/rooted.aligned.serogroups_of_interest_L_codon.fasta.newick"
  ) %>%
  ggtree(
    .,
    ladderize = TRUE,
    color = "black"
  )
  # ) %<+%
  # L_metadata

N23_ML_M <-
  (read.iqtree(
    file = "data/rooted.aligned.serogroups_of_interest_M_codon.fasta.newick"
  ) %>%
  ggtree(
    .,
    ladderize = TRUE,
    color = "black"
  # ) %<+%
  # metadata) %>%
  )) %>%
  rotate(32)

N23_ML_S <-
  read.iqtree(
    file = "data/rooted.aligned.serogroups_of_interest_S_codon.fasta.newick"
  ) %>%
  ggtree(
    .,
    ladderize = TRUE,
    color = "black"
  # ) %<+%
  # metadata
  )

# first tangelgram for panel B --------------------------------------------

a1 <- 
  N23_ML_L$data

a2 <- 
  N23_ML_M$data %>%
  mutate(x, x = max(x) - x + max(a1$x) + 25)

aa1 <- 
  bind_rows(a1, a2) %>% 
  filter(isTip == "TRUE") %>%
  mutate(x_line = case_when(x > 19 ~ 19, x < 20 ~ 9))

tangle1 <-
  N23_ML_L + 
  geom_tree(data = a2) +
  ggnewscale::new_scale_fill() +
  geom_tiplab(
    aes(label = str_replace(str_remove(label, "\\|[0-9][0-9][0-9][0-9]"), "_", " "), subset = isTip),
    geom = "label", 
    data = a1, 
    hjust = 1, 
    vjust = 0.5,
    offset = 7, 
    size = 4, 
    align = TRUE,
    label.padding = unit(0.1, "lines"),
    label.size = 0
  ) +
  geom_tiplab(
    aes(label = str_replace(str_remove(label, "\\|[0-9][0-9][0-9][0-9]"), "_", " "), subset = isTip),
    geom = "label",
    data = a2, 
    hjust = 0, 
    vjust = 0.5,
    offset = -9, 
    size = 4,
    align = TRUE,
    label.padding = unit(0.1, "lines"),
    label.size = 0
  ) +
  # scale_color_manual(values = c("black", "red"), guide = "none") +
  geom_line(
    data = aa1,
    mapping = aes(
      x_line,
      y,
      group = label
    ),
    linetype = 2,
    linewidth = 0.5
  ) +
  ggtitle("L vs. M") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

# second tanglegram for panel B -------------------------------------------

b1 <- 
  N23_ML_L$data

b2 <- 
  N23_ML_S$data %>%
  mutate(x, x = max(x) - x + max(b1$x) + 25)

bb1 <- 
  bind_rows(b1, b2) %>% 
  filter(isTip == "TRUE") %>%
  mutate(x_line = case_when(x > 19 ~ 19, x < 20 ~ 9))

tangle2 <-
  N23_ML_L + 
  geom_tree(data = b2) +
  ggnewscale::new_scale_fill() +
  geom_tiplab(
    aes(label = str_replace(str_remove(label, "\\|[0-9][0-9][0-9][0-9]"), "_", " "), subset = isTip),
    geom = "label", 
    data = b1, 
    hjust = 1, 
    vjust = 0.5,
    offset = 7, 
    size = 4, 
    align = TRUE,
    label.padding = unit(0.1, "lines"),
    label.size = 0
  ) +
  geom_tiplab(
    aes(label = str_replace(str_remove(label, "\\|[0-9][0-9][0-9][0-9]"), "_", " "), subset = isTip),
    geom = "label",
    data = b2, 
    hjust = 0, 
    vjust = 0.5,
    offset = -10, 
    size = 4,
    align = TRUE,
    label.padding = unit(0.1, "lines"),
    label.size = 0
  ) +
  # scale_color_manual(values = c("black", "red"), guide = "none") +
  geom_line(
    data = bb1,
    mapping = aes(
      x_line,
      y,
      group = label
    ),
    linetype = 2,
    linewidth = 0.5
  ) +
  ggtitle("L vs. S") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

# third tanglegram for panel B --------------------------------------------

c1 <- 
  N23_ML_M$data

c2 <- 
  N23_ML_S$data %>%
  mutate(x, x = max(x) - x + max(c1$x) + 25)

cc1 <- 
  bind_rows(c1, c2) %>% 
  filter(isTip == "TRUE") %>%
  mutate(x_line = case_when(x > 19 ~ 19, x < 20 ~ 10))

tangle3 <-
  N23_ML_M + 
  geom_tree(data = c2) +
  ggnewscale::new_scale_fill() +
  geom_tiplab(
    aes(label = str_replace(str_remove(label, "\\|[0-9][0-9][0-9][0-9]"), "_", " "), subset = isTip),
    geom = "label", 
    data = c1, 
    hjust = 1, 
    vjust = 0.5,
    offset = 7.5, 
    size = 4, 
    align = TRUE,
    label.padding = unit(0.1, "lines"),
    label.size = 0
  ) +
  geom_tiplab(
    aes(label = str_replace(str_remove(label, "\\|[0-9][0-9][0-9][0-9]"), "_", " "), subset = isTip),
    geom = "label",
    data = c2, 
    hjust = 0, 
    vjust = 0.5,
    offset = -10, 
    size = 4,
    align = TRUE,
    label.padding = unit(0.1, "lines"),
    label.size = 0
  ) +
  # scale_color_manual(values = c("black", "red"), guide = "none") +
  geom_line(
    data = cc1,
    mapping = aes(
      x_line,
      y,
      group = label
    ),
    linetype = 2,
    linewidth = 0.5
  ) +
  ggtitle("M vs. S") +
  theme(plot.title = element_text(hjust = 0.5, size = 18))

# tangle3

# create panels -----------------------------------------------------------

Panel_A <-
  cowplot::plot_grid(
    tree1,
    tree4,
    labels = c("", ""),
    ncol = 2
  )

Panel_B <-
  cowplot::plot_grid(
    tangle1,
    tangle2,
    tangle3,
    labels = c("", "", ""),
    ncol = 3
  )

# create and save figure --------------------------------------------------

Fig_02 <-
  cowplot::plot_grid(
    Panel_A,
    Panel_B,
    labels = c("A", "B"),
    label_size = 24,
    ncol = 1,
    rel_heights = c(0.6, 0.4)
  )

# Fig_02

ggsave(
  plot = Fig_02,
  filename = "figures/Fig_02.svg",
  device = "svg",
  width = 7,
  height = 7,
  unit = "in",
  scale = 2
)

# supplementary figure S01 ------------------------------------------------

M_metadata <-
  read_delim(
    file = "data/M_metadata.tsv",
    delim = "\t",
    na = c("NA")
  )

M_tree1 <-
  (read.iqtree(
    file = "data/rooted.aligned.orthobunyavirus_M_AA.fasta.treefile"
  ) %>%
    groupClade(., MRCA(., "QLA47024|Tanga_virus", "ASY08210|Tacaiuma_orthobunyavirus")) %>%
    ggtree(
      .,
      ladderize = TRUE,
      mapping = aes(
        color = group
      )
    ) %<+%
    M_metadata) +
  scale_color_manual(
    values = c("gray80", "black"),
    guide = element_blank()
  ) +
  geom_rootedge(
    0.05,
    color = "gray80"
  ) +
  vexpand(
    ratio = 0.05,
    direction = -1
  ) +
  expand_limits(
    y = 950
  ) +
  geom_treescale(
    fontsize = 4,
    linesize = 1,
    offset = 10,
    x = 0.3,
    y = 200
  ) +
  geom_tippoint(
    mapping = aes(
      fill = serogroup
    ),
    shape = 21,
    size = 3,
    color = "black",
    na.rm = TRUE
  ) +
  guides(
    color = "none",
    fill = "none"
  ) +
  theme(
    legend.position = c(0.2, 0.5)
  )

M_tree3 <-
  read.iqtree(
    file = "data/rooted.aligned.orthobunyavirus_M_AA.fasta.treefile"
  ) %>%
  ggtree(
    .,
    ladderize = TRUE,
    color = "black"
  ) %<+%
  M_metadata

M_tree4 <-
  viewClade(
    tree_view = M_tree3, 
    node = MRCA(M_tree3, "QLA47024|Tanga_virus", "ASY08210|Tacaiuma_orthobunyavirus"),
    xmax_adjust = 0.75
  ) +
  geom_tiplab(
    mapping = aes(label = str_replace(str_replace_all(label, "_", " "), "orthobunya", "")),
    offset = 0.02
  ) +
  geom_nodelab(
    mapping = aes(label = UFboot),
    geom = "text",
    nudge_x = -0.095,
    nudge_y = 0.5
  ) +
  geom_tippoint(
    mapping = aes(
      fill = serogroup
    ),
    shape = 21,
    size = 5,
    color = "black",
    na.rm = TRUE
  ) +
  theme(
    legend.position = "bottom",
    plot.margin = unit(c(0, 0, 2, 0), "cm"),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  guides(
    fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5, title = "Serogroup")
  )

N_metadata <-
  read_delim(
    file = "data/N_metadata.tsv",
    delim = "\t",
    na = c("NA")
  )

N_tree1 <-
  (read.iqtree(
    file = "data/rooted.aligned.orthobunyavirus_N_AA.fasta.treefile"
  ) %>%
    groupClade(., MRCA(., "QLA47023|Tanga_virus", "ASY08208|Tacaiuma_orthobunyavirus")) %>%
    ggtree(
      .,
      ladderize = TRUE,
      mapping = aes(
        color = group
      )
    ) %<+%
    N_metadata) +
  scale_color_manual(
    values = c("gray80", "black"),
    guide = element_blank()
  ) +
  geom_rootedge(
    0.05,
    color = "gray80"
  ) +
  vexpand(
    ratio = 0.05,
    direction = -1
  ) +
  expand_limits(
    y = 1450
  ) +
  geom_treescale(
    fontsize = 4,
    linesize = 1,
    offset = 10,
    x = 0.3,
    y = 300
  ) +
  geom_tippoint(
    mapping = aes(
      fill = serogroup
    ),
    shape = 21,
    size = 3,
    color = "black",
    na.rm = TRUE
  ) +
  guides(
    color = "none",
    fill = "none"
  ) +
  theme(
    legend.position = c(0.2, 0.5)
  )

N_tree3 <-
  read.iqtree(
    file = "data/rooted.aligned.orthobunyavirus_N_AA.fasta.treefile"
  ) %>%
  ggtree(
    .,
    ladderize = TRUE,
    color = "black"
  ) %<+%
  N_metadata

N_tree4 <-
  viewClade(
    tree_view = N_tree3, 
    node = MRCA(N_tree3, "QLA47023|Tanga_virus", "ASY08208|Tacaiuma_orthobunyavirus"),
    xmax_adjust = 0.5
  ) +
  geom_tiplab(
    mapping = aes(label = str_replace(str_replace_all(label, "_", " "), "orthobunya", "")),
    offset = 0.02
  ) +
  geom_nodelab(
    mapping = aes(label = UFboot),
    geom = "text",
    nudge_x = -0.075,
    nudge_y = 0.5
  ) +
  geom_tippoint(
    mapping = aes(
      fill = serogroup
    ),
    shape = 21,
    size = 5,
    color = "black",
    na.rm = TRUE
  ) +
  theme(
    legend.position = "top",
    plot.margin = unit(c(0, 0, 2, 0), "cm"),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  guides(
    fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5, title = "Serogroup")
  )

Supp_Panel_A <-
  cowplot::plot_grid(
    M_tree1,
    M_tree4,
    labels = c("", ""),
    ncol = 2
  )

Supp_Panel_B <-
  cowplot::plot_grid(
    N_tree1,
    N_tree4,
    labels = c("", ""),
    ncol = 2
  )

Fig_S01 <-
  cowplot::plot_grid(
    Supp_Panel_A,
    Supp_Panel_B,
    labels = c("A", "B"),
    label_size = 24,
    hjust = 0,
    ncol = 1,
    rel_heights = c(0.5, 0.5)
  )

# Fig_S01

ggsave(
  plot = Fig_S01,
  filename = "figures/Fig_S01.svg",
  device = "svg",
  width = 7,
  height = 9,
  unit = "in",
  scale = 2
)

# supplementary figure S02 ------------------------------------------------

# No_Nyangole_M_tree1 <-
#   (read.iqtree(
#     file = "data/rooted.no_nyangole.aligned.orthobunyavirus_M_AA.fasta.treefile"
#   ) %>%
#     groupClade(., MRCA(., "QLA47024|Tanga_virus", "ASY08210|Tacaiuma_orthobunyavirus")) %>%
#     ggtree(
#       .,
#       ladderize = TRUE,
#       mapping = aes(
#         color = group
#       )
#     ) %<+%
#     M_metadata) +
#   scale_color_manual(
#     values = c("gray80", "black"),
#     guide = element_blank()
#   ) +
#   geom_rootedge(
#     0.05,
#     color = "gray80"
#   ) +
#   vexpand(
#     ratio = 0.05,
#     direction = -1
#   ) +
#   expand_limits(
#     y = 950
#   ) +
#   geom_treescale(
#     fontsize = 4,
#     linesize = 1,
#     offset = 10,
#     x = 0.3,
#     y = 200
#   ) +
#   geom_tippoint(
#     mapping = aes(
#       fill = serogroup
#     ),
#     shape = 21,
#     size = 3,
#     color = "black",
#     na.rm = TRUE
#   ) +
#   guides(
#     color = "none",
#     fill = "none"
#   ) +
#   theme(
#     legend.position = c(0.2, 0.5)
#   )

No_Nyangole_M_tree3 <-
  read.iqtree(
    file = "data/rooted.no_nyangole.aligned.orthobunyavirus_M_AA.fasta.treefile"
  ) %>%
  ggtree(
    .,
    ladderize = TRUE,
    color = "black"
  ) %<+%
  M_metadata

No_Nyangole_M_tree4 <-
  viewClade(
    tree_view = No_Nyangole_M_tree3, 
    node = MRCA(No_Nyangole_M_tree3, "QLA47024|Tanga_virus", "ASY08210|Tacaiuma_orthobunyavirus"),
    xmax_adjust = 0.75
  ) +
  geom_tiplab(
    mapping = aes(label = str_replace(str_replace_all(label, "_", " "), "orthobunya", "")),
    offset = 0.02
  ) +
  geom_nodelab(
    mapping = aes(label = UFboot),
    geom = "text",
    nudge_x = -0.125,
    nudge_y = 0.5
  ) +
  geom_tippoint(
    mapping = aes(
      fill = serogroup
    ),
    shape = 21,
    size = 5,
    color = "black",
    na.rm = TRUE
  ) +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 2, 0), "cm"),
    legend.background = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold")
  ) +
  guides(
    fill = guide_legend(ncol = 2, title.position = "top", title.hjust = 0.5, title = "Serogroup")
  )

Fig_S01_inset <-
  No_Nyangole_M_tree4

# Fig_S01_inset

ggsave(
  plot = Fig_S01_inset,
  filename = "figures/Fig_S01_inset.svg",
  device = "svg",
  width = 3.5,
  height = 2.5,
  unit = "in",
  scale = 2
)
