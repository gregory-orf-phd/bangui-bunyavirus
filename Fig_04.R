# Gregory S. Orf, Ph.D.
# Virus Discovery Group, Abbott Diagnostics Division, Abbott Laboratories, U.S.A.
# Purifying selection decreases the potential for Bangui orthobunyavirus outbreaks in humans
# DOI: https://doi.org/10.1093/ve/vead018

# load libraries ----------------------------------------------------------

library(tidyverse)
library(treeio)
library(ggtree)
library(tidytree)
library(hexbin)
library(cowplot)

# import data -------------------------------------------------------------

L_heterochronous_MCC_tree <-
  read.beast(
    "data/rooted.mcc.L_codon_heterochronous.treefile"
  )

L_isochronous_MCC_tree <- 
  read.beast(
    "data/rooted.mcc.L_codon_isochronous.treefile"
  )

M_heterochronous_MCC_tree <- 
  read.beast(
    "data/rooted.mcc.M_codon_heterochronous.treefile"
  )

M_isochronous_MCC_tree <- 
  read.beast(
    "data/rooted.mcc.M_codon_isochronous.treefile"
  )

S_heterochronous_MCC_tree <- 
  read.beast(
    "data/rooted.mcc.S_codon_heterochronous.treefile"
  )

S_isochronous_MCC_tree <- 
  read.beast(
    "data/rooted.mcc.S_codon_isochronous.treefile"
  )

L_codeml_m0 <- 
  read.codeml(
    "data/rst_L_m0.log",
    "data/mlc_L_m0.log"
  )

L_codeml_m2 <- 
  read.codeml(
    "data/rst_L_m2.log",
    "data/mlc_L_m2.log"
  )

L_codeml_m8 <- 
  read.codeml(
    "data/rst_L_m8.log",
    "data/mlc_L_m8.log"
  )

M_codeml_m0 <- 
  read.codeml(
    "data/rst_M_m0.log",
    "data/mlc_M_m0.log"
  )

M_codeml_m2 <- 
  read.codeml(
    "data/rst_M_m2.log",
    "data/mlc_M_m2.log"
  )

M_codeml_m8 <- 
  read.codeml(
    "data/rst_M_m8.log",
    "data/mlc_M_m8.log"
  )

S_codeml_m0 <- 
  read.codeml(
    "data/rst_S_m0.log",
    "data/mlc_S_m0.log"
  )

S_codeml_m2 <- 
  read.codeml(
    "data/rst_S_m2.log",
    "data/mlc_S_m2.log"
  )

S_codeml_m8 <- 
  read.codeml(
    "data/rst_S_m8.log",
    "data/mlc_S_m8.log"
  )

L_merged_m0 <- 
  merge_tree(L_heterochronous_MCC_tree, L_codeml_m0)

L_merged_m2 <- 
  merge_tree(L_heterochronous_MCC_tree, L_codeml_m2)

L_merged_m8 <- 
  merge_tree(L_heterochronous_MCC_tree, L_codeml_m8)

M_merged_m0 <- 
  merge_tree(M_heterochronous_MCC_tree, M_codeml_m0)

M_merged_m2 <- 
  merge_tree(M_heterochronous_MCC_tree, M_codeml_m2)

M_merged_m8 <- 
  merge_tree(M_heterochronous_MCC_tree, M_codeml_m8)

S_merged_m0 <- 
  merge_tree(S_heterochronous_MCC_tree, S_codeml_m0)

S_merged_m2 <- 
  merge_tree(S_heterochronous_MCC_tree, S_codeml_m2)

S_merged_m8 <- 
  merge_tree(S_heterochronous_MCC_tree, S_codeml_m8)

# manipulate data ---------------------------------------------------------

df_L_m0 <-
  L_merged_m0 %>%
  tidytree::as_tibble() %>%
  select(dN_vs_dS, dN, dS, rate) %>%
  filter(dN_vs_dS >= 0 & dN_vs_dS <= 1.5) %>%
  pivot_longer(cols = c("dN_vs_dS", "dN", "dS"), names_to = "type", values_to = "value") %>%
  drop_na() %>%
  mutate(segment = "L Segment") %>%
  mutate(site_model = "M0")

df_L_m2 <-
  L_merged_m2 %>%
  tidytree::as_tibble() %>%
  select(dN_vs_dS, dN, dS, rate) %>%
  filter(dN_vs_dS >= 0 & dN_vs_dS <= 1.5) %>%
  pivot_longer(cols = c("dN_vs_dS", "dN", "dS"), names_to = "type", values_to = "value") %>%
  drop_na() %>%
  mutate(segment = "L Segment") %>%
  mutate(site_model = "M2")

df_L_m8 <-
  L_merged_m8 %>%
  tidytree::as_tibble() %>%
  select(dN_vs_dS, dN, dS, rate) %>%
  filter(dN_vs_dS >= 0 & dN_vs_dS <= 1.5) %>%
  pivot_longer(cols = c("dN_vs_dS", "dN", "dS"), names_to = "type", values_to = "value") %>%
  drop_na() %>%
  mutate(segment = "L Segment") %>%
  mutate(site_model = "M8")

df_M_m0 <-
  M_merged_m0 %>%
  tidytree::as_tibble() %>%
  select(dN_vs_dS, dN, dS, rate) %>%
  filter(dN_vs_dS >= 0 & dN_vs_dS <= 1.5) %>%
  pivot_longer(cols = c("dN_vs_dS", "dN", "dS"), names_to = "type", values_to = "value") %>%
  drop_na() %>%
  mutate(segment = "M Segment") %>%
  mutate(site_model = "M0")

df_M_m2 <-
  M_merged_m2 %>%
  tidytree::as_tibble() %>%
  select(dN_vs_dS, dN, dS, rate) %>%
  filter(dN_vs_dS >= 0 & dN_vs_dS <= 1.5) %>%
  pivot_longer(cols = c("dN_vs_dS", "dN", "dS"), names_to = "type", values_to = "value") %>%
  drop_na() %>%
  mutate(segment = "M Segment") %>%
  mutate(site_model = "M2")

df_M_m8 <-
  M_merged_m8 %>%
  tidytree::as_tibble() %>%
  select(dN_vs_dS, dN, dS, rate) %>%
  filter(dN_vs_dS >= 0 & dN_vs_dS <= 1.5) %>%
  pivot_longer(cols = c("dN_vs_dS", "dN", "dS"), names_to = "type", values_to = "value") %>%
  drop_na() %>%
  mutate(segment = "M Segment") %>%
  mutate(site_model = "M8")

df_S_m0 <-
  S_merged_m0 %>%
  tidytree::as_tibble() %>%
  select(dN_vs_dS, dN, dS, rate) %>%
  filter(dN_vs_dS >= 0 & dN_vs_dS <= 1.5) %>%
  pivot_longer(cols = c("dN_vs_dS", "dN", "dS"), names_to = "type", values_to = "value") %>%
  drop_na() %>%
  mutate(segment = "S Segment") %>%
  mutate(site_model = "M0")

df_S_m2 <-
  S_merged_m2 %>%
  tidytree::as_tibble() %>%
  select(dN_vs_dS, dN, dS, rate) %>%
  filter(dN_vs_dS >= 0 & dN_vs_dS <= 1.5) %>%
  pivot_longer(cols = c("dN_vs_dS", "dN", "dS"), names_to = "type", values_to = "value") %>%
  drop_na() %>%
  mutate(segment = "S Segment") %>%
  mutate(site_model = "M2")

df_S_m8 <-
  S_merged_m8 %>%
  tidytree::as_tibble() %>%
  select(dN_vs_dS, dN, dS, rate) %>%
  filter(dN_vs_dS >= 0 & dN_vs_dS <= 1.5) %>%
  pivot_longer(cols = c("dN_vs_dS", "dN", "dS"), names_to = "type", values_to = "value") %>%
  drop_na() %>%
  mutate(segment = "S Segment") %>%
  mutate(site_model = "M8")


df_LMS_m028 <-
  bind_rows(
    df_L_m0,
    df_L_m2,
    df_L_m8,
    df_M_m0,
    df_M_m2,
    df_M_m8,
    df_S_m0,
    df_S_m2,
    df_S_m8
  ) %>%
  mutate(type = replace(type, type %in% c("dN_vs_dS"), "dN/dS")) %>%
  mutate(type = fct_relevel(type, c("dN/dS", "dN", "dS")))

# Panel A -----------------------------------------------------------------

Panel_A <-
  ggplot(
    data = df_LMS_m028,
    mapping = aes(
      x = rate / 10^-3,
      y = value,
      # fill = site_model
      color = site_model
    )
  ) +
  theme_bw() +
  scale_x_continuous(
    # trans = rate * 10^-3,
    name = expression("rate ("*"\u00D7"*10^-3*")")
  ) +
  scale_y_continuous(
    position = "right"
  ) +
  geom_point(
    # color = "black",
    size = 2,
    shape = 16
    # shape = 21
  ) +
  geom_hline(
    data = filter(df_LMS_m028, type == "dN/dS"),
    mapping = aes(yintercept = 1),
    color = "black", 
    linetype = 2
  ) +
  scale_color_manual(
    name = "Site Model",
    values = c(
      "M0" = "black",
      "M2" = "red",
      "M8" = "blue"
    )
  ) +
  facet_grid(
    rows = vars(type),
    cols = vars(segment),
    scale = "free",
    switch = "y"
  ) +
  theme(
    legend.position = "top",
    plot.margin = unit(c(0, 0, 1, 1), "cm"),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

# ggdraw(Panel_A)

# tree options ------------------------------------------------------------

tree_labels <-
  list(
    theme_tree2(),
    # theme_tree(),
    # geom_treescale(
    #   x = 0,
    #   y = 23,
    #   width = 0.2,
    #   fontsize = 4,
    #   offset = 0.5,
    #   linesize = 1,
    #   color = "black"
    # ),
    geom_rootedge(
      rootedge = 0.05,
      size = 1,
      color = "grey70"
    ),
    vexpand(
      ratio = 0.05,
      direction = -1
    ),
    hexpand(
      ratio = 0.15,
      direction = 1
    ),
    expand_limits(
      y = 24
    ),
    geom_tiplab(
      color = "black",
      size = 2
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
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(hjust = 0.5),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

# Panel B -----------------------------------------------------------------

Panel_B <-
  ggtree(
    L_isochronous_MCC_tree,
    color = "grey70",
    size = 1
  ) +
  scale_x_continuous(
    name = expression("Substitutions/site"),
    limits = c(-0.1, 2.2),
    breaks = seq(from = 0, to = 2, by = 0.5),
    labels = seq(from = 0, to = 2, by = 0.5)
  ) +
  ggtitle("L Segment") +
  tree_labels +
  tree_theme

# ggdraw(Panel_B)

# Panel C -----------------------------------------------------------------

Panel_C <-
  ggtree(
    M_isochronous_MCC_tree,
    color = "grey70",
    size = 1
  ) +
  scale_x_continuous(
    name = expression("Substitutions/site"),
    limits = c(-0.1, 3.25),
    breaks = seq(from = 0, to = 3.25, by = 0.5),
    labels = seq(from = 0, to = 3.25, by = 0.5)
  ) +
  ggtitle("M Segment") +
  tree_labels +
  tree_theme

# ggdraw(Panel_C)

# Panel D -----------------------------------------------------------------

Panel_D <-
  ggtree(
    S_isochronous_MCC_tree, 
    color = "grey70",
    size = 1
  ) +
  scale_x_continuous(
    name = expression("Substitutions/site"),
    limits = c(-0.1, 3.6),
    breaks = seq(from = 0, to = 3.6, by = 0.5),
    labels = seq(from = 0, to = 3.6, by = 0.5)
  ) +
  ggtitle("S Segment") +
  tree_labels +
  tree_theme

# ggdraw(Panel_D)

# arrange and save figure -------------------------------------------------

Fig_04 <-
  cowplot::plot_grid(
    Panel_A,
    Panel_B,
    Panel_C,
    Panel_D,
    nrow = 2,
    align = c("hv"),
    labels = c("A", "B", "C", "D"),
    label_size = 24
  )

# ggdraw(Fig_04)

ggsave(
  filename = "figures/Fig_04.png",
    plot = Fig_04,
    device = "png",
    dpi = 1200,
    height = 4,
    width = 5,
    bg = "white",
    units = c("in"),
    limitsize = FALSE,
    scale = 2
)
