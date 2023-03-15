# Gregory S. Orf, Ph.D.
# Virus Discovery Group, Abbott Diagnostics Division, Abbott Laboratories, U.S.A.
# Purifying selection decreases the potential for Bangui orthobunyavirus outbreaks in humans
# DOI: https://doi.org/10.1093/ve/vead018

# load libraries ----------------------------------------------------------

library(tidyverse)
library(cowplot)
library(ggpubr)

# import data -------------------------------------------------------------

L_RTT <-
  read_delim(
    "data/L_tempest_RTT_regression_data.tsv",
    delim = "\t"
  )
  
M_RTT <-
  read_delim(
    "data/M_tempest_RTT_regression_data.tsv",
    delim = "\t"
  )

S_RTT <-
  read_delim(
    "data/S_tempest_RTT_regression_data.tsv",
    delim = "\t"
  )

joined_data <-
  bind_rows(
    L_RTT %>% select(-residual),
    M_RTT %>% select(-residual),
    S_RTT %>% select(-residual),
    .id = "segment"
  ) %>%
  mutate(
    segment = str_replace(segment, "1", "L Segment"),
    segment = str_replace(segment, "2", "M Segment"),
    segment = str_replace(segment, "3", "S Segment")
  )
  
# create plots ------------------------------------------------------------

plots <-
  ggplot(
    data = joined_data,
    mapping = aes(x = date, y = distance)
  ) +
  theme_bw() +
  geom_point(
    mapping = aes(color = segment),
    size = 2,
    show.legend = FALSE
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "black"
  ) +
  scale_x_continuous(
    name = "Sampling Date",
    limits = c(1940, 2020)
  ) +
  scale_y_continuous(
    name = "Root-to-Tip Distance",
  ) +
  stat_regline_equation(
    label.x = 1985,
    label.y = 3.4,
    aes(label = after_stat(eq.label)),
    size = 3
  ) +
  stat_regline_equation(
    label.x = 1985,
    label.y = 3.2,
    aes(label = ..rr.label..),
    size = 3
  ) +
  facet_wrap(
    ~segment,
    nrow = 1
    # ,scales = "free"
  ) +
  theme(
    strip.text = element_text(size = 10),
    strip.background = element_blank(),
    strip.placement = "outside", 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  
# save figure -------------------------------------------------------------

ggsave(
  plot = plots,
  filename = "figures/Fig_S04.svg",
  device = "svg",
  height = 450,
  width = 1200,
  scale = 2,
  unit = "px"
)
