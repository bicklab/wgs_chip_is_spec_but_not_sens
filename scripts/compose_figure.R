library(patchwork)
library(ggplot2)
library(readr)

p1a <- read_rds("figures/figure_data/dts_fig1a.rds") + labs(x = "VAF by DTS", fill = "min(VAF)")
p1b <- read_rds("figures/figure_data/dts_fig1b.rds") + labs(y = "WGS PPV", fill = "min(VAF)")
p1c <- read_rds("sim_data/p1c_2025-07-31.rds") + labs(tag = "C", fill = "min(VAF)")
p1d <- read_rds("sim_data/p1d_2025-07-31.rds") + labs(tag = "D", fill = "min(VAF)")
p1e <- read_rds("sim_data/p1e_2025-07-31.rds") + labs(tag = "E", fill = "min(VAF)")
p1f <- read_rds("sim_data/p1f_2025-07-31.rds") + labs(tag = "F", fill = "min(VAF)")

(p1a + p1b) /
   (p1c + p1d) /
   (p1e + p1f) +
   plot_layout(guides = "collect") &
   theme(legend.position = "bottom")

ggsave(filename = "figures/fig1_JCI.pdf", height = 8, width = 8)



p1c <- read_rds("p1c_pos_ctrl.rds")
p1d <- read_rds("p1d_pos_ctrl.rds")
p1e <- read_rds("p1e_pos_ctrl.rds")
p1f <- read_rds("p1f_pos_ctrl.rds")

(p1c + p1d) /
   (p1e + p1f) +
   plot_layout(guides = "collect") &
   theme(legend.position = "bottom")

ggsave(filename = "fig1_JCI_pos_ctrl.pdf", height = 8, width = 8)
