library(glue)
library(tidyverse)
library(broom)
library(survival)
library(scales)
library(colorspace)

library(parallel)
library(foreach)
library(doParallel)

num_cores <- detectCores() - 1 # Leave one core free
registerDoParallel(num_cores)

default_colors <- hue_pal()(3)
lighter_colors <- lighten(default_colors, 0.5)
num_reps <- 200
do_read_write <- FALSE

chip_rate <- function(age) {
   # Parameters tuned to match: ~1% at 40, ~10% at 60, ~20% at 75, ~30% at 85
   return(pmax(0, 1e-6 * (age - 20)^3))
}

chip_af <- function(n, age, min_af = 0.02, max_af = 0.4) {
   stopifnot(all(between(age, 0, 100)))

   # TODO: compare this with empiric distributions
   alpha <- 1
   beta <- 20 - 18 * (age / 100)

   af_range <- max_af - min_af

   return(min_af + af_range * rbeta(n, shape1 = alpha, shape2 = beta))
}

wgs_chip_detector <- function(af) {
   stopifnot(all(between(af, 0, 1)))

   # empirical values from WGS vs targeted seq jupyter nb
   detection_prob <- case_when(
      af == 0 ~ 0.06,
      between(af, 0.02, 0.05) ~ 0.09,
      between(af, 0.05, 0.10) ~ 0.32,
      between(af, 0.10, 0.20) ~ 0.67,
      between(af, 0.20, 1.00) ~ 0.86
   )

   return(runif(n = length(af)) < detection_prob)
}

{
   false_pos_wgs_afs <-
      c(
         0.152, 0.103, 0.113, 0.102, 0.111, 0.077, 0.13, 0.162, 0.147,
         0.114, 0.194, 0.1, 0.103, 0.137, 0.125, 0.104, 0.119, 0.075,
         0.091, 0.343, 0.107, 0.06, 0.673, 0.094, 0.184, 0.109, 0.111,
         0.292, 0.347, 0.113, 0.125, 0.105, 0.163, 0.147, 0.088, 0.115,
         0.479, 0.103, 0.22, 0.165, 0.118, 0.069, 0.148, 0.506, 0.179,
         0.088, 0.07, 0.1, 0.133, 0.045, 0.259, 0.086, 0.167, 0.127, 0.217,
         0.08, 0.097, 0.082, 0.138, 0.109, 0.068, 0.093, 0.09, 0.119,
         0.175, 0.115, 0.193, 0.131, 0.12, 0.285, 0.442, 0.125, 0.116,
         0.18, 0.172, 0.125, 0.176, 0.322, 0.1, 0.319, 0.095, 0.241, 0.263,
         0.228, 0.207, 0.218, 0.135, 0.235, 0.13, 0.099, 0.18, 0.072,
         0.106, 0.216, 0.302, 0.21, 0.336, 0.083, 0.118, 0.086, 0.139,
         0.116, 0.075, 0.15, 0.143, 0.113, 0.072, 0.138, 0.112, 0.187,
         0.133, 0.087, 0.237, 0.128, 0.13, 0.104, 0.115, 0.143, 0.096,
         0.097, 0.118, 0.091, 0.229, 0.274, 0.241, 0.206, 0.206, 0.219,
         0.125, 0.111, 0.13, 0.291, 0.103, 0.095, 0.101, 0.395, 0.094,
         0.479, 0.115, 0.199, 0.054, 0.138, 0.111, 0.511, 0.104, 0.061,
         0.113, 0.062, 0.248, 0.08, 0.143, 0.152, 0.414, 0.071, 0.092,
         0.103, 0.108, 0.075, 0.23, 0.153, 0.096, 0.055, 0.26, 0.102,
         0.128, 0.3, 0.093, 0.35, 0.132, 0.081, 0.11, 0.115, 0.077, 0.093,
         0.074, 0.101, 0.07, 0.076, 0.115, 0.204, 0.143, 0.062, 0.258,
         0.137, 0.071, 0.12, 0.187, 0.291, 0.249, 0.133, 0.091, 0.159,
         0.312, 0.552, 0.1, 0.06, 0.088, 0.119, 0.107, 0.176, 0.35, 0.078,
         0.2, 0.111, 0.137, 0.182, 0.159, 0.095, 0.093, 0.2, 0.102, 0.107,
         0.21, 0.279, 0.303, 0.1, 0.467, 0.093, 0.093, 0.194, 0.256, 0.108,
         0.231, 0.311, 0.161, 0.222, 0.481, 0.577, 0.194, 0.077, 0.075,
         0.403, 0.135, 0.322, 0.17, 0.3, 0.243, 0.2, 0.207, 0.294, 0.075,
         0.467, 0.075, 0.162, 0.123, 0.094, 0.175, 0.134, 0.081, 0.198,
         0.096, 0.097, 0.103, 0.083, 0.147, 0.195, 0.206, 0.107, 0.151,
         0.107, 0.093, 0.066, 0.263, 0.111, 0.297, 0.125, 0.103, 0.083,
         0.125, 0.073, 0.107, 0.161, 0.126, 0.137, 0.333, 0.192, 0.079,
         0.185, 0.111, 0.175, 0.071, 0.193, 0.175, 0.156, 0.131, 0.099,
         0.086, 0.095, 0.125, 0.073, 0.272, 0.116, 0.219, 0.088, 0.309,
         0.423, 0.101, 0.119, 0.114, 0.278, 0.088, 0.091, 0.103, 0.107,
         0.267, 0.135, 0.107, 0.084, 0.107, 0.128, 0.166, 0.138, 0.147,
         0.097, 0.093, 0.116, 0.133, 0.159, 0.12, 0.314, 0.128, 0.1, 0.109,
         0.142, 0.096, 0.1, 0.111, 0.2, 0.295, 0.058, 0.1, 0.057, 0.111,
         0.2, 0.167, 0.09, 0.067, 0.12, 0.193, 0.336, 0.121, 0.266, 0.115,
         0.094, 0.373, 0.136, 0.081, 0.141, 0.166, 0.107, 0.057, 0.091,
         0.092, 0.157, 0.214, 0.526, 0.101, 0.103, 0.139, 0.133, 0.124,
         0.129, 0.235, 0.149, 0.131, 0.097, 0.083, 0.165, 0.195, 0.061,
         0.072, 0.074, 0.23, 0.2, 0.081, 0.129, 0.08, 0.511, 0.103, 0.441,
         0.067, 0.158, 0.152, 0.075, 0.172, 0.088, 0.115, 0.075, 0.185,
         0.114, 0.325, 0.157, 0.304, 0.121, 0.066, 0.233, 0.306, 0.12,
         0.146, 0.055, 0.113, 0.134, 0.099, 0.147, 0.107, 0.091, 0.455,
         0.138, 0.146, 0.09, 0.241, 0.359, 0.205, 0.28, 0.091, 0.087,
         0.138, 0.167, 0.088, 0.094, 0.161, 0.22, 0.175, 0.096, 0.055,
         0.35, 0.079, 0.247, 0.137, 0.104, 0.405, 0.125, 0.226, 0.072,
         0.382, 0.1, 0.075, 0.114, 0.35, 0.222, 0.067, 0.332, 0.103, 0.102,
         0.088, 0.218, 0.143, 0.143, 0.231, 0.094, 0.308, 0.097, 0.068,
         0.079, 0.112, 0.315, 0.116, 0.134, 0.089, 0.133, 0.097, 0.369,
         0.12, 0.327, 0.077, 0.179, 0.178, 0.088, 0.122, 0.089, 0.058,
         0.101, 0.159, 0.404, 0.113, 0.135, 0.285, 0.104, 0.448, 0.228,
         0.114, 0.29, 0.086, 0.256, 0.31, 0.106, 0.249, 0.09, 0.105, 0.161,
         0.112, 0.166, 0.229, 0.138, 0.179, 0.348, 0.123, 0.093, 0.114,
         0.088, 0.093, 0.1, 0.089, 0.075, 0.208, 0.135, 0.093, 0.091,
         0.295, 0.103, 0.471, 0.073, 0.201, 0.267, 0.098, 0.081, 0.103,
         0.273, 0.123, 0.072, 0.171, 0.512, 0.309, 0.163, 0.11, 0.178,
         0.122, 0.1, 0.131, 0.094, 0.418, 0.062, 0.138, 0.117, 0.311,
         0.455
      )
}

# Shared components
sim_baseline_data <- function(n) {
   tibble(
      sex = sample(c(0, 1), size = n, replace = TRUE),
      age = sample(40:79, size = n, replace = TRUE),
      p_chip = chip_rate(age),
      chip_status = rbinom(n = n, size = 1, prob = p_chip),
      af = ifelse(chip_status, chip_af(n, age), 0),
      wgs_chip_status = wgs_chip_detector(af),
      false_pos_af = sample(false_pos_wgs_afs, size = n, replace = TRUE),
      wgs_af = case_when(
         chip_status & wgs_chip_status ~ af * 0.43 + 0.1,
         !chip_status & wgs_chip_status ~ false_pos_af,
         .default = 0
      )
   )
}

sim_incident_data <- function(chip_hr, n, beta_age, beta_sex, lambda0) {
   sim_baseline_data(n) %>%
      mutate(
         lp = beta_age * age + beta_sex * sex + log(chip_hr) * chip_status,
         hazard = lambda0 * exp(lp + rnorm(n = n, sd = 1)),
         surv_time = -log(runif(n)) / hazard,
         censoring_time = rexp(n, rate = 0.005 + 0.002 * (age - 40)),
         obs_time = pmin(surv_time, censoring_time),
         event = as.integer(surv_time < censoring_time)
      )
}

sim_prevalent_data <- function(chip_or, n, beta_age, beta_sex) {
   sim_baseline_data(n) %>%
      mutate(
         lp = beta_age * age + beta_sex * sex + log(chip_or) * chip_status - 3,
         disease_prob = plogis(lp),
         disease_status = rbinom(n, 1, disease_prob)
      )
}

sim_chip_study_incident_dz <- function(chip_hr,
                                       n = 1e4,
                                       beta_age = 0.03,
                                       beta_sex = 0.4,
                                       lambda0 = 0.01,
                                       af_thresholds = c(0.02, 0.05, 0.10)) {
   sd <- sim_incident_data(chip_hr, n, beta_age, beta_sex, lambda0)

   tibble(chip_min_af = af_thresholds) |>
      mutate(
         data = map(chip_min_af, ~ filter(sd, wgs_af == 0 | wgs_af > .x)),
         # model = map(data, ~ coxph(Surv(obs_time, event) ~ age + sex + wgs_chip_status, data = .x)),
         model = map(data, ~ coxph(Surv(obs_time, event) ~ age + sex + chip_status, data = .x)),
         results = map(model, ~ tidy(.x) |> slice(3))
      ) |>
      unnest(results) |>
      transmute(chip_min_af, HR = exp(estimate), p = p.value)
}


sim_chip_study_prevalent_dz <- function(chip_or,
                                        n = 1e4,
                                        beta_age = 0.03,
                                        beta_sex = 0.4,
                                        af_thresholds = c(0.02, 0.05, 0.10)) {
   sd <- sim_prevalent_data(chip_or, n, beta_age, beta_sex)
   tibble(chip_min_af = af_thresholds) |>
      mutate(
         data = map(chip_min_af, ~ filter(sd, wgs_af == 0 | wgs_af > .x)),
         # model = map(data, ~ glm(disease_status ~ age + sex + wgs_chip_status,
         model = map(data, ~ glm(disease_status ~ age + sex + chip_status,
            family = binomial, data = .x
         )),
         results = map(model, ~ tidy(.x) |> slice(4)) # 4th row for wgs_chip_status
      ) |>
      unnest(results) |>
      transmute(chip_min_af, OR = exp(estimate), p = p.value)
}

prev_sim_results_df <- foreach(sim_idx = 1:num_reps, .combine = rbind) %dopar% {
   set.seed(sim_idx)
   tibble(chip_or = c(1, 1.5, 2, 2.5, 3)) |>
      mutate(result = map(chip_or, sim_chip_study_prevalent_dz)) |>
      unnest(result) |>
      mutate(sim_idx = sim_idx)
}

incid_sim_results_df <- foreach(sim_idx = 1:num_reps, .combine = rbind) %dopar% {
   set.seed(sim_idx)
   tibble(chip_hr = c(1, 1.5, 2, 2.5, 3)) |>
      mutate(result = map(chip_hr, sim_chip_study_incident_dz)) |>
      unnest(result) |>
      mutate(sim_idx = sim_idx)
}

if (do_read_write) {
   write_rds(prev_sim_results_df, glue("sim_data/prev_sim_results_df_{today()}.rds"))
   write_rds(incid_sim_results_df, glue("sim_data/incid_sim_results_df_{today()}.rds"))
   prev_sim_results_df <- read_rds(glue("sim_data/prev_sim_results_df_{today()}.rds")) |> glimpse()
   incid_sim_results_df <- read_rds(glue("sim_data/incid_sim_results_df_{today()}.rds")) |> glimpse()
}

# prev_sim_results_df = read_rds(glue('sim_data/prev_sim_results_df_2025-07-31.rds')) |> glimpse()
# incid_sim_results_df = read_rds(glue('sim_data/incid_sim_results_df_2025-07-31.rds')) |> glimpse()

prev_sim_results_df |>
   filter(chip_or > 1) |>
   group_by(chip_or, chip_min_af) |>
   summarise(power = mean(p < 0.01)) |>
   print() |>
   ggplot(aes(x = chip_or, y = power)) +
   geom_bar(aes(fill = factor(chip_min_af)),
      stat = "identity",
      position = position_dodge()
   ) +
   scale_y_continuous(
      labels = scales::percent,
      breaks = seq(0, 1, 0.1)
   ) +
   labs(x = "CHIP's odds ratio for disease", fill = "min(AF)", y = "power") +
   scale_fill_manual(values = lighter_colors) +
   theme_bw() +
   theme(legend.position = "bottom") ->
p1c

p1c

incid_sim_results_df |>
   filter(chip_hr > 1) |>
   group_by(chip_hr, chip_min_af) |>
   summarise(power = mean(p < 0.01)) |>
   print() |>
   ggplot(aes(x = chip_hr, y = power)) +
   geom_bar(aes(fill = factor(chip_min_af)),
      stat = "identity",
      position = position_dodge()
   ) +
   scale_y_continuous(
      labels = scales::percent,
      breaks = seq(0, 1, 0.1)
   ) +
   labs(x = "CHIP's hazard ratio for disease", fill = "min(AF)", y = "power") +
   scale_fill_manual(values = lighter_colors) +
   theme_bw() +
   theme(legend.position = "bottom") ->
p1e

p1e

prev_sim_results_df |>
   lm(formula = OR ~ chip_or, data = _) |>
   tidy() |>
   print() ->
prev_coefs

prev_sim_results_df |>
   group_by(chip_or, chip_min_af) |>
   summarise(est_or = mean(OR))

prev_sim_results_df |>
   ggplot(aes(x = chip_or, y = OR)) +
   geom_abline(slope = 1, intercept = 0, linetype = 2) +
   geom_smooth(method = "lm", color = "black") +
   geom_violin(
      aes(
         fill = factor(chip_min_af),
         group = interaction(chip_or, factor(chip_min_af))
      ),
      color = NA, width = 0.4, alpha = 0.5
   ) +
   geom_point(
      aes(
         color = factor(chip_min_af),
         group = interaction(chip_or, factor(chip_min_af))
      ),
      alpha = 0.3,
      size = 0.05,
      position = position_jitterdodge(
         jitter.width = 0.1,
         dodge.width = 0.4
      )
   ) +
   theme_bw() +
   scale_shape_manual(values = c(1, 4)) +
   scale_fill_manual(values = lighter_colors, guide = "none") +
   labs(fill = "min(AF)", x = "odds ratio", y = "estimated odds ratio") +
   guides(color = "none") +
   theme(legend.position = "bottom") ->
p1d

p1d

incid_sim_results_df |>
   group_by(chip_hr) |>
   summarise(mean_HR = mean(HR)) |>
   lm(formula = log(mean_HR) ~ log(chip_hr), data = _) |>
   tidy() |>
   print() ->
incid_coefs

incid_sim_results_df |>
   ggplot(aes(x = chip_hr, y = HR)) +
   geom_abline(slope = 1, intercept = 0, linetype = 2) +
   geom_smooth(method = "lm", color = "black") +
   geom_violin(
      aes(
         fill = factor(chip_min_af),
         group = interaction(chip_hr, factor(chip_min_af))
      ),
      color = NA, width = 0.4, alpha = 0.5
   ) +
   geom_point(
      aes(
         color = factor(chip_min_af),
         group = interaction(chip_hr, factor(chip_min_af))
      ),
      alpha = 0.3,
      size = 0.05,
      position = position_jitterdodge(
         jitter.width = 0.1,
         dodge.width = 0.4
      )
   ) +
   theme_bw() +
   scale_shape_manual(values = c(1, 4)) +
   scale_fill_manual(values = lighter_colors, guide = "none") +
   labs(x = "hazard ratio", y = "estimated hazard ratio") +
   guides(color = "none") +
   theme(legend.position = "bottom") ->
p1f

p1f

if (do_read_write) {
   write_rds(p1c, glue("sim_data/p1c_{today()}.rds"))
   write_rds(p1d, glue("sim_data/p1d_{today()}.rds"))
   write_rds(p1e, glue("sim_data/p1e_{today()}.rds"))
   write_rds(p1f, glue("sim_data/p1f_{today()}.rds"))
}
