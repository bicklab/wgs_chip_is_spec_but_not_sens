library(glue)
library(tidyverse)
library(broom)
library(survival)
library(scales)
library(colorspace)

library(parallel)
library(foreach)
library(doParallel)

num_cores <- detectCores() - 1  # Leave one core free
registerDoParallel(num_cores)

chip_rate = function(age) {
   # Parameters tuned to match: ~1% at 40, ~10% at 60, ~20% at 75, ~30% at 85
   return(pmax(0, 1e-6 * (age - 20)^3))
}

chip_af = function(n, age, min_af = 0.02, max_af = 0.4) {

   stopifnot(all(between(age, 0, 100)))

   # TODO: compare this with empiric distributions
   alpha <- 1
   beta <- 20 - 18*(age/100)

   af_range = max_af - min_af

   return(min_af + af_range*rbeta(n, shape1 = alpha, shape2 = beta))
}

wgs_chip_detector = function(af) {

   stopifnot(all(between(af, 0, 1)))

   # empirical values from WGS vs targeted seq jupyter nb
   detection_prob <- case_when(
      af == 0   ~ 0.03,
      between(af, 0.02, 0.05) ~ 0.13,
      between(af, 0.05, 0.10) ~ 0.33,
      between(af, 0.10, 0.20) ~ 0.63,
      between(af, 0.20, 0.30) ~ 0.48,
      TRUE      ~ 0.15
   )

   return(runif(n = length(af)) < detection_prob)
}

{
   false_pos_wgs_afs = c(
      0.102, 0.187, 0.229, 0.179, 0.075, 0.075, 0.135, 0.125, 0.096,
      0.22, 0.165, 0.103, 0.086, 0.079, 0.105, 0.075, 0.18, 0.061,
      0.147, 0.08, 0.205, 0.143, 0.143, 0.138, 0.125, 0.162, 0.112,
      0.1, 0.128, 0.088, 0.071, 0.17, 0.115, 0.172, 0.243, 0.199, 0.077,
      0.089, 0.167, 0.148, 0.194, 0.195, 0.159, 0.072, 0.138, 0.104,
      0.073, 0.143, 0.068, 0.204, 0.114, 0.107, 0.171, 0.131, 0.119,
      0.193, 0.134, 0.086, 0.079, 0.137, 0.071, 0.086, 0.097, 0.163,
      0.152, 0.088, 0.119, 0.122, 0.1, 0.235, 0.137, 0.116, 0.241,
      0.142, 0.103, 0.06, 0.108, 0.136, 0.116, 0.12, 0.071, 0.088,
      0.193, 0.097, 0.099, 0.159, 0.201, 0.147, 0.195, 0.108, 0.111,
      0.206, 0.206, 0.219, 0.161, 0.057, 0.157, 0.074, 0.135, 0.091,
      0.113, 0.217, 0.11, 0.091, 0.095, 0.091, 0.134, 0.075, 0.138,
      0.179, 0.1, 0.233, 0.102, 0.13, 0.133, 0.159, 0.055, 0.131, 0.125,
      0.133, 0.103, 0.179, 0.103, 0.107, 0.101, 0.2, 0.128, 0.093,
      0.1, 0.176, 0.249, 0.09, 0.128, 0.133, 0.058, 0.093, 0.114, 0.222,
      0.093, 0.231, 0.125, 0.12, 0.222, 0.121, 0.266, 0.081, 0.182,
      0.2, 0.111, 0.075, 0.088, 0.143, 0.105, 0.161, 0.113, 0.104,
      0.147, 0.088, 0.226, 0.093, 0.074, 0.139, 0.092, 0.1, 0.066,
      0.2, 0.138, 0.2, 0.091, 0.2, 0.121, 0.07, 0.185, 0.079, 0.072,
      0.089, 0.055, 0.115, 0.241, 0.131, 0.103, 0.137, 0.111, 0.077,
      0.102, 0.104, 0.172, 0.13, 0.062, 0.115, 0.107, 0.167, 0.119,
      0.176, 0.062, 0.137, 0.058, 0.125, 0.123, 0.093, 0.116, 0.061,
      0.103, 0.1, 0.141, 0.08, 0.089, 0.147, 0.067, 0.11, 0.083, 0.165,
      0.115, 0.111, 0.22, 0.123, 0.21, 0.1, 0.13, 0.178, 0.101, 0.218,
      0.091, 0.103, 0.114, 0.055, 0.076, 0.095, 0.113, 0.062, 0.075,
      0.095, 0.147, 0.083, 0.101, 0.161, 0.091, 0.1, 0.054, 0.146,
      0.081, 0.135, 0.1, 0.094, 0.138, 0.214, 0.115, 0.097, 0.138,
      0.106, 0.118, 0.069, 0.124, 0.09, 0.084, 0.075, 0.125, 0.127,
      0.073, 0.139, 0.094, 0.175, 0.23, 0.267, 0.123, 0.045, 0.114,
      0.107, 0.175, 0.175, 0.086, 0.1, 0.166, 0.118, 0.093, 0.077,
      0.111, 0.107, 0.241, 0.193, 0.111, 0.113, 0.125, 0.137, 0.096,
      0.088, 0.103, 0.102, 0.162, 0.098, 0.107, 0.175, 0.129, 0.258,
      0.104, 0.12, 0.099, 0.073, 0.133, 0.166, 0.094, 0.094, 0.206,
      0.107, 0.12, 0.159, 0.2, 0.125, 0.091, 0.151, 0.194, 0.088, 0.235,
      0.111, 0.087, 0.133, 0.138
   )
}

sim_data = function(chip_hr, n, beta_age, beta_sex, lambda0) {

   tibble(
      # reality
      sex = sample(c(0, 1), size = n, replace = TRUE),
      age = sample(40:79, size = n, replace = TRUE),
      p_chip = chip_rate(age),
      chip_status = rbinom(n = n, size = 1, prob = p_chip),
      af = ifelse(chip_status, chip_af(n, age), 0),
      # TODO consider higher LP for bigger CHIP
      lp = beta_age*age + beta_sex*sex + log(chip_hr)*chip_status,
      hazard = lambda0*exp(lp + rnorm(n = n, sd = 1)),
      surv_time = -log(runif(n))/hazard,
      censoring_time = rexp(n, rate = 0.005 + 0.002*(age - 40)),
      obs_time = pmin(surv_time, censoring_time),
      event = as.integer(surv_time < censoring_time),
      # experiment
      wgs_chip_status = wgs_chip_detector(af),
      false_pos_af = sample(false_pos_wgs_afs, size = n, replace = TRUE),
      wgs_af = case_when(
         chip_status & wgs_chip_status ~ af*0.43 + 0.1,
         !chip_status & wgs_chip_status ~ false_pos_af,
         .default = 0
      )
   )
}


sim_chip_epi_study = function(chip_hr,
                              chip_min_af,
                              n = 1e5,
                              beta_age = 0.03,
                              beta_sex = 0.4,
                              lambda0 = 0.01) {

   # message('Simulating data...', appendLF = FALSE)
   data_sim = sim_data(chip_hr, n, beta_age, beta_sex, lambda0)
   # message('done')

   # message('Fitting DTS model...', appendLF = FALSE)
   data_sim |>
      filter(af == 0 | af > chip_min_af) |>
      coxph(Surv(obs_time, event) ~ age + sex + chip_status, data = _) |>
      tidy() ->
      dts_fit
   # message('done')

   # message('Fitting WGS model...', appendLF = FALSE)
   data_sim |>
      filter(wgs_af == 0 | wgs_af > chip_min_af) |>
      coxph(Surv(obs_time, event) ~ age + sex + wgs_chip_status, data = _) |>
      tidy() ->
      wgs_fit
   # message('done')

   return(
      tibble(
         `gold standard` = exp(dts_fit$estimate[3]),
         `gold standard p` = dts_fit$p.value[3],
         `WGS-based` = exp(wgs_fit$estimate[3]),
         `WGS-based p` = wgs_fit$p.value[3]
      )
   )
}

sim_results_list = list()
num_reps = 100

expand_grid(
   chip_hr = c(1, 1.1, 1.2, 1.25, 1.3, 1.4, 1.5, 1.75, 2, 2.5, 3),
   chip_min_af = c(0.02, 0.05, 0.1)) |>
   glimpse() ->
   param_grid

# sims that run in sequence
{
   # for (sim_idx in 1:num_reps) {
   #
   #    message(sim_idx)
   #
   #    param_grid |>
   #       mutate(result = map2_df(chip_hr, chip_min_af, sim_chip_epi_study)) |>
   #       unnest_wider(result) ->
   #       sim_results
   #
   #    sim_results_list[[sim_idx]] = tibble(sim_idx = sim_idx,
   #                                         sim_results)
   #
   # }
   #
   # sim_results_list |>
   #    bind_rows() |>
   #    glimpse() ->
   #    sim_results_df
}

# sims that run in parallel
{
   sim_results_df <- foreach(sim_idx = 1:num_reps, .combine = rbind) %dopar% {
      param_grid |>
         mutate(result = map2_df(chip_hr, chip_min_af, sim_chip_epi_study)) |>
         unnest_wider(result) |>
         mutate(sim_idx = sim_idx)
   }
}

write_rds(sim_results_df, glue('sim_results_df_{today()}.rds'))
sim_results_df = read_rds(glue('sim_results_df_{today()}.rds')) |> glimpse()

sim_results_df |>
   group_by(chip_hr, chip_min_af) |>
   summarise(dts_power = mean(`gold standard p` < 0.01),
             wgs_power = mean(`WGS-based p` < 0.01)) |>
   print(n = Inf)



default_colors <- hue_pal()(3)
lighter_colors <- lighten(default_colors, 0.3)

sim_results_df |>
   filter(chip_hr %in% c(1.1, 1.2, 1.3, 1.4)) |>
   group_by(chip_hr, chip_min_af) |>
   summarise(DTS = mean(`gold standard p` < 0.01),
             WGS = mean(`WGS-based p` < 0.01)) |>
   pivot_longer(cols = c('DTS', 'WGS')) |>
   ggplot(aes(x = name, y = value)) +
   geom_bar(aes(fill = factor(chip_min_af)),
            stat = 'identity',
            position = position_dodge()) +
   geom_text(aes(x = name,
                  y = value - 5/100,
                  label = scales::percent(value),
                  group = factor(chip_min_af)),
              position = position_dodge(width = 0.9),
              size = 1.7) +
   facet_grid(cols = vars(chip_hr), labeller = as_labeller(function(x) glue('HR = {x}'))) +
   scale_y_continuous(labels = scales::percent) +
   labs(x = NULL, fill = 'min(AF)', y = 'power') +
   scale_fill_manual(values = lighter_colors) +
   theme_bw() +
   theme(legend.position = 'bottom')

ggsave(filename = glue('power_from_sims{today()}.pdf'),
       height = 3, width = 8)

is_int = function(x) { (x %% 1) == 0 }

bind_cols(
   pivot_longer(data = filter(sim_results_df, is_int(chip_hr * 4)),
                cols = c(`gold standard`, `WGS-based`)) |>
      select(sim_idx, chip_hr, chip_min_af, chip_call_type = name, beta = value),
   pivot_longer(data = filter(sim_results_df, is_int(chip_hr * 4)),
                cols = c(`gold standard p`, `WGS-based p`)) |>
      select(chip_call_type_p = name, p = value)) |>
   glimpse() ->
   toplot

toplot |>
   filter(chip_hr %in% c(1, 1.5, 2, 2.5, 3)) |>
   group_by(chip_call_type) |>
   summarise(
      model = list(lm(beta ~ chip_hr, data = pick(everything()))),
      .groups = 'drop'
   ) |>
   mutate(tidy = map(model, broom::tidy)) |>
   unnest(tidy) |>
   # filter(term == "chip_hr") |>
   print() ->
   chip_coefs

toplot |>
   filter(chip_hr %in% c(1, 1.5, 2, 2.5, 3)) |>
   ggplot(aes(x = chip_hr, y = beta)) +
   geom_abline(slope = 1, intercept = 0, linetype = 2) +
   geom_smooth(method = 'lm', color = 'black') +
   geom_violin(aes(fill = factor(chip_min_af),
                   group = interaction(chip_hr, factor(chip_min_af))),
               color = NA, width = 0.4, alpha = 0.5) +
   geom_point(aes(color = factor(chip_min_af),
                  group = interaction(chip_hr, factor(chip_min_af))),
              alpha = 0.3,
              size = 0.05,
              position = position_jitterdodge(jitter.width = 0.1,
                                              dodge.width = 0.4)) +
   facet_grid(cols = vars(chip_call_type)) +
   theme_bw() +
   scale_shape_manual(values = c(1, 4)) +
   scale_fill_manual(values = lighter_colors) +
   labs(fill = 'min(AF)', x = 'hazard ratio', y = 'estimated hazard ratio') +
   guides(color = 'none')

ggsave(filename = glue('betas_from_sims_violin_and_dots_{today()}.pdf'),
       height = 4, width = 8)







