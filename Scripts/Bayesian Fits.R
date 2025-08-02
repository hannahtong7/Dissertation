# ── Load Libraries ─────────────────────────────────────────────────────────────
install.packages("brms")
library(dplyr)
library(brms)
library(ggplot2)
library(cowplot)
library(posterior)
library(tibble)

#1. PARAMETERS ──────────────────────────────────────────────────────────────
gen_to_run <- "CCMP2929"   #e.g. "CCMP1587" or "CCMP2929"

#2. LOAD & FILTER DATA ──────────────────────────────────────────────────────
#This loads the full dataset then filters for only the chosen genotype
raw_data <- read.csv(
  "~/Documents/temp/Dissertation/Data/Exploratory Data/Temp x CO2_surface_WITH_ADDED_FAKE_DATA_AND_HT_DATA.csv"
)

dat <- raw_data %>%
  select(-CO2) %>%  #remove redundant column if exists
  rename(
    CO2         = pco2_calculated,
    temperature = temp
  ) %>%
  filter(
    genotype == gen_to_run,
    !is.na(CO2),
    !is.na(growth_rate),
    !is.na(temperature)
  ) %>%
  select(genotype, temperature, CO2, growth_rate)

temps_to_run <- sort(unique(dat$temperature))
print(temps_to_run)

#3. DEFINE MODEL ────────────────────────────────────────────────────────────
eilers_peeters_rescaled <- function(CO2, mu_max, alpha, CO2_opt){
  (mu_max * CO2) /
    (
      (mu_max / ((alpha/10000) * (CO2_opt*100)^2)) * CO2^2 +
        (1 - (2*mu_max)/((alpha/10000)*(CO2_opt*100))) * CO2 +
        mu_max/(alpha/10000)
    )
}

prior_eilers_peeters <- #my predetermined priors from the literature
  prior(normal(1.5, 1), nlpar = "mumax", lb = 0) +
  prior(normal(3,   30), nlpar = "alpha", lb = 0) +
  prior(normal(6,   10), nlpar = "CO2opt", lb = 0)

#model formula (from the EP model)
model_formula <- bf(
  growth_rate ~ ((mumax * CO2) /
                   (
                     (mumax / ((alpha/10000) * (CO2opt*100)^2)) * CO2^2 +
                       (1 - (2*mumax)/((alpha/10000)*(CO2opt*100))) * CO2 +
                       mumax/(alpha/10000)
                   )
  ),
  mumax + alpha + CO2opt ~ 1,
  nl = TRUE
)

#4. FIT MODEL FOR EACH TEMPERATURE ─────────────────────────────────────────
fits <- list()

for (T in temps_to_run) {
  message("▶ Fitting ", gen_to_run, " at ", T, " °C…")
  
  growthdat <- dat %>% filter(temperature == T)
  
  #compile once
  warmup <- brm(
    formula      = model_formula,
    data         = growthdat,
    prior        = prior_eilers_peeters,
    iter         = 1, chains = 1, cores = 1,
    sample_prior = "no",
    control      = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 15)
  )
  
  #full run - with 4 chains
  fit <- update(
    warmup,
    newdata = growthdat,
    prior   = prior_eilers_peeters,
    iter    = 5000, chains = 4, cores = 4,
    seed    = 11297, refresh = 0,
    control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 15)
  )
  
  fits[[paste0("fit_", T)]] <- fit
}

###DIAGNOSTICS###
# ── DIAGNOSTICS FOR ALL FITS ─────────────────────────────────────────────────────

for (nm in names(fits)) {
  T <- as.numeric(sub("fit_", "", nm))
  fit <- fits[[nm]]
  
  # 1. Summary diagnostics (R̂, CI, ESS)
  cat("\n\n==================== Summary Diagnostics for", gen_to_run, "@", T, "°C ====================\n")
  print(summary(fit))
  
  # 2. Traceplots: check chain mixing and convergence
  trace_plot <- plot(fit, variable = c("b_mumax_Intercept", "b_alpha_Intercept", "b_CO2opt_Intercept"))
  ggsave(filename = paste0("Traceplot_", gen_to_run, "_", T, "C.png"), plot = trace_plot[[1]], width = 6, height = 5)
  
  # 3. Posterior predictive checks: check model fit to observed data
  pp <- pp_check(fit, ndraws = 100) +
    ggtitle(paste("Posterior Predictive Check –", gen_to_run, "@", T, "°C")) +
    theme_minimal(base_size = 14)
  ggsave(filename = paste0("PPCheck_", gen_to_run, "_", T, "C.png"), plot = pp, width = 6, height = 5)
}


# ── 5. PLOT CURVES FOR EACH FITTED TEMPERATURE ───────────────────────────────
for (T in temps_to_run) {
  fit <- fits[[paste0("fit_", T)]]
  
  ce <- conditional_effects(
    fit,
    effects   = "CO2",
    resolution= 100,
    spaghetti = TRUE,
    ndraws    = 100
  )
  
  p <- plot(ce,
            points     = TRUE,
            point_args = list(size = 3, alpha = 0.7),
            line_args  = list(colour = "black")
  )[[1]] +
    labs(
      x     = "CO₂ (ppm)",
      y     = "Growth rate (d⁻¹)",
      title = paste0("Genotype ", gen_to_run, " @ ", T, " °C")
    ) +
    theme_bw(base_size = 15) +
    theme(panel.grid = element_blank())
  
  print(p)
}

# ── 5b. COMBINED CURVE PLOT ACROSS TEMPERATURES ───────────────────────────────

library(tidyverse)  # in case not already loaded

# Common CO₂ range for prediction
CO2_seq <- seq(min(dat$CO2), max(dat$CO2), length.out = 200)

# Use posterior medians to generate fitted curves
curve_data <- map_dfr(temps_to_run, function(T) {
  fit <- fits[[paste0("fit_", T)]]
  
  posterior <- posterior_summary(fit, variable = c("b_mumax_Intercept", "b_alpha_Intercept", "b_CO2opt_Intercept"))
  mu_max  <- posterior["b_mumax_Intercept", "Estimate"]
  alpha   <- posterior["b_alpha_Intercept", "Estimate"]
  CO2_opt <- posterior["b_CO2opt_Intercept", "Estimate"]
  
  tibble(
    temperature = T,
    CO2         = CO2_seq,
    growth      = eilers_peeters_rescaled(CO2_seq, mu_max, alpha, CO2_opt)
  )
})

# Raw data points for visual overlay
dat_points <- dat %>%
  mutate(temperature = factor(temperature))

# Final combined plot
ggplot() +
  geom_line(data = curve_data, aes(x = CO2, y = growth, color = factor(temperature)), linewidth = 1.2) +
  geom_point(data = dat_points, aes(x = CO2, y = growth_rate, color = factor(temperature)), size = 2, alpha = 0.7) +
  scale_color_viridis_d(name = "Temperature (°C)") +
  labs(
    title = paste("Fitted CO₂–growth curves for genotype", gen_to_run),
    x = "CO₂ (ppm)",
    y = "Growth rate (d⁻¹)"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())

# ── 6. VIOLIN PLOTS FOR mumax & CO2opt ────────────────────────────────────────
all_posts <- bind_rows(
  lapply(names(fits), function(nm) {
    T   <- as.numeric(sub("fit_", "", nm))
    fit <- fits[[nm]]
    d   <- as_draws_df(fit, regex = "b_(mumax|CO2opt)_Intercept")
    tibble(
      temperature = T,
      mumax       = d$b_mumax_Intercept,
      CO2opt      = d$b_CO2opt_Intercept * 100
    )
  }),
  .id = "which"
)

lapply(fits, summary)


# μₘₐₓ
ggplot(all_posts, aes(factor(temperature), mumax)) +
  geom_violin(fill = "lightblue", alpha = 0.6) +
  stat_summary(fun = median, geom = "point", size = 2, color = "red") +
  theme_minimal() +
  labs(x = "Temperature (°C)", y = expression(mu[max]))

# CO₂ₒₚₜ
ggplot(all_posts, aes(factor(temperature), CO2opt)) +
  geom_violin(fill = "lightgreen", alpha = 0.6) +
  stat_summary(fun = median, geom = "point", size = 2, color = "red") +
  theme_minimal() +
  labs(x = "Temperature (°C)", y = expression(CO[2]~"opt (ppm)"))

# α extraction (unit: μmol⁻¹)
all_posts <- bind_rows(
  lapply(names(fits), function(nm) {
    T   <- as.numeric(sub("fit_", "", nm))
    fit <- fits[[nm]]
    d   <- as_draws_df(fit, regex = "b_(mumax|CO2opt|alpha)_Intercept")  # ✅ now includes alpha
    tibble(
      temperature = T,
      mumax       = d$b_mumax_Intercept,
      CO2opt      = d$b_CO2opt_Intercept * 100,
      alpha       = d$b_alpha_Intercept
    )
  }),
  .id = "which"
)

ggplot(all_posts, aes(factor(temperature), alpha)) +
  geom_violin(fill = "lightcoral", alpha = 0.5) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  theme_minimal() +
  labs(x = "Temperature (°C)", y = expression(alpha))


#save individual conditional effect plots
for (T in temps_to_run) {
  fit <- fits[[paste0("fit_", T)]]
  ce  <- conditional_effects(fit, effects = "CO2", resolution = 100, spaghetti = TRUE, ndraws = 100)
  plot_ce <- plot(ce, points = TRUE)[[1]] +
    labs(title = paste(gen_to_run, "-", T, "°C"), x = "CO₂ (ppm)", y = "Growth rate (d⁻¹)") +
    theme_minimal()
  
  ggsave(filename = paste0("Fit_", gen_to_run, "_", T, "C.png"), plot = plot_ce, width = 6, height = 5)
}



#Posterior comparison (posterior probabilities) - CO2 opt
all_posts <- bind_rows(
  lapply(names(fits), function(nm) {
    T <- as.numeric(sub("fit_", "", nm))
    fit <- fits[[nm]]
    d <- as_draws_df(fit, regex = "b_CO2opt_Intercept")
    tibble(
      temperature = T,
      CO2opt = d$b_CO2opt_Intercept * 100  # to convert to ppm
    )
  }),
  .id = "which"
)
# Extract posterior draws for each temperature
co2opt_draws <- all_posts %>%
  group_by(temperature) %>%
  summarise(draws = list(CO2opt)) %>%
  arrange(temperature)  # Ensure ordered: 26, 30, 35

# Assign to named variables for clarity
draws_26 <- co2opt_draws$draws[[1]]
draws_30 <- co2opt_draws$draws[[2]]
draws_35 <- co2opt_draws$draws[[3]]

# Calculate pairwise probabilities
prob_30_gt_26 <- mean(draws_30 > draws_26)
prob_35_gt_30 <- mean(draws_35 > draws_30)
prob_35_gt_26 <- mean(draws_35 > draws_26)

# Print results
prob_30_gt_26
prob_35_gt_30
prob_35_gt_26

# Extract μmax and alpha posterior draws
all_posts_other_params <- bind_rows(
  lapply(names(fits), function(nm) {
    T <- as.numeric(sub("fit_", "", nm))
    fit <- fits[[nm]]
    d <- as_draws_df(fit, regex = "b_(mumax|alpha)_Intercept")
    tibble(
      temperature = T,
      mumax = d$b_mumax_Intercept,
      alpha = d$b_alpha_Intercept
    )
  }),
  .id = "which"
)

# μmax posterior comparison
mumax_draws <- all_posts_other_params %>%
  group_by(temperature) %>%
  summarise(draws = list(mumax)) %>%
  arrange(temperature)

draws_mumax_26 <- mumax_draws$draws[[1]]
draws_mumax_30 <- mumax_draws$draws[[2]]
draws_mumax_35 <- mumax_draws$draws[[3]]

prob_mumax_30_gt_26 <- mean(draws_mumax_30 > draws_mumax_26)
prob_mumax_35_gt_30 <- mean(draws_mumax_35 > draws_mumax_30)
prob_mumax_35_gt_26 <- mean(draws_mumax_35 > draws_mumax_26)

# alpha posterior comparison
alpha_draws <- all_posts_other_params %>%
  group_by(temperature) %>%
  summarise(draws = list(alpha)) %>%
  arrange(temperature)

draws_alpha_26 <- alpha_draws$draws[[1]]
draws_alpha_30 <- alpha_draws$draws[[2]]
draws_alpha_35 <- alpha_draws$draws[[3]]

prob_alpha_30_gt_26 <- mean(draws_alpha_30 > draws_alpha_26)
prob_alpha_35_gt_30 <- mean(draws_alpha_35 > draws_alpha_30)
prob_alpha_35_gt_26 <- mean(draws_alpha_35 > draws_alpha_26)

# Print results
list(
  mumax = c(
    "P(30°C > 26°C)" = prob_mumax_30_gt_26,
    "P(35°C > 30°C)" = prob_mumax_35_gt_30,
    "P(35°C > 26°C)" = prob_mumax_35_gt_26
  ),
  alpha = c(
    "P(30°C > 26°C)" = prob_alpha_30_gt_26,
    "P(35°C > 30°C)" = prob_alpha_35_gt_30,
    "P(35°C > 26°C)" = prob_alpha_35_gt_26
  )
)


# Function to summarise posterior difference with CI
summarise_effect_size <- function(high, low) {
  diff <- high - low
  median <- median(diff)
  ci <- quantile(diff, probs = c(0.025, 0.975))
  paste0(round(median, 2), " [", round(ci[1], 2), ", ", round(ci[2], 2), "]")
}

# EFFECT SIZE TABLE: μmax
mumax_effects <- tibble(
  Parameter = "μmax",
  Comparison = c("30°C - 26°C", "35°C - 30°C", "35°C - 26°C"),
  EffectSize = c(
    summarise_effect_size(draws_mumax_30, draws_mumax_26),
    summarise_effect_size(draws_mumax_35, draws_mumax_30),
    summarise_effect_size(draws_mumax_35, draws_mumax_26)
  )
)

# EFFECT SIZE TABLE: CO2opt
co2opt_effects <- tibble(
  Parameter = "CO₂opt (ppm)",
  Comparison = c("30°C - 26°C", "35°C - 30°C", "35°C - 26°C"),
  EffectSize = c(
    summarise_effect_size(draws_30, draws_26),
    summarise_effect_size(draws_35, draws_30),
    summarise_effect_size(draws_35, draws_26)
  )
)

# EFFECT SIZE TABLE: α
alpha_effects <- tibble(
  Parameter = "α",
  Comparison = c("30°C - 26°C", "35°C - 30°C", "35°C - 26°C"),
  EffectSize = c(
    summarise_effect_size(draws_alpha_30, draws_alpha_26),
    summarise_effect_size(draws_alpha_35, draws_alpha_30),
    summarise_effect_size(draws_alpha_35, draws_alpha_26)
  )
)

# Combine into one table
effect_size_table <- bind_rows(mumax_effects, co2opt_effects, alpha_effects)

print(effect_size_table)



#####Final plots######

#1. Bayesian growth rate plots
# Extract conditional effects with uncertainty from each fit
curve_data_ribbon <- map_dfr(temps_to_run, function(T) {
  fit <- fits[[paste0("fit_", T)]]
  
  # Extract conditional effects for CO2
  ce <- conditional_effects(fit, effects = "CO2", resolution = 100)[["CO2"]]
  
  ce %>%
    mutate(
      temperature = factor(T, levels = c(26, 30, 35)),
      genotype = gen_to_run
    ) %>%
    rename(
      CO2 = CO2,
      fit = estimate__,
      lower = lower__,
      upper = upper__
    )
})
# Ensure raw data is properly formatted
dat_points <- dat %>%
  mutate(temperature = factor(temperature, levels = c(26, 30, 35)))
# Assign temperature as factor for legend ordering
curve_data_ribbon <- curve_data_ribbon %>%
  mutate(temperature = factor(temperature, levels = c(26, 30, 35)))

dat_points <- dat_points %>%
  mutate(temperature = factor(temperature, levels = c(26, 30, 35)))

# Custom colors
temp_colors <- c(
  "26" = "#4477AA",
  "30" = "#66CCAA",
  "35" = "#CC6677"
)

ggplot() +
  # 95% credible interval ribbons
  geom_ribbon(
    data = curve_data_ribbon,
    aes(x = CO2, ymin = lower, ymax = upper, fill = temperature, group = temperature),
    alpha = 0.2
  ) +
  
  # Posterior median lines
  geom_line(
    data = curve_data_ribbon,
    aes(x = CO2, y = fit, color = temperature, group = temperature),
    linewidth = 1.1
  ) +
  
  # Raw data points
  geom_point(
    data = dat_points,
    aes(x = CO2, y = growth_rate, color = temperature),
    size = 2, alpha = 0.75
  ) +
  
  # Manual color and fill scales
  scale_color_manual(values = temp_colors, name = "Temperature (°C)") +
  scale_fill_manual(values = temp_colors, name = "Temperature (°C)") +
  
  # Axis and title labels
  labs(
    title = paste("", gen_to_run),
    x = expression("pCO₂ ("*mu*"atm)"),
    y = expression("Specific growth rate ("*mu*", d"^{-1}*")")
  ) +
  
  # Nature-style theme with gridlines and refined font
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "grey85", size = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "grey30"),
    plot.title = element_text(size = 13, hjust = 0.5, face = "plain"),
    legend.position = "right",  # <- ✅ Legend to the side
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

write.csv(curve_data_ribbon, paste0("curve_data_", gen_to_run, ".csv"), row.names = FALSE)
write.csv(dat_points, paste0("dat_points_", gen_to_run, ".csv"), row.names = FALSE)


####Do not run until all genotypes run
# Load both genotypes
curve_1587 <- read.csv("curve_data_CCMP1587.csv")
curve_2929 <- read.csv("curve_data_CCMP2929.csv")
data_1587  <- read.csv("dat_points_CCMP1587.csv")
data_2929  <- read.csv("dat_points_CCMP2929.csv")

# Combine into one dataframe
curve_all <- bind_rows(curve_1587, curve_2929)
data_all  <- bind_rows(data_1587, data_2929)
# Ensure temperature is a factor for plotting
curve_all <- curve_all %>%
  mutate(temperature = factor(temperature, levels = c(26, 30, 35)))

data_all <- data_all %>%
  mutate(temperature = factor(temperature, levels = c(26, 30, 35)))

# Custom colors
temp_colors <- c("26" = "#4477AA", "30" = "#66CCAA", "35" = "#CC6677")

# Final combined plot with facets
ggplot() +
  geom_ribbon(data = curve_all, aes(x = CO2, ymin = lower, ymax = upper,
                                    fill = temperature, group = interaction(temperature, genotype)),
              alpha = 0.2) +
  geom_line(data = curve_all, aes(x = CO2, y = fit,
                                  color = temperature, group = interaction(temperature, genotype)),
            linewidth = 1.1) +
  geom_point(data = data_all, aes(x = CO2, y = growth_rate,
                                  color = temperature),
             size = 2, alpha = 0.75) +
  scale_color_manual(values = temp_colors, name = "Temperature (°C)") +
  scale_fill_manual(values = temp_colors, name = "Temperature (°C)") +
  facet_wrap(~ genotype, ncol = 2, scales = "fixed") +
  labs(
    x = expression("pCO₂ ("*mu*"atm)"),
    y = expression("Specific growth rate ("*mu*", d"^{-1}*")"),
    title = ""
  ) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    panel.grid.major = element_line(color = "grey85", size = 0.3),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.4),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    strip.text = element_text(size = 12, face = "bold")
  )


###For violin plots 

# AFTER FITTING CCMP2929
posteriors_2929 <- bind_rows(
  lapply(names(fits), function(nm) {
    T <- as.numeric(sub("fit_", "", nm))
    d <- as_draws_df(fits[[nm]], regex = "b_(mumax|alpha|CO2opt)_Intercept")
    tibble(
      genotype    = "CCMP2929",
      temperature = factor(T, levels = c(26, 30, 35)),
      mumax       = d$b_mumax_Intercept,
      alpha       = d$b_alpha_Intercept,
      CO2opt      = d$b_CO2opt_Intercept * 100
    )
  })
)
write.csv(posteriors_2929, "posterior_draws_2929.csv", row.names = FALSE)
post_1587 <- read.csv("posterior_draws_1587.csv")
post_2929 <- read.csv("posterior_draws_2929.csv")

posterior_all <- bind_rows(post_1587, post_2929)
library(tidyverse)

