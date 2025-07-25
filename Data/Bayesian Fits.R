# ── Load Libraries ─────────────────────────────────────────────────────────────
install.packages("brms")
library(dplyr)
library(brms)
library(ggplot2)
library(cowplot)
library(posterior)
library(tibble)

# ── 1. PARAMETERS ──────────────────────────────────────────────────────────────
gen_to_run <- "CCMP2929"   # e.g. "CCMP1587" or "CCMP2929"

# ── 2. LOAD & FILTER DATA ──────────────────────────────────────────────────────
raw_data <- read.csv(
  "~/Documents/temp/Dissertation/Data/Exploratory Data/Temp x CO2_surface_WITH_ADDED_FAKE_DATA_AND_HT_DATA.csv"
)

dat <- raw_data %>%
  select(-CO2) %>%  # remove redundant column if exists
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

# ── 2b. Visualise Priors ───────────────────────────────────────────────────────
# Function to plot normal priors (μ, α, CO2opt)
plot_prior <- function(x, mean, sd, title) {
  y <- dnorm(x, mean = mean, sd = sd)
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_area(fill = "blue", alpha = 0.4) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal()
}

# Plot prior distributions (unit conversion included)
p1 <- plot_prior(seq(0, 3, 0.01), mean = 1.5, sd = 1, title = "Prior: μₘₐₓ (day⁻¹)")
p2 <- plot_prior(seq(0, 0.03, 0.0005), mean = 0.003, sd = 0.01, title = "Prior: α (μmol CO₂⁻¹ m² s)")
p3 <- plot_prior(seq(0, 2000, 10), mean = 600, sd = 800, title = "Prior: CO₂ₒₚₜ (ppm)")

cowplot::plot_grid(p1, p2, p3, nrow = 1)



# ── 3. DEFINE MODEL ────────────────────────────────────────────────────────────
eilers_peeters_rescaled <- function(CO2, mu_max, alpha, CO2_opt){
  (mu_max * CO2) /
    (
      (mu_max / ((alpha/10000) * (CO2_opt*100)^2)) * CO2^2 +
        (1 - (2*mu_max)/((alpha/10000)*(CO2_opt*100))) * CO2 +
        mu_max/(alpha/10000)
    )
}

prior_eilers_peeters <- 
  prior(normal(1.5, 1), nlpar = "mumax", lb = 0) +
  prior(normal(3,   30), nlpar = "alpha", lb = 0) +
  prior(normal(6,   10), nlpar = "CO2opt", lb = 0)

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

# ── 4. FIT MODEL FOR EACH TEMPERATURE ─────────────────────────────────────────
fits <- list()

for (T in temps_to_run) {
  message("▶ Fitting ", gen_to_run, " at ", T, " °C…")
  
  growthdat <- dat %>% filter(temperature == T)
  
  # compile once
  warmup <- brm(
    formula      = model_formula,
    data         = growthdat,
    prior        = prior_eilers_peeters,
    iter         = 1, chains = 1, cores = 1,
    sample_prior = "no",
    control      = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 15)
  )
  
  # full run
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



# ── 1. PARAMETERS ──────────────────────────────────────────────────────────────
gen_to_run <- "CCMP1587"   # e.g. "CCMP1587" or "CCMP2929"

# ── 2. LOAD & FILTER DATA ──────────────────────────────────────────────────────
raw_data <- read.csv(
  "~/Documents/temp/Dissertation/Data/Exploratory Data/Temp x CO2_surface_WITH_ADDED_FAKE_DATA_AND_HT_DATA.csv"
)

dat <- raw_data %>%
  select(-CO2) %>%  # remove redundant column if exists
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

# ── 3. DEFINE MODEL ────────────────────────────────────────────────────────────
eilers_peeters_rescaled <- function(CO2, mu_max, alpha, CO2_opt){
  (mu_max * CO2) /
    (
      (mu_max / ((alpha/10000) * (CO2_opt*100)^2)) * CO2^2 +
        (1 - (2*mu_max)/((alpha/10000)*(CO2_opt*100))) * CO2 +
        mu_max/(alpha/10000)
    )
}

prior_eilers_peeters <- 
  prior(normal(1.5, 1), nlpar = "mumax", lb = 0) +
  prior(normal(3,   30), nlpar = "alpha", lb = 0) +
  prior(normal(6,   10), nlpar = "CO2opt", lb = 0)

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

# ── 4. FIT MODEL FOR EACH TEMPERATURE ─────────────────────────────────────────
fits <- list()

for (T in temps_to_run) {
  message("▶ Fitting ", gen_to_run, " at ", T, " °C…")
  
  growthdat <- dat %>% filter(temperature == T)
  
  # compile once
  warmup <- brm(
    formula      = model_formula,
    data         = growthdat,
    prior        = prior_eilers_peeters,
    iter         = 1, chains = 1, cores = 1,
    sample_prior = "no",
    control      = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 15)
  )
  
  # full run
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
ggplot(all_posts, aes(factor(temperature), alpha)) +
  geom_violin(fill = "lightcoral", alpha = 0.5) +
  stat_summary(fun = median, geom = "point", size = 2, color = "black") +
  theme_minimal() +
  labs(x = "Temperature (°C)", y = expression(alpha))

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


# Optional: save individual conditional effect plots
for (T in temps_to_run) {
  fit <- fits[[paste0("fit_", T)]]
  ce  <- conditional_effects(fit, effects = "CO2", resolution = 100, spaghetti = TRUE, ndraws = 100)
  plot_ce <- plot(ce, points = TRUE)[[1]] +
    labs(title = paste(gen_to_run, "-", T, "°C"), x = "CO₂ (ppm)", y = "Growth rate (d⁻¹)") +
    theme_minimal()
  
  ggsave(filename = paste0("Fit_", gen_to_run, "_", T, "C.png"), plot = plot_ce, width = 6, height = 5)
}


#-----------------------------------

library(tidyverse)

# Generate a common CO₂ sequence
CO2_seq <- seq(min(dat$CO2), max(dat$CO2), length.out = 200)

# Define Eilers–Peeters function again
eilers_peeters_rescaled <- function(CO2, mu_max, alpha, CO2_opt) {
  (mu_max * CO2) /
    (
      (mu_max / ((alpha/10000) * (CO2_opt*100)^2)) * CO2^2 +
        (1 - (2*mu_max)/((alpha/10000)*(CO2_opt*100))) * CO2 +
        mu_max/(alpha/10000)
    )
}

# Generate predicted curves for each fitted temperature using median parameter estimates
curve_data <- map_dfr(temps_to_run, function(T) {
  fit <- fits[[paste0("fit_", T)]]
  
  # Extract median parameter estimates
  posterior <- posterior_summary(fit, pars = c("b_mumax_Intercept", "b_alpha_Intercept", "b_CO2opt_Intercept"))
  mu_max  <- posterior["b_mumax_Intercept", "Estimate"]
  alpha   <- posterior["b_alpha_Intercept", "Estimate"]
  CO2_opt <- posterior["b_CO2opt_Intercept", "Estimate"]
  
  # Predict growth across CO2 range
  tibble(
    temperature = T,
    CO2 = CO2_seq,
    growth = eilers_peeters_rescaled(CO2_seq, mu_max, alpha, CO2_opt)
  )
})

# Actual data for scatter
dat_points <- dat %>%
  filter(genotype == gen_to_run) %>%
  mutate(temperature = as.factor(temperature))

# Plot all curves on one plot
ggplot() +
  geom_line(data = curve_data, aes(x = CO2, y = growth, color = factor(temperature)), linewidth = 1.2) +
  geom_point(data = dat_points, aes(x = CO2, y = growth_rate, color = factor(temperature)), size = 2, alpha = 0.7) +
  scale_color_viridis_d(name = "Temperature (°C)") +
  labs(
    title = paste("Fitted CO₂-growth curves for genotype", gen_to_run),
    x = "CO₂ (ppm)",
    y = "Growth rate (d⁻¹)"
  ) +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())






