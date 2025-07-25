install.packages("brms")
install.packages("cowplot")
library(dplyr)
library(brms) #for bayesian modelling
library(ggplot2)
library(cowplot)

raw_data <- read.csv("~/Documents/temp/Dissertation/Data/Exploratory Data/Temp x CO2_surface_WITH_ADDED_FAKE_DATA_2025_04_29.csv")

#raw_data <-TxCO2_all_data_long_table_SG_2025_05_16

head(raw_data)

#renames variables for clarity/filters out rows of N/A or missing values
dat <- raw_data %>%
  rename(
    CO2 = pco2_calculated,
    temperature = temp
  ) %>%
  select(CO2, temperature, growth_rate, ph_measured, ph_calculated, dic_calculated) %>%
  filter(!is.na(CO2), !is.na(growth_rate), !is.na(temperature))

head(dat)


#dat <- raw_data %>%
#rename(
# CO2 = pco2_calculated,
# temperature = temp
#) %>%
# select(CO2, temperature, growth_rate) %>%
#filter(!is.na(CO2), !is.na(growth_rate), !is.na(temperature))

#Looping over each temperature 
#Filters the dataset for a specific temperature
#The subset is then modeled individually
target_temp <- 26  # <<< choose your temp from 16, 18, 22, 26, 30 or 33
growthdat <- filter(dat, temperature == target_temp)


#Defines the CO2-growth response model that accounts for - this is the nonlinear CO2 model used to fit CO2-growth curves:
#mu_max: maximum growth rate
#alpha: shape/efficiency of CO₂ use
#CO2_opt: CO₂ level where growth is optimal
eilers_peeters_rescaled <- function(CO2, mu_max, alpha, CO2_opt){
  gr <- ((mu_max * CO2) / (((mu_max / ((alpha/10000) * (CO2_opt*100)^2)) * (CO2^2)) +
                             ((1 - ((2 * mu_max)/((alpha/10000) * (CO2_opt*100)))) * CO2) +
                             (mu_max / (alpha/10000))))
  gr
}

#Visualise priors
#Plots the prior distributions used in Bayesian fitting for: mu_max, alpha, CO2_opt
#Helps check if the priors are sensible
plot_prior <- function(x, mean, sd, title) {
  y <- dnorm(x, mean = mean, sd = sd)
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_area(fill = "blue", alpha = 0.5) +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL)
}

p <- plot_prior(seq(0, 3, 0.01), mean = 0.5, sd = 0.5, title = "mu_max")
q <- plot_prior(seq(0, 30/1000, 0.01/1000), mean = 1/1000, sd = 10/1000, title = "alpha")
r <- plot_prior(seq(0, 20*100, 0.01*100), mean = 5*100, sd = 8*100, title = "CO2_opt")

cowplot::plot_grid(p, q, r, nrow = 1)

#Define the priors for bayesian fitting
#This specifies prior beliefs (normal distributions) about the three parameters
#They are all bounded at 0 (lb = 0) since negative values are biologically implausible
prior_eilers_peeters_rescaled <- 
  prior(normal(1.5, 1), nlpar = "mumax", lb = 0) +
  prior(normal(3, 30), nlpar = "alpha", lb = 0) +
  prior(normal(6, 10), nlpar = "CO2opt", lb = 0)

#Define the Bayesian nonlinear formula 
#Defines a nonlinear Bayesian formula for modeling growth_rate as a function of CO₂
#All parameters are treated as constants (intercepts only, no predictors)
eilers_peeters_formula_rescaled <- bf(
  growth_rate ~  ((mumax * CO2) / (((mumax / ((alpha/10000) * (CO2opt*100)^2)) * (CO2^2)) +
                                     ((1 - ((2 * mumax)/((alpha/10000) * (CO2opt*100)))) * CO2) +
                                     (mumax / (alpha/10000)))),
  mumax + alpha + CO2opt ~ 1, nl = TRUE
)

#Fit the model with brms 
compiled_model <- brm(
  formula = eilers_peeters_formula_rescaled,
  data = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 1,
  chains = 1,
  cores = 1,
  sample_prior = "no",
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

fit26 <- update(
  compiled_model,
  newdata = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 5000, #5000 interations
  chains = 4, 
  cores = 4, 
  seed = 11297,
  refresh = 0,
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

print(fit26)
plot(fit26)


#next temp

target_temp <- 16  # <<< choose your temp from 16, 18, 22, 26, 30 or 33
growthdat <- filter(dat, temperature == target_temp)

eilers_peeters_rescaled <- function(CO2, mu_max, alpha, CO2_opt){
  gr <- ((mu_max * CO2) / (((mu_max / ((alpha/10000) * (CO2_opt*100)^2)) * (CO2^2)) +
                             ((1 - ((2 * mu_max)/((alpha/10000) * (CO2_opt*100)))) * CO2) +
                             (mu_max / (alpha/10000))))
  gr
}


plot_prior <- function(x, mean, sd, title) {
  y <- dnorm(x, mean = mean, sd = sd)
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_area(fill = "blue", alpha = 0.5) +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL)
}

p <- plot_prior(seq(0, 3, 0.01), mean = 0.5, sd = 0.5, title = "mu_max")
q <- plot_prior(seq(0, 30/1000, 0.01/1000), mean = 1/1000, sd = 10/1000, title = "alpha")
r <- plot_prior(seq(0, 20*100, 0.01*100), mean = 5*100, sd = 8*100, title = "CO2_opt")

cowplot::plot_grid(p, q, r, nrow = 1)

prior_eilers_peeters_rescaled <- 
  prior(normal(1.5, 1), nlpar = "mumax", lb = 0) +
  prior(normal(3, 30), nlpar = "alpha", lb = 0) +
  prior(normal(6, 10), nlpar = "CO2opt", lb = 0)

eilers_peeters_formula_rescaled <- bf(
  growth_rate ~  ((mumax * CO2) / (((mumax / ((alpha/10000) * (CO2opt*100)^2)) * (CO2^2)) +
                                     ((1 - ((2 * mumax)/((alpha/10000) * (CO2opt*100)))) * CO2) +
                                     (mumax / (alpha/10000)))),
  mumax + alpha + CO2opt ~ 1, nl = TRUE
)

compiled_model <- brm(
  formula = eilers_peeters_formula_rescaled,
  data = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 1,
  chains = 1,
  cores = 1,
  sample_prior = "no",
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

fit16 <- update(
  compiled_model,
  newdata = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 5000,
  chains = 4, 
  cores = 4, 
  seed = 11297,
  refresh = 0,
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

print(fit16)
plot(fit16)


#next temp

target_temp <- 18  # <<< choose your temp from 16, 18, 22, 26, 30 or 33
growthdat <- filter(dat, temperature == target_temp)

eilers_peeters_rescaled <- function(CO2, mu_max, alpha, CO2_opt){
  gr <- ((mu_max * CO2) / (((mu_max / ((alpha/10000) * (CO2_opt*100)^2)) * (CO2^2)) +
                             ((1 - ((2 * mu_max)/((alpha/10000) * (CO2_opt*100)))) * CO2) +
                             (mu_max / (alpha/10000))))
  gr
}


plot_prior <- function(x, mean, sd, title) {
  y <- dnorm(x, mean = mean, sd = sd)
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_area(fill = "blue", alpha = 0.5) +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL)
}

p <- plot_prior(seq(0, 3, 0.01), mean = 0.5, sd = 0.5, title = "mu_max")
q <- plot_prior(seq(0, 30/1000, 0.01/1000), mean = 1/1000, sd = 10/1000, title = "alpha")
r <- plot_prior(seq(0, 20*100, 0.01*100), mean = 5*100, sd = 8*100, title = "CO2_opt")

cowplot::plot_grid(p, q, r, nrow = 1)

prior_eilers_peeters_rescaled <- 
  prior(normal(1.5, 1), nlpar = "mumax", lb = 0) +
  prior(normal(3, 30), nlpar = "alpha", lb = 0) +
  prior(normal(6, 10), nlpar = "CO2opt", lb = 0)

eilers_peeters_formula_rescaled <- bf(
  growth_rate ~  ((mumax * CO2) / (((mumax / ((alpha/10000) * (CO2opt*100)^2)) * (CO2^2)) +
                                     ((1 - ((2 * mumax)/((alpha/10000) * (CO2opt*100)))) * CO2) +
                                     (mumax / (alpha/10000)))),
  mumax + alpha + CO2opt ~ 1, nl = TRUE
)

compiled_model <- brm(
  formula = eilers_peeters_formula_rescaled,
  data = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 1,
  chains = 1,
  cores = 1,
  sample_prior = "no",
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

fit18 <- update(
  compiled_model,
  newdata = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 5000,
  chains = 4, 
  cores = 4, 
  seed = 11297,
  refresh = 0,
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

print(fit18)
plot(fit18)


#next temp

target_temp <- 22  # <<< choose your temp from 16, 18, 22, 26, 30 or 33
growthdat <- filter(dat, temperature == target_temp)

eilers_peeters_rescaled <- function(CO2, mu_max, alpha, CO2_opt){
  gr <- ((mu_max * CO2) / (((mu_max / ((alpha/10000) * (CO2_opt*100)^2)) * (CO2^2)) +
                             ((1 - ((2 * mu_max)/((alpha/10000) * (CO2_opt*100)))) * CO2) +
                             (mu_max / (alpha/10000))))
  gr
}


plot_prior <- function(x, mean, sd, title) {
  y <- dnorm(x, mean = mean, sd = sd)
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_area(fill = "blue", alpha = 0.5) +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL)
}

p <- plot_prior(seq(0, 3, 0.01), mean = 0.5, sd = 0.5, title = "mu_max")
q <- plot_prior(seq(0, 30/1000, 0.01/1000), mean = 1/1000, sd = 10/1000, title = "alpha")
r <- plot_prior(seq(0, 20*100, 0.01*100), mean = 5*100, sd = 8*100, title = "CO2_opt")

cowplot::plot_grid(p, q, r, nrow = 1)

prior_eilers_peeters_rescaled <- 
  prior(normal(1.5, 1), nlpar = "mumax", lb = 0) +
  prior(normal(3, 30), nlpar = "alpha", lb = 0) +
  prior(normal(6, 10), nlpar = "CO2opt", lb = 0)

eilers_peeters_formula_rescaled <- bf(
  growth_rate ~  ((mumax * CO2) / (((mumax / ((alpha/10000) * (CO2opt*100)^2)) * (CO2^2)) +
                                     ((1 - ((2 * mumax)/((alpha/10000) * (CO2opt*100)))) * CO2) +
                                     (mumax / (alpha/10000)))),
  mumax + alpha + CO2opt ~ 1, nl = TRUE
)

compiled_model <- brm(
  formula = eilers_peeters_formula_rescaled,
  data = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 1,
  chains = 1,
  cores = 1,
  sample_prior = "no",
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

fit22 <- update(
  compiled_model,
  newdata = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 5000,
  chains = 4, 
  cores = 4, 
  seed = 11297,
  refresh = 0,
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

print(fit22)
plot(fit22)

#next temp

target_temp <- 26  # <<< choose your temp from 16, 18, 22, 26, 30 or 33
growthdat <- filter(dat, temperature == target_temp)

eilers_peeters_rescaled <- function(CO2, mu_max, alpha, CO2_opt){
  gr <- ((mu_max * CO2) / (((mu_max / ((alpha/10000) * (CO2_opt*100)^2)) * (CO2^2)) +
                             ((1 - ((2 * mu_max)/((alpha/10000) * (CO2_opt*100)))) * CO2) +
                             (mu_max / (alpha/10000))))
  gr
}


plot_prior <- function(x, mean, sd, title) {
  y <- dnorm(x, mean = mean, sd = sd)
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_area(fill = "blue", alpha = 0.5) +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL)
}

p <- plot_prior(seq(0, 3, 0.01), mean = 0.5, sd = 0.5, title = "mu_max")
q <- plot_prior(seq(0, 30/1000, 0.01/1000), mean = 1/1000, sd = 10/1000, title = "alpha")
r <- plot_prior(seq(0, 20*100, 0.01*100), mean = 5*100, sd = 8*100, title = "CO2_opt")

cowplot::plot_grid(p, q, r, nrow = 1)

prior_eilers_peeters_rescaled <- 
  prior(normal(1.5, 1), nlpar = "mumax", lb = 0) +
  prior(normal(3, 30), nlpar = "alpha", lb = 0) +
  prior(normal(6, 10), nlpar = "CO2opt", lb = 0)

eilers_peeters_formula_rescaled <- bf(
  growth_rate ~  ((mumax * CO2) / (((mumax / ((alpha/10000) * (CO2opt*100)^2)) * (CO2^2)) +
                                     ((1 - ((2 * mumax)/((alpha/10000) * (CO2opt*100)))) * CO2) +
                                     (mumax / (alpha/10000)))),
  mumax + alpha + CO2opt ~ 1, nl = TRUE
)

compiled_model <- brm(
  formula = eilers_peeters_formula_rescaled,
  data = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 1,
  chains = 1,
  cores = 1,
  sample_prior = "no",
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

fit26 <- update(
  compiled_model,
  newdata = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 5000,
  chains = 4, 
  cores = 4, 
  seed = 11297,
  refresh = 0,
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

print(fit26)
plot(fit26)

#next temp

target_temp <- 30  # <<< choose your temp from 16, 18, 22, 26, 30 or 33
growthdat <- filter(dat, temperature == target_temp)

eilers_peeters_rescaled <- function(CO2, mu_max, alpha, CO2_opt){
  gr <- ((mu_max * CO2) / (((mu_max / ((alpha/10000) * (CO2_opt*100)^2)) * (CO2^2)) +
                             ((1 - ((2 * mu_max)/((alpha/10000) * (CO2_opt*100)))) * CO2) +
                             (mu_max / (alpha/10000))))
  gr
}


plot_prior <- function(x, mean, sd, title) {
  y <- dnorm(x, mean = mean, sd = sd)
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_area(fill = "blue", alpha = 0.5) +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL)
}

p <- plot_prior(seq(0, 3, 0.01), mean = 0.5, sd = 0.5, title = "mu_max")
q <- plot_prior(seq(0, 30/1000, 0.01/1000), mean = 1/1000, sd = 10/1000, title = "alpha")
r <- plot_prior(seq(0, 20*100, 0.01*100), mean = 5*100, sd = 8*100, title = "CO2_opt")

cowplot::plot_grid(p, q, r, nrow = 1)

prior_eilers_peeters_rescaled <- 
  prior(normal(1.5, 1), nlpar = "mumax", lb = 0) +
  prior(normal(3, 30), nlpar = "alpha", lb = 0) +
  prior(normal(6, 10), nlpar = "CO2opt", lb = 0)

eilers_peeters_formula_rescaled <- bf(
  growth_rate ~  ((mumax * CO2) / (((mumax / ((alpha/10000) * (CO2opt*100)^2)) * (CO2^2)) +
                                     ((1 - ((2 * mumax)/((alpha/10000) * (CO2opt*100)))) * CO2) +
                                     (mumax / (alpha/10000)))),
  mumax + alpha + CO2opt ~ 1, nl = TRUE
)

compiled_model <- brm(
  formula = eilers_peeters_formula_rescaled,
  data = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 1,
  chains = 1,
  cores = 1,
  sample_prior = "no",
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

fit30 <- update(
  compiled_model,
  newdata = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 5000,
  chains = 4, 
  cores = 4, 
  seed = 11297,
  refresh = 0,
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

print(fit30)
plot(fit30)


#next temp

target_temp <- 33  # <<< choose your temp from 16, 18, 22, 26, 30 or 33
growthdat <- filter(dat, temperature == target_temp)

eilers_peeters_rescaled <- function(CO2, mu_max, alpha, CO2_opt){
  gr <- ((mu_max * CO2) / (((mu_max / ((alpha/10000) * (CO2_opt*100)^2)) * (CO2^2)) +
                             ((1 - ((2 * mu_max)/((alpha/10000) * (CO2_opt*100)))) * CO2) +
                             (mu_max / (alpha/10000))))
  gr
}


plot_prior <- function(x, mean, sd, title) {
  y <- dnorm(x, mean = mean, sd = sd)
  ggplot(data.frame(x = x, y = y), aes(x, y)) +
    geom_area(fill = "blue", alpha = 0.5) +
    theme_bw() +
    labs(title = title, x = NULL, y = NULL)
}

p <- plot_prior(seq(0, 3, 0.01), mean = 0.5, sd = 0.5, title = "mu_max")
q <- plot_prior(seq(0, 30/1000, 0.01/1000), mean = 1/1000, sd = 10/1000, title = "alpha")
r <- plot_prior(seq(0, 20*100, 0.01*100), mean = 5*100, sd = 8*100, title = "CO2_opt")

cowplot::plot_grid(p, q, r, nrow = 1)

prior_eilers_peeters_rescaled <- 
  prior(normal(1.5, 1), nlpar = "mumax", lb = 0) +
  prior(normal(3, 30), nlpar = "alpha", lb = 0) +
  prior(normal(6, 10), nlpar = "CO2opt", lb = 0)

eilers_peeters_formula_rescaled <- bf(
  growth_rate ~  ((mumax * CO2) / (((mumax / ((alpha/10000) * (CO2opt*100)^2)) * (CO2^2)) +
                                     ((1 - ((2 * mumax)/((alpha/10000) * (CO2opt*100)))) * CO2) +
                                     (mumax / (alpha/10000)))),
  mumax + alpha + CO2opt ~ 1, nl = TRUE
)

compiled_model <- brm(
  formula = eilers_peeters_formula_rescaled,
  data = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 1,
  chains = 1,
  cores = 1,
  sample_prior = "no",
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

fit33 <- update(
  compiled_model,
  newdata = growthdat,
  prior = prior_eilers_peeters_rescaled, 
  iter = 5000,
  chains = 4, 
  cores = 4, 
  seed = 11297,
  refresh = 0,
  control = list(adapt_delta = 0.95,
                 stepsize = 0.1,
                 max_treedepth = 15)
)

print(fit33)
plot(fit33)


#some pretty violin plots. try first with just one

CO2opt16 <- as_draws_df(fit16, variable = "b_CO2opt_Intercept")$b_CO2opt_Intercept * 100
CO2opt18 <- as_draws_df(fit18, variable = "b_CO2opt_Intercept")$b_CO2opt_Intercept * 100
CO2opt22 <- as_draws_df(fit22, variable = "b_CO2opt_Intercept")$b_CO2opt_Intercept * 100
CO2opt26 <- as_draws_df(fit26, variable = "b_CO2opt_Intercept")$b_CO2opt_Intercept * 100
CO2opt30 <- as_draws_df(fit30, variable = "b_CO2opt_Intercept")$b_CO2opt_Intercept * 100
CO2opt33 <- as_draws_df(fit33, variable = "b_CO2opt_Intercept")$b_CO2opt_Intercept * 100


CO2opt_posterior <- tibble(
  CO2opt       = c(CO2opt16, CO2opt18, CO2opt22, CO2opt26, CO2opt30, CO2opt33),
  Temperature  = rep(unique(dat$temperature),
                     each = length(CO2opt16))
)


ggplot(CO2opt_posterior, aes(x = Temperature, y = CO2opt, group = Temperature)) +
  geom_violin(fill = "lightblue", alpha = 0.6) +
  stat_summary(fun = median, geom = "point", size = 2, color = "red") +
  theme_minimal() +
  labs(x = "Temperature",
       y = "Optimum CO2 for growth (ppm)") 





mumax16 <- as_draws_df(fit16, variable = "b_mumax_Intercept")$b_mumax_Intercept
mumax18 <- as_draws_df(fit18, variable = "b_mumax_Intercept")$b_mumax_Intercept
mumax22 <- as_draws_df(fit22, variable = "b_mumax_Intercept")$b_mumax_Intercept
mumax26 <- as_draws_df(fit26, variable = "b_mumax_Intercept")$b_mumax_Intercept
mumax30 <- as_draws_df(fit30, variable = "b_mumax_Intercept")$b_mumax_Intercept
mumax33 <- as_draws_df(fit33, variable = "b_mumax_Intercept")$b_mumax_Intercept


mumax_posterior <- tibble(
  mumax       = c(mumax16, mumax18, mumax22, mumax26, mumax30, mumax33),
  Temperature  = rep(unique(dat$temperature),
                     each = length(mumax16))
)


ggplot(mumax_posterior, aes(x = Temperature, y = mumax, group = Temperature)) +
  geom_violin(fill = "lightblue", alpha = 0.6) +
  stat_summary(fun = median, geom = "point", size = 2, color = "red") +
  theme_minimal() +
  labs(x = "Temperature",
       y = "Mu max") 



#individual plots if you want them

target_temp <-26
p1_prep <- conditional_effects(fit26, resolution = 100, 
                               spaghetti = TRUE, ndraws = 100)
plot(p1_prep, points = TRUE,
     point_args = c(size = 3, alpha = 0.7), 
     line_args = c(colour = 'black'))[[1]] +
  xlab('CO2 (ppm)') +
  ylab('Growth rate (per day)') +
  ggtitle(paste("Temp =", target_temp, "°C")) +
  theme_bw(base_size = 15) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

target_temp <-30
p1_prep <- conditional_effects(fit26, resolution = 100, 
                               spaghetti = TRUE, ndraws = 100)
plot(p1_prep, points = TRUE,
     point_args = c(size = 3, alpha = 0.7), 
     line_args = c(colour = 'black'))[[1]] +
  xlab('CO2 (ppm)') +
  ylab('Growth rate (per day)') +
  ggtitle(paste("Temp =", target_temp, "°C")) +
  theme_bw(base_size = 15) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


target_temp <-33
p1_prep <- conditional_effects(fit26, resolution = 100, 
                               spaghetti = TRUE, ndraws = 100)
plot(p1_prep, points = TRUE,
     point_args = c(size = 3, alpha = 0.7), 
     line_args = c(colour = 'black'))[[1]] +
  xlab('CO2 (ppm)') +
  ylab('Growth rate (per day)') +
  ggtitle(paste("Temp =", target_temp, "°C")) +
  theme_bw(base_size = 15) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())








