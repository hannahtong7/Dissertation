# Load libraries
library(rTPC)
library(nls.multstart)
library(tidyverse)
library(broom)
setwd("~/Documents/temp/Dissertation/Data")
# Load data and fit 7D
df_7D <- read.csv("TPC_data.csv") %>%
  filter(transfer == "tr4", strain == "7D") %>%
  rename(rate = gr, temp_C = temp)

start_vals_7D <- get_start_vals(df_7D$temp_C, df_7D$rate, model_name = "sharpeschoolhigh_1981")
low_lims_7D   <- get_lower_lims(df_7D$temp_C, df_7D$rate, model_name = "sharpeschoolhigh_1981")
upper_lims_7D <- get_upper_lims(df_7D$temp_C, df_7D$rate, model_name = "sharpeschoolhigh_1981")

fit_7D <- nls_multstart(
  rate ~ sharpeschoolhigh_1981(temp = temp_C, r_tref, e, eh, th, tref = 15),
  data = df_7D,
  iter = 500,
  start_lower = start_vals_7D - 10,
  start_upper = start_vals_7D + 10,
  lower = low_lims_7D,
  upper = upper_lims_7D,
  supp_errors = "Y"
)

new_data_7D <- data.frame(temp_C = seq(min(df_7D$temp_C), max(df_7D$temp_C), 0.5))
preds_7D <- augment(fit_7D, newdata = new_data_7D) %>% mutate(strain = "7D")
df_7D <- df_7D %>% mutate(strain = "7D")

# Load data and fit 7E
df_7E <- read.csv("TPC_data.csv") %>%
  filter(transfer == "tr4", strain == "7E") %>%
  rename(rate = gr, temp_C = temp)

start_vals_7E <- get_start_vals(df_7E$temp_C, df_7E$rate, model_name = "sharpeschoolhigh_1981")
low_lims_7E   <- get_lower_lims(df_7E$temp_C, df_7E$rate, model_name = "sharpeschoolhigh_1981")
upper_lims_7E <- get_upper_lims(df_7E$temp_C, df_7E$rate, model_name = "sharpeschoolhigh_1981")

fit_7E <- nls_multstart(
  rate ~ sharpeschoolhigh_1981(temp = temp_C, r_tref, e, eh, th, tref = 15),
  data = df_7E,
  iter = 500,
  start_lower = start_vals_7E - 10,
  start_upper = start_vals_7E + 10,
  lower = low_lims_7E,
  upper = upper_lims_7E,
  supp_errors = "Y"
)

new_data_7E <- data.frame(temp_C = seq(min(df_7E$temp_C), max(df_7E$temp_C), 0.5))
preds_7E <- augment(fit_7E, newdata = new_data_7E) %>% mutate(strain = "7E")
df_7E <- df_7E %>% mutate(strain = "7E")

# Combine datasets
all_data <- bind_rows(df_7D, df_7E)
all_preds <- bind_rows(preds_7D, preds_7E)

# Plot both curves
ggplot(all_data, aes(x = temp_C, y = rate, color = strain)) +
  geom_point(size = 2) +
  geom_line(data = all_preds, aes(x = temp_C, y = .fitted, color = strain), linewidth = 1.2) +
  scale_color_manual(values = c("7D" = "blue", "7E" = "darkorange")) +
  theme_bw(base_size = 14) +
  labs(
    x = "Temperature (ºC)",
    y = "Growth rate (per day)",
    title = "Sharpe–Schoolfield TPCs for Strains 7D and 7E (Transfer 4)",
    color = "Strain"
  )

