# ─────────────────────────────────────────────────────────────────────────────
# 0) Setup
# ─────────────────────────────────────────────────────────────────────────────
setwd("~/Documents/temp/Dissertation/Data")

library(tidyverse)
library(minpack.lm)
library(purrr)

temp_levels <- c("26","30","35")
temp_colors <- c("26"="#4477AA","30"="#66CCAA","35"="#CC6677")

# ─────────────────────────────────────────────────────────────────────────────
# 1) Import & clean
# ─────────────────────────────────────────────────────────────────────────────
ros <- read_csv("ROS data.csv", show_col_types = FALSE) %>%
  rename_with(~ trimws(.)) %>%
  filter(!is.na(temperature)) %>%
  mutate(
    temperature = factor(as.integer(temperature), levels = as.integer(temp_levels)) |> fct_drop(),
    genotype    = factor(genotype),
    co2         = as.numeric(co2)
  ) %>%
  drop_na(co2, mean_median_corrected)

# Remove unreliable panel
ros <- ros %>% filter(!(genotype == "CCMP2929" & temperature == "26"))

# ─────────────────────────────────────────────────────────────────────────────
# 2) Log-normal peak fits (descriptive)
#    DCF(x) = a * exp( - (log(x/x0))^2 / (2*b^2) )
# ─────────────────────────────────────────────────────────────────────────────
fit_lognorm <- function(df) {
  tryCatch({
    # sensible starts
    a0  <- max(df$mean_median_corrected, na.rm = TRUE)
    x00 <- df$co2[which.max(df$mean_median_corrected)]
    b0  <- 1
    mod <- nlsLM(mean_median_corrected ~ a * exp(- (log(co2/x0)^2) / (2*b^2)),
                 data = df,
                 start = list(a = a0, x0 = x00, b = b0),
                 control = nls.lm.control(maxiter = 200))
    # predictions on a fine grid
    grid <- tibble(co2 = seq(min(df$co2), max(df$co2), length.out = 200))
    pred <- predict(mod, newdata = grid)
    # simple diagnostics
    rmse <- sqrt(mean((df$mean_median_corrected - predict(mod))^2))
    pars <- coef(mod)
    tibble(a = pars["a"], x0 = pars["x0"], b = pars["b"],
           rmse = rmse,
           pred_ok = TRUE,
           grid_n = nrow(df)) %>%
      mutate(pred_df = list(bind_cols(df[1, c("genotype","temperature")], grid, fit = pred)))
  }, error = function(e) {
    tibble(a = NA_real_, x0 = NA_real_, b = NA_real_,
           rmse = NA_real_, pred_ok = FALSE, grid_n = nrow(df),
           pred_df = list(NULL))
  })
}

fits <- ros %>%
  group_by(genotype, temperature) %>%
  group_split() %>%
  map(~ mutate(.x, .panel_id = paste(first(genotype), first(temperature)))) %>%
  set_names(map_chr(., ~ paste(first(.$genotype), first(.$temperature), sep = "_"))) %>%
  map(~ list(data = .x, fit = fit_lognorm(.x)))

# Extract param table
ros_params <- imap_dfr(fits, function(obj, nm) {
  f <- obj$fit
  rng <- range(obj$data$co2)
  edge_flag <- isTRUE(f$pred_ok) && !is.na(f$x0) &&
    (f$x0 <= (rng[1] + 0.10*diff(rng)) || f$x0 >= (rng[2] - 0.10*diff(rng)))
  tibble(
    genotype    = obj$data$genotype[1],
    temperature = obj$data$temperature[1],
    a           = f$a,
    CO2opt      = f$x0,
    b           = f$b,
    rmse        = f$rmse,
    fit_ok      = f$pred_ok,
    peak_near_edge = edge_flag
  )
})

# Prediction DF (only for panels with usable fits and non-edge peaks)
pred_df <- map_dfr(fits, ~ .x$fit) %>%
  filter(pred_ok) %>%
  pull(pred_df) %>%
  bind_rows() %>%
  inner_join(ros_params %>% filter(fit_ok, !peak_near_edge),
             by = c("genotype","temperature"))

# Optional: peak marker points (x0, a) where fits are usable
peak_pts <- ros_params %>%
  filter(fit_ok, !is.na(CO2opt), !peak_near_edge) %>%
  select(genotype, temperature, co2 = CO2opt, mean_median_corrected = a)

# ─────────────────────────────────────────────────────────────────────────────
# 3) Plot (publication)
# ─────────────────────────────────────────────────────────────────────────────
p_ros <- ggplot(ros, aes(co2, mean_median_corrected, color = temperature)) +
  geom_point(size = 2.4, alpha = 0.9) +
  geom_line(
    data = pred_df,
    aes(x = co2, y = fit,
        color = temperature,
        group = interaction(genotype, temperature)),
    linewidth = 0.9, linetype = "22", inherit.aes = FALSE
  ) +
  geom_point(
    data = peak_pts, shape = 24, size = 3, stroke = 0.6, fill = NA
  ) +
  facet_wrap(~ genotype) +
  scale_color_manual("Temperature (°C)", values = temp_colors, drop = FALSE) +
  labs(
    x = expression("pCO"[2]*" (µatm)"),
    y = expression("DCF fluorescence (a.u. cell"^-1*")")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(size = 13, face = "bold"),
    legend.position  = "right"
  )


print(p_ros)
ggsave("Figure_Sy_ROS.png", p_ros, width = 7.5, height = 3.6, dpi = 600)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Write appendix table
# ─────────────────────────────────────────────────────────────────────────────
write_csv(ros_params, "Table_Sy_ROS_fit_parameters.csv")
ros_params
