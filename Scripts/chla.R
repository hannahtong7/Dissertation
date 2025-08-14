# ─────────────────────────────────────────────────────────────────────────────
# 0) Setup
# ─────────────────────────────────────────────────────────────────────────────
setwd("~/Documents/temp/Dissertation/Data")

library(tidyverse)    # ggplot2, dplyr, readr, tidyr, purrr

# Colour palette (fixed order)
temp_levels <- c("26","30","35")
temp_colors <- c("26"="#4477AA","30"="#66CCAA","35"="#CC6677")

# ─────────────────────────────────────────────────────────────────────────────
# 1) Data import & prep
# ─────────────────────────────────────────────────────────────────────────────
df <- read_csv("chla.csv", show_col_types = FALSE) %>%
  mutate(
    temp     = factor(temp, levels = temp_levels),
    genotype = factor(genotype)
  ) %>%
  # keep only rows with all needed fields
  drop_na(pco2_calculated, Chl_a_per_cell_in_picograms, temp, genotype)

# ─────────────────────────────────────────────────────────────────────────────
# 2) Quadratic summaries (descriptive only)
#    Chl ~ poly(pCO2, 2) per genotype × temp, predict on a fine grid
# ─────────────────────────────────────────────────────────────────────────────
quad_preds <- df %>%
  group_by(genotype, temp) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Chl_a_per_cell_in_picograms ~ poly(pco2_calculated, 2),
                           data = .x)),
    grid  = map(data, ~ tibble(pco2_calculated =
                                 seq(min(.x$pco2_calculated),
                                     max(.x$pco2_calculated), length.out = 200))),
    pred  = map2(model, grid, ~ mutate(.y,
                                       fit = as.numeric(predict(.x, newdata = .y))))
  ) %>%
  select(genotype, temp, pred) %>%
  unnest(pred) %>%
  ungroup()

# ─────────────────────────────────────────────────────────────────────────────
# 3) Observed extrema per panel (for text + appendix table)
# ─────────────────────────────────────────────────────────────────────────────
extrema <- df %>%
  group_by(genotype, temp) %>%
  summarise(
    trough_pCO2 = pco2_calculated[which.min(Chl_a_per_cell_in_picograms)],
    trough_Chl  = min(Chl_a_per_cell_in_picograms),
    crest_pCO2  = pco2_calculated[which.max(Chl_a_per_cell_in_picograms)],
    crest_Chl   = max(Chl_a_per_cell_in_picograms),
    .groups = "drop"
  )

# For plotting triangles
crests <- extrema %>%
  transmute(genotype, temp,
            pco2_calculated = crest_pCO2,
            Chl_a_per_cell_in_picograms = crest_Chl)

troughs <- extrema %>%
  transmute(genotype, temp,
            pco2_calculated = trough_pCO2,
            Chl_a_per_cell_in_picograms = trough_Chl)

# ─────────────────────────────────────────────────────────────────────────────
# 4) Figure: points + thin lines + dashed quadratics + crest/trough markers
# ─────────────────────────────────────────────────────────────────────────────
p_chla <- ggplot(df,
                 aes(x = pco2_calculated,
                     y = Chl_a_per_cell_in_picograms,
                     color = temp, group = temp)) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.6) +
  geom_line(data = quad_preds,
            aes(y = fit),
            linewidth = 0.9, linetype = "22") +  # dashed descriptive fit
  # crest (▲) and trough (▼) markers
  geom_point(data = crests, shape = 24, size = 3, stroke = 0.6) +
  geom_point(data = troughs, shape = 25, size = 3, stroke = 0.6) +
  facet_wrap(~ genotype, scales = "fixed") +
  scale_color_manual(name = "Temperature (°C)",
                     values = temp_colors, drop = FALSE) +
  labs(
    x = expression(paste("pCO"[2], " (µatm, reconstructed)")),
    y = expression(paste(italic("Chl"), italic(" a"),
                         " per cell (pg ", cell^{-1}, ")"))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    strip.text       = element_text(face = "bold", size = 14),
    strip.background = element_blank(),
    legend.position  = "right",
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11),
    axis.title       = element_text(size = 13),
    axis.text        = element_text(size = 12)
  )

print(p_chla)

# Save high-res export for manuscript
ggsave("Figure_Sx_Chla_per_cell.png", p_chla, width = 7.5, height = 3.6, dpi = 600)

# ─────────────────────────────────────────────────────────────────────────────
# 5) Write appendix table (CSV) with crests/troughs
# ─────────────────────────────────────────────────────────────────────────────
write_csv(extrema, "Table_Sx_Chla_extrema.csv")

# (Optional) show table in console
extrema

p_chla_clean <- ggplot(df,
                       aes(x = pco2_calculated,
                           y = Chl_a_per_cell_in_picograms,
                           color = temp)) +
  geom_point(size = 2.6, alpha = 0.9) +                  # points only
  geom_line(data = quad_preds,
            aes(y = fit, group = interaction(genotype, temp)),
            linewidth = 0.9, linetype = "22") +          # dashed descriptive fit
  geom_point(data = crests,  shape = 24, size = 3, stroke = 0.6, fill = NA) +
  geom_point(data = troughs, shape = 25, size = 3, stroke = 0.6, fill = NA) +
  facet_wrap(~ genotype, scales = "fixed") +
  scale_color_manual("Temperature (°C)", values = temp_colors, drop = FALSE) +
  labs(
    x = expression(paste("pCO"[2], " (µatm, reconstructed)")),
    y = expression(paste(italic("Chl"), italic(" a"), " per cell (pg ", cell^{-1}, ")"))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    strip.text       = element_text(face = "bold", size = 14),
    legend.position  = "right"
  )

print(p_chla_clean)




