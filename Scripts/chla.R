# Set Working Directory
setwd("~/Documents/temp/Dissertation/Data")

# Load Required Libraries
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# Colorblind-Friendly Color Palette for Temperatures
temp_colors <- c(
  "26" = "#4477AA",  # Blue
  "30" = "#66CCAA",  # Teal-green
  "35" = "#CC6677"   # Muted red/orange
)

# 1. Load data (already contains calculated chlorophyll-a per cell in pg)
df <- read_csv("chla.csv")

# Ensure temp is a factor
df$temp <- as.factor(df$temp)

# 2. Raw Line Plot (Observed Treatments)
p_raw <- ggplot(df, aes(x = pco2_calculated, 
                        y = Chl_a_per_cell_in_picograms,
                        color = temp, group = temp)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ genotype) +
  scale_color_manual(name = "Temperature (°C)", values = temp_colors) +
  labs(
    title = "Raw Chlorophyll-a per Cell (lines connect treatments)",
    x = expression(paste("pCO"[2], " (µatm, calculated)")),
    y = "Chlorophyll-a per cell (pg)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    strip.text   = element_text(face = "bold")
  )

# 3. LOESS Smoothed Plot (Visual trend only)
p_smooth <- ggplot(df, aes(x = pco2_calculated, 
                           y = Chl_a_per_cell_in_picograms,
                           color = temp)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_smooth(method = "loess", span = 1.0, se = FALSE, linewidth = 1) +
  facet_wrap(~ genotype) +
  scale_color_manual(name = "Temperature (°C)", values = temp_colors) +
  labs(
    title = "Smoothed Trends (LOESS fit, no error bars)",
    x = expression(paste("pCO"[2], " (µatm, calculated)")),
    y = "Chlorophyll-a per cell (pg)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    strip.text   = element_text(face = "bold")
  )

# 4. Combine Plots Vertically
combined_plot <- p_raw / p_smooth + plot_layout(heights = c(1, 1.05))
print(combined_plot)




# --- Add this after you load df and set factors ---

library(tidyr)

# Heuristic params per genotype × temp panel
# Model form: y = C + A * exp( - (x - x0) / k )
# where x0 is the panel's xmin, k ~ one-third of the x-range
rng <- df %>%
  group_by(genotype, temp) %>%
  summarise(
    xmin = min(pco2_calculated, na.rm = TRUE),
    xmax = max(pco2_calculated, na.rm = TRUE),
    ymin = min(Chl_a_per_cell_in_picograms, na.rm = TRUE),
    ymax = max(Chl_a_per_cell_in_picograms, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    A  = pmax(ymax - ymin, 0),                 # amplitude
    C  = ymin,                                 # lower asymptote
    x0 = xmin,                                 # start point
    k  = pmax((xmax - xmin) / 3, .Machine$double.eps)  # decay scale
  )

# Build a curve grid per panel
curve_df <- rng %>%
  group_by(genotype, temp) %>%
  do({
    xs <- seq(.$xmin, .$xmax, length.out = 200)
    data.frame(
      genotype = .$genotype,
      temp     = .$temp,
      pco2_calculated = xs,
      y_pred   = .$C + .$A * exp( - (xs - .$x0) / .$k )
    )
  }) %>%
  ungroup()

# Overlay the unfitted exponential-decay curve on your existing plots
p_raw_exp <- p_raw +
  geom_line(
    data = curve_df,
    aes(x = pco2_calculated, y = y_pred, color = temp, group = temp),
    linewidth = 1, linetype = "22"  # dashed
  ) +
  labs(subtitle = "Dashed lines: unfitted exponential-decay reference")

p_smooth_exp <- p_smooth +
  geom_line(
    data = curve_df,
    aes(x = pco2_calculated, y = y_pred, color = temp, group = temp),
    linewidth = 1, linetype = "22"
  ) +
  labs(subtitle = "Dashed lines: unfitted exponential-decay reference")

combined_plot <- p_raw_exp / p_smooth_exp + plot_layout(heights = c(1, 1.05))
print(combined_plot)


