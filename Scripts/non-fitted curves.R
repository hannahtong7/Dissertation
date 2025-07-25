#1. Setup ───────────────────────────────────────────────────────────────────

#Set working directory
setwd("~/Documents/temp/Dissertation/Data/Exploratory Data")

#Load required libraries
library(ggplot2)
library(dplyr)
library(readr)
library(scales)

#2. Load & Clean Data ───────────────────────────────────────────────────────

#Read CSV file
df <- read_csv("Temp x CO2_surface_WITH_ADDED_FAKE_DATA_AND_HT_DATA_NO_ZERO.csv")

#Filter to relevant genotypes and temperatures
target_gens  <- c("CCMP1587", "CCMP2929")
target_temps <- c(26, 30, 35)

df_clean <- df %>%
  filter(genotype %in% target_gens,
         temp %in% target_temps) %>%
  mutate(
    pco2_calculated = as.numeric(pco2_calculated),
    growth_rate     = as.numeric(growth_rate)
  ) %>%
  filter(is.finite(pco2_calculated), is.finite(growth_rate))

#3. Define Colorblind-Friendly Temperature Palette ──────────────────────────

temp_colors <- c(
  "26" = "#4477AA",  # Blue
  "30" = "#66CCAA",  # Teal-green
  "35" = "#CC6677"   # Muted red-orange
)

#4. Plot 1: One Plot per Genotype (CO₂ vs Growth at Each Temp) ──────────────

make_genotype_plot <- function(gen) {
  ggplot(filter(df_clean, genotype == gen),
         aes(x = pco2_calculated,
             y = growth_rate,
             color = factor(temp),
             group = temp)) +
    geom_line(size = 0.8) +
    geom_point(size = 2) +
    scale_color_manual(values = temp_colors, name = "Temperature (°C)") +
    labs(
      title = paste("CO2–Growth Response,", gen),
      x = expression(paste("pCO2", " (µatm)")),
      y = expression(paste("Growth rate (d"^{-1}, ")"))
    ) +
    theme_minimal(base_size = 14)
}

p1 <- make_genotype_plot("CCMP1587")
p2 <- make_genotype_plot("CCMP2929")

print(p1)
print(p2)

#5. Plot 2: Faceted by Temperature, Colored by Genotype ─────────────────────

p_by_temp <- ggplot(df_clean,
                    aes(x = pco2_calculated,
                        y = growth_rate,
                        color = genotype,
                        group = genotype)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ temp, ncol = 2,
             labeller = labeller(temp = function(t) paste0(t, " °C"))) +
  scale_color_brewer(palette = "Dark2", name = "Genotype") +
  labs(
    title = "CO2–Growth Response by Genotype at Each Temperature",
    x = expression(paste("pCO2", " (µatm)")),
    y = expression(paste("Growth rate (d"^{-1}, ")"))
  ) +
  theme_minimal(base_size = 14)

print(p_by_temp)

#6. Plot 3: Genotype × Temperature Grid with Consistent Axes ────────────────

p_facet_grid <- ggplot(df_clean,
                       aes(x = pco2_calculated,
                           y = growth_rate,
                           color = factor(temp),
                           group = temp)) +
  geom_line(size = 0.8) +
  geom_point(size = 2, shape = 21, fill = "white", stroke = 1) +
  facet_grid(rows = vars(genotype), cols = vars(temp),
             labeller = labeller(
               temp     = function(t) paste0(t, " °C"),
               genotype = label_value
             )) +
  scale_color_manual(name = "Temperature (°C)", values = temp_colors) +
  scale_y_continuous(
    limits = c(1.25, 2.25),
    breaks = seq(1.25, 2.25, by = 0.25),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_x_continuous(
    breaks = pretty_breaks(n = 5),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "Growth Rate Responses to CO2 Across Temperatures and Genotypes",
    x = expression(paste("pCO2", " (µatm)")),
    y = expression(paste("Growth rate (d"^{-1}, ")"))
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.background    = element_blank(),
    strip.text          = element_text(size = 13, face = "bold"),
    legend.position     = "right",
    legend.title        = element_text(size = 12),
    legend.text         = element_text(size = 11),
    axis.text           = element_text(size = 11),
    axis.title          = element_text(size = 13),
    panel.grid.major.y  = element_line(color = "grey85", size = 0.4),
    panel.grid.major.x  = element_line(color = "grey90", size = 0.3),
    panel.grid.minor.x  = element_blank(),
    panel.grid.minor.y  = element_blank(),
    plot.title          = element_text(size = 16, face = "bold", hjust = 0.5)
  )

print(p_facet_grid)
