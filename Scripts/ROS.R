# Set Working Directory
setwd("~/Documents/temp/Dissertation/Data")

# Load Required Libraries
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

# Colorblind-Friendly Palette
temp_colors <- c(
  "26" = "#4477AA",  # Blue
  "30" = "#66CCAA",  # Teal-green
  "35" = "#CC6677"   # Muted red/orange
)

# Step 1: Load and clean ROS data
ros_raw <- read_csv("ROS data.csv") %>%
  rename_with(~ trimws(.)) %>% # Remove any whitespace from column names
  filter(!is.na(temperature)) %>%
  mutate(
    temperature = as.factor(as.integer(temperature)), # Treat as factor for plotting
    genotype = as.factor(genotype),
    co2 = as.numeric(`co2`) # Ensure numeric CO₂
  )

# Step 2: Filter to exclude CCMP2929 at 26ºC (instrument error)
ros_clean <- ros_raw %>%
  filter(!(genotype == "CCMP2929" & temperature == "26"))
# Step 3: Raw Line Plot for Corrected ROS
p_ros <- ggplot(ros_clean, aes(x = co2, y = mean_median_corrected,
                               color = temperature, group = temperature)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ genotype) +
  scale_color_manual(name = "Temperature (°C)", values = temp_colors) +
  labs(
    title = "Corrected ROS Fluorescence vs CO₂",
    x = expression("Reconstructed pCO"[2]~"(ppm)"),
    y = "Mean median DCF fluorescence (corrected, per cell)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x    = element_text(size = 13),
    axis.title.y    = element_text(size = 13),
    strip.text      = element_text(size = 13, face = "bold"),
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11)
  )

# Print the plot
p_ros

