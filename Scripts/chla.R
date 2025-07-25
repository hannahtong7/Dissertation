#Set Working Directory
setwd("~/Documents/temp/Dissertation/Data")

#Load Required Libraries
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

#Colorblind-Friendly Color Palette for Temperatures
temp_colors <- c(
  "26" = "#4477AA",  # Blue
  "30" = "#66CCAA",  # Teal-green
  "35" = "#CC6677"   # Muted red/orange
)

###1. Load and Process Raw Chlorophyll Data
df_raw <- read_csv("Chlorophyll a data.csv")

#chl_a_ug_per_ml is already calculated in the spreadsheet using:
#13.2654 * (A665 - A750) - 2.6839 * (A632 - A750)

#Step 1: Convert cell concentration from µL to mL (provides Chla per cell)
df_raw <- df_raw %>%
  mutate(
    cell_conc_per_mL = cell_conc_per_uL * 1000, #Convert cell concentration to per mL
    cells_filtered = cell_conc_per_mL * chla_extraction_vol_ml, #Calculate total cells filtered from known volumes filtered
    chl_a_total_ug = chla_ug_per_ml * 5,  #Assuming 5 mL methanol - total Chla extracted
    chl_a_per_cell_pg = (chl_a_total_ug / cells_filtered) * 1e6 #Chlorophyll-a per cell in picograms
  )

#Step 2: Save output
write_csv(df_raw, "Chlorophyll_a_per_cell_output.csv")

###2. Load Processed Data for Plotting
df <- read_csv("Chlorophyll_a_per_cell_output.csv")

#Ensure temp is treated as a factor
df$temp <- as.factor(df$temp)

###3. Raw Line Plot (Observed Treatments)
p_raw <- ggplot(df, aes(x = pco2_calculated, y = chl_a_per_cell_pg,
                        color = temp, group = temp)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  facet_wrap(~ genotype) +
  scale_color_manual(name = "Temperature (°C)", values = temp_colors) +
  labs(
    title = "Raw Chlorophyll-a per Cell (lines connect treatments)",
    x = NULL,
    y = "Chlorophyll-a per cell (pg)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title.x    = element_blank(),
    strip.text      = element_text(size = 13, face = "bold"),
    legend.title    = element_text(size = 12),
    legend.text     = element_text(size = 11)
  )

###4. LOESS Smoothed Plot (No Replication, Visual Trend Only)
p_smooth <- ggplot(df, aes(x = pco2_calculated, y = chl_a_per_cell_pg,
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
    plot.title   = element_text(face = "bold", size = 14, hjust = 0.5),
    strip.text   = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 11)
  )

###5. Combine Raw and Smoothed Plots Vertically
combined_plot <- p_raw / p_smooth + plot_layout(heights = c(1, 1.05))
print(combined_plot)


#########Alternative plots#########
df %>%
  group_by(genotype, temp) %>%
  summarise(mean_chla = mean(chl_a_per_cell_pg, na.rm = TRUE)) %>%
  ggplot(aes(x = temp, y = mean_chla, fill = temp)) +
  geom_col() +
  facet_wrap(~ genotype) +
  scale_fill_manual(values = temp_colors) +
  labs(
    title = "Mean Chlorophyll-a per Cell by Temperature",
    x = "Temperature (°C)",
    y = "Mean Chlorophyll-a per cell (pg)"
  ) +
  theme_minimal(base_size = 14)






