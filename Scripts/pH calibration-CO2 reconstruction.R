#CODE TO RECONSTRUCT PCO2 FROM PH AND TA WITH SEACARB#
#Notes on conventions used here:
#TA units are mol/kg-SW (e.g., 0.0022 ≈ 2200 µmol/kg-SW).
#pH is on the Free scale (pHscale = "F").
#k1k2 = "r" uses Roy et al. (1993) constants for K1 and K2.
#P is total pressure in bar; P = 1 ≈ atmospheric pressure.
#S is practical salinity (unitless).
#Temperatures are in °C.

#Core packages
library(seacarb) #carbonate system calculations
library(writexl) #export results to excel

#############################
#PART 1: PCO2 RECONSTRUCTION
#############################

#Input vectors: pH measured at 20 °C (start of P3) and TA per treatment (n = 6 here)
#Make sure TA is in mol/kg-SW. Values ~0.0022 correspond to 2200 µmol/kg-SW.


pH_values <- c(8.293,
               8.143,
               8.090,
               8.076,
               7.851,
               7.640)

TA_values <- c(0.0022259,
               0.00191,
               0.002071,
               0.0021835,
               0.0021157,
               0.0022676)

#Initialise vectors to store results
DIC_results <- numeric(length(pH_values)) #DIC reconstructed at 20 °C from pH+TA
pH2_results <- numeric(length(pH_values)) #pH recomputed at treatment T from TA+DIC
pCO2_results <- numeric(length(pH_values)) #final pCO2 at treatment T from pH2+TA

#Loop over each set of pH and TA
for (i in seq_along(pH_values)) {
  #Step 1: Calculate DIC from pH and TA (T=20 where pH measured at this temp)
  #Using flag 8 - computes full carbonate system including DIC.
  #set S = 35 (constant), P = 1 bar (surface), pHscale = "F".
  #Pt and Sit are total inorganic P and Si (mol/kg-SW); use small background values.
  result1 <- carb(
    flag = 8,  #input pair: pH + TA
    var1 = pH_values[i], #pH at 20 °C
    var2 = TA_values[i], #TA (mol/kg-SW)
    T = 20, #pH measurement temperature (ºC)
    S = 35, #salinity constant
    P = 1, #pressure (bar)
    Pt = 1e-6, #total phosphate (mol/kg-SW)
    Sit = 10e-6, #total silicate (mol/kg-SW)
    k1k2 = "r", #K1/K2 constants: Roy et al. (1993)
    pHscale = "F" #free pH scale
  )
  DIC1 <- result1$DIC #DIC returned in mol/kg-SW
  DIC_results[i] <- DIC1
  
  #Step 2: Calculate pH at treatment temperature from TA and DIC
  #T=treatment temp
  #Using flag = 15: var1 = TA, var2 = DIC
  #Keep S, P, Pt, Sit, constants, and pH scale consistent.
  
  result2 <- carb(
    flag = 15, #input pair: TA + DIC
    var1 = TA_values[i], #TA (mol/kg-SW)
    var2 = DIC1, #DIC (mol/kg-SW) from Step 1
    T = 26, #treatment temperature
    S = 35, 
    P = 1, 
    Pt = 1e-6, 
    Sit = 10e-6, 
    k1k2 = "r", 
    pHscale = "F"
  )
  pH2 <- result2$pH #pH at treatment temp on the free scale
  pH2_results[i] <- pH2
  
  #Step 3: Calculate final pCO2 from pH2 and TA at treatment temp
  #Using flag 8 again
  result3 <- carb(
    flag = 8, 
    var1 = pH2, #pH at treatment temp from step 2. 
    var2 = TA_values[i], 
    T = 26, #treatment temp to edit
    S = 35, 
    P = 1, 
    Pt = 1e-6, 
    Sit = 10e-6, 
    k1k2 = "r", 
    pHscale = "F"
  )
  pCO2_final <- result3$pCO2 #pCO2 at µatm
  pCO2_results[i] <- pCO2_final
}

#Combine results into a data frame
results_df <- data.frame(
  Initial_pH = pH_values,
  TA = TA_values,
  DIC = DIC_results,
  pH_at_35D = pH2_results,
  Final_pCO2 = pCO2_results
)

#Print the results
print(results_df)

#Export to Excel file (in working directory)
write_xlsx(results_df, path = "carbonate_results_26D_start.xlsx")


##########################################
#PART 2: DATA ANALYSIS AND FIGURE CREATION
##########################################

library(ggplot2)
library(dplyr)

#Nominal (bubbled) CO2 vs reconstructed start pCO2

df <- read.csv("co2_reconstruction_start_and_end.csv") %>%
  mutate(
    temp_f = factor(temp, levels = c(26, 30, 35)),
    gen    = factor(genotype)
  )

ggplot(df, aes(x = co2, y = startpCO2, color = temp_f, shape = genotype)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("26" = "#4477AA", "30" = "#66CCAA", "35" = "#CC6677"),
    name   = "Temperature (°C)"
  ) +
  scale_shape_discrete(name = "Genotype") +
  scale_x_continuous(
    name   = expression("Nominal bubbled CO"[2]~"(ppm)"),
    limits = c(100, 1000),
    breaks = seq(100, 1000, by = 100)
  ) +
  scale_y_continuous(
    name   = expression("Reconstructed pCO"[2]~"(µatm)"),
    limits = c(100, 2500),
    breaks = seq(100, 2500, by = 400)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )


library(ggplot2)
library(dplyr)
# If you want to patch them together, install.packages("patchwork")
library(patchwork)

# 1. Load & prepare
df <- read.csv("co2_reconstruction_start_and_end.csv", check.names=TRUE) %>%
  mutate(
    temp_f   = factor(temp, levels = c(26,30,35)),
    genotype = factor(genotype),
    nominal  = co2,             # ppm
    start    = startpCO2        # µatm
  )

#A. Combined plot
p_combined <- ggplot(df, aes(x = nominal, y = start, color = temp_f, shape = genotype)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  geom_point(size=3) +
  scale_color_manual(values = c("26"="#4477AA","30"="#66CCAA","35"="#CC6677"),
                     name="Temperature (°C)") +
  scale_shape_manual(values = c(16,17), name="Genotype") +
  scale_x_continuous(limits=c(100,1000), breaks=c(100,200,300,400,600,1000),
                     name="Nominal bubbled CO₂ (ppm)") +
  labs(y=expression("Reconstructed start pCO"[2]~"(µatm)"),
       title="Agreement between Nominal and Reconstructed pCO₂") +
  theme_minimal() +
  theme(legend.position="bottom")

#B. Genotype-specific plots
p_1587 <- ggplot(filter(df, genotype=="CCMP1587"),
                 aes(x = nominal, y = start, color = temp_f)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  geom_point(size=3) +
  scale_color_manual(values = c("26"="#4477AA","30"="#66CCAA","35"="#CC6677"),
                     name="Temperature (°C)") +
  scale_x_continuous(limits=c(100,1000), breaks=c(100,200,300,400,600,1000)) +
  labs(
    x="Nominal bubbled CO₂ (ppm)",
    y=expression("Reconstructed start pCO"[2]~"(µatm)"),
    title="CCMP1587: Nominal vs. Reconstructed pCO₂"
  ) +
  theme_minimal() +
  theme(legend.position="none")

p_2929 <- ggplot(filter(df, genotype=="CCMP2929"),
                 aes(x = nominal, y = start, color = temp_f)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey60") +
  geom_point(size=3) +
  scale_color_manual(values = c("26"="#4477AA","30"="#66CCAA","35"="#CC6677"),
                     name="Temperature (°C)") +
  scale_x_continuous(limits=c(100,1000), breaks=c(100,200,300,400,600,1000)) +
  labs(
    x="Nominal bubbled CO₂ (ppm)",
    y=NULL,
    title="CCMP2929: Nominal vs. Reconstructed pCO₂"
  ) +
  theme_minimal() +
  theme(legend.position="none")

#C. Display 
#Combined plus the two separate panels side‑by‑side:
p_combined / (p_1587 | p_2929) + 
  plot_annotation(tag_levels = 'A')  #requires patchwork





#2. Looking at change in pCO2 over the course of the experiment 
setwd("~/Documents/temp/Dissertation/Data")
library(ggplot2)
library(dplyr)

#Custom palette
temp_colors <- c(
  "26" = "#4477AA",  
  "30" = "#66CCAA",  
  "35" = "#CC6677"   
)

#Load data
df <- read.csv("co2_reconstruction_start_and_end.csv")
head(df)
df <- df %>%
  mutate(
    delta_pCO2 = startpCO2 - endpCO2,
    co2_f = factor(co2, levels = c(100,200,300,400,600,1000)),
    temp_f = factor(temp, levels = c(26,30,35))
  )


df2 <- df %>%
  pivot_longer(
    cols = c(startpCO2, endpCO2),
    names_to = "time", values_to = "pCO2"
  ) %>%
  mutate(
    time = factor(time, levels = c("startpCO2","endpCO2"),
                  labels = c("Start","End")),
    co2_f = factor(co2, levels = c(100,200,300,400,600,1000))
  )

ggplot(df2, aes(
  x = pCO2, y = co2_f, group = interaction(co2,temp),
  color = temp_f, shape = time
)) +
  geom_line(size = 0.6) +
  geom_point(size = 3) +
  facet_wrap(~ genotype, nrow = 1) +
  scale_color_manual(values = temp_colors, name = "Temp (°C)") +
  scale_y_discrete(name = "Nominal CO₂ (ppm)") +
  scale_x_continuous(name = expression(pCO[2]~(µatm))) +
  theme_minimal() +
  theme(legend.position="bottom")

library(ggplot2)
library(dplyr)
library(scales)

#Load and prepare in one go
df <- read.csv("co2_reconstruction_start_and_end.csv", stringsAsFactors = FALSE) %>%
  mutate(
    delta_pCO2 = startpCO2 - endpCO2,
    co2_f      = factor(co2, levels = c(100, 200, 300, 400, 600, 1000)),
    temp_f     = factor(temp, levels = c(26, 30, 35))
  )

#Custom diverging palette (blue→white→red)
div_palette <- c("#4477AA", "white", "#CC6677")

#Plot heatmap
ggplot(df, aes(x = co2_f, y = temp_f, fill = delta_pCO2)) +
  geom_tile(color = "grey70") +
  geom_text(aes(label = round(delta_pCO2)), size = 3) +
  facet_wrap(~ genotype, nrow = 1) +
  scale_fill_gradient2(
    low      = div_palette[1],
    mid      = div_palette[2],
    high     = div_palette[3],
    midpoint = 0,
    name     = expression(Delta~pCO[2]~(µatm))
  ) +
  scale_x_discrete(name = "Nominal CO₂ (ppm)") +
  scale_y_discrete(name = "Temperature (°C)") +
  labs(
    title = "Heatmap of CO₂ Drawdown (Start – End)\nby Genotype, Temperature, and CO₂ Level"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(vjust = 0.5),
    legend.position  = "bottom"
  )
