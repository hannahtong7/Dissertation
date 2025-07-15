library(seacarb)

# Input vectors
pH_values <- c(8.387, 8.169, 8.096, 7.986, 7.832, 7.702)
TA_values <- c(2225.9e-6,2405e-6, 2484.1e-6, 2563.2e-6, 2721.4e-6, 3037.8e-6)

# Initialize vectors to store results
DIC_results <- numeric(length(pH_values))
pH2_results <- numeric(length(pH_values))
pCO2_results <- numeric(length(pH_values))

# Loop over each set of pH and TA
for (i in seq_along(pH_values)) {
  # Step 1: Calculate DIC from pH and TA
  result1 <- carb(
    flag = 8, 
    var1 = pH_values[i], 
    var2 = TA_values[i], 
    T = 20, 
    S = 35, 
    P = 1, 
    Pt = 1e-6, 
    Sit = 10e-6, 
    k1k2 = "r", 
    pHscale = "F"
  )
  DIC1 <- result1$DIC
  DIC_results[i] <- DIC1
  
  # Step 2: Calculate pH at 16°C from TA and DIC
  result2 <- carb(
    flag = 15, 
    var1 = TA_values[i], 
    var2 = DIC1, 
    T = 16, 
    S = 35, 
    P = 1, 
    Pt = 1e-6, 
    Sit = 10e-6, 
    k1k2 = "r", 
    pHscale = "F"
  )
  pH2 <- result2$pH
  pH2_results[i] <- pH2
  
  # Step 3: Calculate final pCO2 from pH2 and TA at 16°C
  result3 <- carb(
    flag = 8, 
    var1 = pH2, 
    var2 = TA_values[i], 
    T = 16, 
    S = 35, 
    P = 1, 
    Pt = 1e-6, 
    Sit = 10e-6, 
    k1k2 = "r", 
    pHscale = "F"
  )
  pCO2_final <- result3$pCO2
  pCO2_results[i] <- pCO2_final
}

# Combine results into a data frame
results_df <- data.frame(
  Initial_pH = pH_values,
  TA = TA_values,
  DIC = DIC_results,
  pH_at_33C = pH2_results,
  Final_pCO2 = pCO2_results
)

# Print the results
print(results_df)

