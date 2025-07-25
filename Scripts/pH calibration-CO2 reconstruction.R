library(seacarb)
library(writexl)

#Input vectors - pH (start of P3) and TA values per treatment n=36

pH_values <- c(7.713,
               7.723,
               7.825,
               7.732,
               7.716,
               7.573)
TA_values <- c(1687.4e-6, 
               1830.2e-6, 
               1995.7e-6, 
               2018.6e-6, 
               2130.5e-6, 
               2015.1e-6)

#Initialize vectors to store results
DIC_results <- numeric(length(pH_values))
pH2_results <- numeric(length(pH_values))
pCO2_results <- numeric(length(pH_values))

#Loop over each set of pH and TA
for (i in seq_along(pH_values)) {
  #Step 1: Calculate DIC from pH and TA (T=20 where pH measured at this temp)
  result1 <- carb(
    flag = 8, 
    var1 = pH_values[i], 
    var2 = TA_values[i], 
    T = 20, 
    S = 35, #salinity constant
    P = 1, 
    Pt = 1e-6, 
    Sit = 10e-6, 
    k1k2 = "r", 
    pHscale = "F"
  )
  DIC1 <- result1$DIC
  DIC_results[i] <- DIC1
  
  #Step 2: Calculate pH at treatment temperature from TA and DIC
  #T=treatment temp
  
  result2 <- carb(
    flag = 15, 
    var1 = TA_values[i], 
    var2 = DIC1, 
    T = 35, 
    S = 35, 
    P = 1, 
    Pt = 1e-6, 
    Sit = 10e-6, 
    k1k2 = "r", 
    pHscale = "F"
  )
  pH2 <- result2$pH
  pH2_results[i] <- pH2
  
  #Step 3: Calculate final pCO2 from pH2 and TA at 26°C
  result3 <- carb(
    flag = 8, 
    var1 = pH2, 
    var2 = TA_values[i], 
    T = 35, #treatment temp to edit
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
write_xlsx(results_df, path = "carbonate_results_35E_end.xlsx")


