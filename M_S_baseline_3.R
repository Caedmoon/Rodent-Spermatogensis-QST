#packages----
library("dplyr")
library("FME")
library("ggplot2")
library("deSolve")
library("gridExtra")

#notes for development----
#k_breakdown must be equal to the rate of replacement?
#fill in rates , new initial values 
#CSV files for the values so model can be calibrated
#consider that the whole process takes 34 days

#28/05/2025
#Have to alter Nakata_data_transformation so that instead of Cell Count 



#Model Assumptions----
#Post puberty
#34 day length (lectures)
#Levi (2022) found apoptosis was found in dividing germ cells such as spermatogonia and spermatocytes
#Nakata (2015) and Levi (2022) both have cell count per sminiferous tubule section
#Levi - Ki67 assay for proliferating cells
#Nakata (2015) counts different cell types - I will split them into the 4 categories instead of 12 subsections
#In Takashima (2011)-, the number of tubules in a control mouse was 1,321 in mice - will need to multiply by this for total volume
#Ray (2014) contains kinetic parameters identified from dynamic imaigng , cell kinetic studies, trying these first
#Age can have a large effect on production may need a modifier
#middle ground model focusing on Nakata reports
#assumption that every spermatid forms a spermatozoa
#Import observed data----
source("Nakata_data_transformation.R")

# Define the ODE system + wrap----
Model <- function(parms){
  derivs <- function(times, y, parms, fixed) {
    with(as.list(c(y, parms, fixed)), {
      d <- length(initial)
      #calculations
      #feedback to maintain constant levels
      fdbc <- 1
      #Sertoli cells?
      
      # Rate equations
      #Spermatagonium - type A
      d[1] <- k_mitosis * y["Sg_A"] * fdbc - k_trans_A * y["Sg_A"]
      #Spermatagonium - type intermediate
      d[2] <- k_trans_A * y["Sg_A"] * fdbc + k_mitosis * y["Sg_I"] - k_trans_I * y["Sg_I"]
      #Spermatagonium - type B
      d[3] <- k_trans_I * y["Sg_I"] - k_GP * y["Sg_B"]
      #Spermatocytes -PI-L-Z
      d[4] <- k_GP * y["Sg_B"] - k_PILZ * y["Sc_PILZ"]
      #Spermatocytes -P-Di
      d[5] <- k_PILZ * y["Sc_PILZ"] - k_Di * y["Sc_PDi"]
      #Meiotic spermatocytes (Primary and secondary)
      d[6] <- k_Di * y["Sc_PDi"] - k_M * y["Sc_M"]
      #Spermatids 1 - 12
      d[7] <- - k_M * y["Sc_M"] - k_112 * y["St_112"]
      #Spermatids 13-16
      d[8] <- k_112 * y["St_112"] - k_maturation* y["St_1316"]
      #Spermatozoa
      d[9] <- k_maturation * y["St_1316"] - k_breakdown * y["S_zoa"]
      # Return the rates of change
      return(list(d))
    })
  }
  
  # Initial values
  #Spermatagonium A
  Sg_A_0 <- 3.30 #average from Nakata (2015)
  #Spermatagonium I
  Sg_I_0 <- 11.80 #average from Nakata (2015)
  #Spermatagonium B
  Sg_B_0 <- 26.10 #average from Nakata (2015)
  #Spermatocytes -PI-L-Z
  Sc_PILZ_0 <- 51.40 #Nakata (2015) - P1-L-Z avg
  #Spermatocytes -P-Di
  Sc_PDi_0 <- 53.01 #Nakata (2015) - estimated
  #Meiotic spermatocytes (Primary and secondary)
  Sc_M_0 <- 96.30 #Nakata (2015)
  #Spermatids 1 - 12
  St_112_0 <- 144.60 #Nakata (2015)
  #Spermatids 13-16
  St_1316_0 <- 108.80 #Nakata (2015)
  #Spermatozoa
  S_zoa_0 <- 108.8 #assumption that all spermatids form spermatozoa
  initial <- c(Sg_A = Sg_A_0, Sg_I = Sg_I_0, Sg_B = Sg_B_0, Sc_PILZ = Sc_PILZ_0, Sc_PDi = Sc_PDi_0,
               Sc_M = Sc_M_0, St_112 = St_112_0, St_1316 = St_1316_0, S_zoa = S_zoa_0)
  
  #Fixed parameters
  fixed <- list(
  )
  
  # Time sequence for the simulation
  times <- seq(0, 35, by = 0.01)
  
  # Solve the ODE system
  output <- ode(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  #Total per testis? singular
  #output_df[ , paste0(names(output_df)[2:6], ".total")] <- output_df[ , 2:6] * 1321
  return(output_df)
}

#Parameters to fit -----
Time.MOD = 24 #time modifier in hours
parms <- list(
  #all currently from Ray (2014)
  k_mitosis <- 0.002 * Time.MOD, #/h mitosis - determined by division of SCC
  k_trans_A <- 0.004878049 * Time.MOD, #/h S-gonium A --> s-gonium I
  k_trans_I <- 0.01694915 * Time.MOD, #s-gonium I --> S-gonium B 
  k_GP <- 0.04545455 * Time.MOD, #S-gonium B --> S-cyte PI-L-Z 
  k_PILZ <- 0.007142857 * Time.MOD, #/h S-cyte PI-L-Z --> S-cyte P-Di 
  k_Di <- 0.004651163 * Time.MOD, #S-cyte P-Di --> S-cytes M
  k_M <- 0.007142857 * Time.MOD, #/h S-cytes M --> early s-tids via differentiation
  k_112 <- 0.007142857 * Time.MOD, #/h S-tids --> late s-tids via differentiation
  k_maturation <- 0.08 * Time.MOD, #Maturation of spermatids to spermatozoa
  k_breakdown <- 0.002 * Time.MOD #breakdown of spermatozoa
)



#initial estimate -----
Initial_out <- Model(parms = parms)


#Conversion of Observed data for plot ----
obs_all <- bind_rows(obs_v, .id = "Cell")

#plot ----
test <- ggplot() +
  
  # Model output lines
  geom_line(data = Initial_out, aes(x = time, y = Sg_A, colour = "Sg_A")) +
  geom_line(data = Initial_out, aes(x = time, y = Sg_I, colour = "Sg_I")) +
  geom_line(data = Initial_out, aes(x = time, y = Sg_B, colour = "Sg_B")) +
  geom_line(data = Initial_out, aes(x = time, y = Sc_PILZ, colour = "Sc_PILZ")) +
  geom_line(data = Initial_out, aes(x = time, y = Sc_PDi, colour = "Sc_PDi")) +
  geom_line(data = Initial_out, aes(x = time, y = Sc_M, colour = "Sc_M")) +
  geom_line(data = Initial_out, aes(x = time, y = St_112, colour = "St_112")) +
  geom_line(data = Initial_out, aes(x = time, y = St_1316, colour = "St_1316")) +
  geom_line(data = Initial_out, aes(x = time, y = S_zoa, colour = "S_zoa")) +

  # Observed data points
  geom_point(data = obs_all, aes(x = time, y = Count, colour = Cell), shape = 21, fill = "black", size = 2) +
  geom_errorbar(data = obs_all, aes(x = time, ymin = Count - S.D, ymax = Count + S.D, colour = Cell), width = 0.2) +
  
  
  # Color mapping for all elements
  scale_colour_manual(values = c(
    "Sg_A" = "#aec7e8",
    "Sg_I" = "#7f9dc9",
    "Sg_B" = "#4c6cb3",
    "Sc_PILZ" = "#fdb87d",
    "Sc_PDi" = "#e47c19",
    "Sc_M" = "#98df8a",
    "St_112" = "#ff9896",
    "St_1316" = "#d62728",
    "S_zoa" = "#9467bd"  
  )) +
  
  xlab("Time (days)") +
  ylab("Cell count (per tubule)") +
  ggtitle("Mouse Spermatogenesis Model (Detailed)") +
  theme_minimal()

print(test)

#cost function-----
cost <- function(parms) {
  # Run the model with current parameters
  out <- Model(parms = parms)
  
  for (x in seq_along(obs_v_cost)) {
    obs_df <- obs_v_cost[[x]]
    
    # Compute cost
    if (x == 1) {
      cost_value <- modCost(model = out, obs = obs_df, err = "S.D")
    } else {
      cost_value <- modCost(model = out, obs = obs_df, err = "S.D", cost = cost_value)
    }
  }
  
  return(cost_value)
}


# Example of calling the cost function
Modelcost <- cost(parms = parms)

#choose colours for obs data weight

colours <- c(    "Sg_A" = "#aec7e8",
                 "Sg_I" = "#7f9dc9",
                 "Sg_B" = "#4c6cb3",
                 "Sc_PILZ" = "#fdb87d",
                 "Sc_PDi" = "#e47c19",
                 "Sc_M" = "#98df8a",
                 "St_112" = "#ff9896",
                 "St_1316" = "#d62728",
                 "S_zoa" = "#9467bd" 
                 )
plot(Modelcost$residuals$x, Modelcost$residuals$weight, 
     col = colours[Modelcost$residuals$name], 
     pch = 19,  # Solid circle for points
     main = "Weighted Residuals vs. Time", 
     xlab = "Time", 
     ylab = "Weighted Residuals", 
     ylim = range(Modelcost$residuals$weight, na.rm = TRUE))


