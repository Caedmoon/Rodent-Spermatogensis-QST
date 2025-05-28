#packages----
library("dplyr")
library("FME")
library("ggplot2")
library("deSolve")
library("gridExtra")

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
#Import observed data----
# obs_A <- readRDS("Mathematical modelling/obs_A_example")
# obs_B <- readRDS("Mathematical modelling/obs_B_example")
# obs_C <- readRDS("Mathematical modelling/obs_C_example")
# obs_V <- list(obs_A,obs_B,obs_C)
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
      #Spermatagonium
      d[1] <- k_mitosis * y["Sg"] * fdbc - k_GP * y["Sg"]
      #Primary spermatocytes
      d[2] <- k_GP * y["Sg"] - k_meiosisI * y["P_cyte"]
      #Secondary spermatocytes
      d[3] <- k_meiosisI * y["P_cyte"] - k_meiosisII * y["S_cyte"]
      #Spermatids
      d[4] <- k_meiosisII * y["S_cyte"] - k_maturation * y["S_tid"]
      #Mature sperm
      d[5] <- k_maturation * y["S_tid"] - k_breakdown * y["S_zoa"]
      
      
      # Return the rates of change
      return(list(d))
    })
  }
  
  # Initial values
  #Spermatagonium
  Sg_0 <- 41.2 #average from Nakata (2015)
  #Primary spermatocytes
  P_cyte_0 <- 51.4 #Nakata (2015) - P1-L-Z avg
  #Secondary spermatocytes
  S_cyte_0 <- 44.9 #meiotic cell count includes P1 and P2 (96.3 - P1 count) Nakata (2015)
  #Spermatids
  S_tid_0 <- 144.6 # round and elongating spermatid number Nakata (2015)
  #Mature spermatozoa
  S_zoa_0 <- 108.8 #Final maturation stages, best representative Nakata (2015)
  initial <- c(Sg = Sg_0, P_cyte = P_cyte_0, S_cyte = S_cyte_0, S_tid = S_tid_0, S_zoa = S_zoa_0)
  
  #Fixed parameters
  fixed <- list(
    #all currently from Ray (2014)
    k_mitosis <- 0.004545455 * parms[["Time.MOD"]], #/h mitosis - determined by division of SCC
    k_GP <- 0.004878049 * parms[["Time.MOD"]], #/h S-gonium --> P-cyte differentiation
    k_meiosisI <- 0.01694915 * parms[["Time.MOD"]], #P-cyte --> S-cyte via meiosis I, may need to be calibrated as it is using average of P,L,Z,P,D 
    k_meiosisII <- 0.04545455 * parms[["Time.MOD"]], #S-cyte --> S-tid via meiosis II
    k_maturation <- 0.007142857 * parms[["Time.MOD"]], #/h S-tid --> s-zoa via differentiation
    k_breakdown <- 0.004651163 * parms[["Time.MOD"]] #average lifespan of elongated spermatid
  )
  
  # Time sequence for the simulation
  times <- seq(0, 75, by = 0.01)
  
  # Solve the ODE system
  output <- ode(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  #Total per testis? singular
  output_df[ , paste0(names(output_df)[2:6], ".total")] <- output_df[ , 2:6] * 1321
  return(output_df)
}

#Parameters to fit -----
parms <- list(
  Time.MOD = 24
)



#initial estimate -----
Initial_out <- Model(parms = parms)


test <- ggplot() +
  geom_line(data = Initial_out, aes(x = time, y = Sg, colour = "Spermatogonia")) +
  geom_line(data = Initial_out, aes(x = time, y = P_cyte, colour = "Primary spermatocytes")) +
  geom_line(data = Initial_out, aes(x = time, y = S_cyte, colour = "Secondary spermatocytes")) +
  geom_line(data = Initial_out, aes(x = time, y = S_tid, colour = "Spermatids")) +
  geom_line(data = Initial_out, aes(x = time, y = S_zoa, colour = "Spermatozoa")) +
  scale_colour_manual(values = c(
    "Spermatogonia" = "#1f77b4",
    "Primary spermatocytes" = "#ff7f0e",
    "Secondary spermatocytes" = "#2ca02c",
    "Spermatids" = "#d62728",
    "Spermatozoa" = "#9467bd"
  )) +
  xlab("Time (days)") +
  ylab("Cell count (per tubule)") +
  ggtitle("Mouse Spermatogenesis Model") +
  theme_minimal()

print(test)

# 
# #cost function-----
# cost <- function(parms) {
#   # Run the model with the provided parameters
#   out <- Model(parms = parms)
#   # Initialize cost using observations from obs_A
#   cost_value <- modCost(out, obs = obs_A, err = "sd")
# 
#   # Add costs from obs_B
#   cost_value <- modCost(out, obs = obs_B, err = "sd", cost = cost_value)
# 
#   # Add costs from obs_C and return the final cost
#   cost_value <- modCost(out, obs = obs_C, err = "sd", cost = cost_value)
# 
#   return(cost_value)
# }
# 
# # Example of calling the cost function
# Modelcost <- cost(parms = parms)
# 
# #choose colours for obs data weight
# 
# colours <- c("A" = "blue", "B" = "red", "C" = "green")
# 
# plot(Modelcost$residuals$x, Modelcost$residuals$weight, 
#      col = colours[Modelcost$residuals$name], 
#      pch = 19,  # Solid circle for points
#      main = "Weighted Residuals vs. Time", 
#      xlab = "Time", 
#      ylab = "Weighted Residuals", 
#      ylim = range(Modelcost$residuals$weight, na.rm = TRUE))
# 
# 
# #Local sensitivity Analysis ----
# 
# LSA <- sensFun(Model, parms = parms)
# 
# #summary plot
# plot(summary(LSA, vars = TRUE))
# 
# #Bivariate sensitivity 
# pairs(LSA, which = c(LSA$var), col = colours)
# 
# #Sensitivity to X parameters overtime 
# par(mfrow = c(1,3))
# #A
# plot(LSA$x[LSA$var == "A"], LSA$k1[LSA$var == "A"], type = "l", col = "blue", 
#      ylim = range(-1,1), 
#      ylab = "Sensitivity of A to Parameters", xlab = "Time", 
#      main = "Sensitivity of Measured Variable A",
#      lwd = 8)
# lines(LSA$x[LSA$var == "A"], LSA$k2[LSA$var == "A"], col = "red", lwd = 8)
# legend("topright", legend = c("K1", "K2"), 
#        col = c("blue", "red"), lty = 1, lwd = 8)
# #B
# plot(LSA$x[LSA$var == "B"], LSA$k1[LSA$var == "B"], type = "l", col = "blue", 
#      ylim = range(-1,1), 
#      ylab = "Sensitivity of B to Parameters", xlab = "Time", 
#      main = "Sensitivity of Measured Variable B", lwd = 8)
# lines(LSA$x[LSA$var == "B"], LSA$k2[LSA$var == "B"], col = "red", lwd = 8)
# legend("topright", legend = c("K1", "K2"), 
#        col = c("blue", "red"), lty = 1, lwd = 8)
# #C
# plot(LSA$x[LSA$var == "C"], LSA$k1[LSA$var == "C"], type = "l", col = "blue", 
#      ylim = range(-1,1), 
#      ylab = "Sensitivity of C to Parameters", xlab = "Time", 
#      main = "Sensitivity of Measured Variable C", lwd = 8)
# lines(LSA$x[LSA$var == "C"], LSA$k2[LSA$var == "C"], col = "red", lwd = 8)
# legend("topright", legend = c("K1", "K2"), 
#        col = c("blue", "red"), lty = 1, lwd = 8)
# par(mfrow = c(1,1))
# 
# #parameter identifiability analysis----
# 
# ident <- collin(sensfun = LSA)
# 
# print(ident)
# 
# plot(ident, log = "y")
# 
# #Model parameterisation ----
# Modelcost2 <- function(parms) {
#   # Fit the model using the cost function
#   Fit <- modFit(f = cost, p = unlist(parms))
#   
#   # Return the fit object or any relevant output
#   return(Fit)
# }
# 
# Fit <- Modelcost2(parms = parms)
# 
# Fitted_parms <- Fit$par
# 
# #Calibrated model output----
# Fit_out <- Model(parms = Fitted_parms)
# 
# #Plotting comparison of original vs calibrated output-----
# variable <- c("A", "B", "C")
# Fit_plot_v <- vector("list", length(variable))
# l_size = 1.5  # Adjust this value for line thickness
# d_size = 3    # Adjust this value for point size
#   Fit_plot_v[[1]] <- ggplot() +
#     geom_line(data = Initial_out, 
#               aes(x = time, y = A, colour = "A", linetype = "Initial"), size = l_size) +
#     geom_line(data = Fit_out, 
#               aes(x = time, y = A, colour = "A", linetype = "Calibrated"), size = l_size) +
#     
#     geom_point(data = obs_A, 
#                aes(x = time, y = A, colour = "A"), size = d_size) +
#     geom_errorbar(data = obs_A, 
#                   aes(x = time, 
#                       ymin = A - sd, 
#                       ymax = A + sd, 
#                       colour = "A"), 
#                   width = 0.2, size = l_size) +  # Adjust size for error bars too
#     ggtitle(paste("Change in chemical A over time")) +  # Concatenating the title into a single string
#     xlab("Time") +
#     ylab("Amount") +
#     scale_colour_manual(values = c("A" = "blue", "B" = "red", "C" = "green")) +
#     scale_linetype_manual(values = c("Initial" = "dashed", "Calibrated" = "solid")) +
#     labs(colour = "Chemical")  # Set the color legend title
# 
#   Fit_plot_v[[2]] <- ggplot() +
#     geom_line(data = Initial_out, 
#               aes(x = time, y = B, colour = "B", linetype = "Initial"), size = l_size) +
#     geom_line(data = Fit_out, 
#               aes(x = time, y = B, colour = "B", linetype = "Calibrated"), size = l_size) +
#     
#     geom_point(data = obs_B, 
#                aes(x = time, y = B, colour = "B"), size = d_size) +
#     geom_errorbar(data = obs_B, 
#                   aes(x = time, 
#                       ymin = B - sd, 
#                       ymax = B + sd, 
#                       colour = "B"), 
#                   width = 0.2, size = l_size) +  # Adjust size for error bars too
#     ggtitle(paste("Change in chemical B over time")) +  # Concatenating the title into a single string
#     xlab("Time") +
#     ylab("Amount") +
#     scale_colour_manual(values = c("A" = "blue", "B" = "red", "C" = "green")) +
#     scale_linetype_manual(values = c("Initial" = "dashed", "Calibrated" = "solid")) +
#     labs(colour = "Chemical")  # Set the color legend title
#   
#   Fit_plot_v[[3]] <- ggplot() +
#     geom_line(data = Initial_out, 
#               aes(x = time, y = C, colour = "C", linetype = "Initial"), size = l_size) +
#     geom_line(data = Fit_out, 
#               aes(x = time, y = C, colour = "C", linetype = "Calibrated"), size = l_size) +
#     
#     geom_point(data = obs_C, 
#                aes(x = time, y = C, colour = "C"), size = d_size) +
#     geom_errorbar(data = obs_C, 
#                   aes(x = time, 
#                       ymin = C - sd, 
#                       ymax = C + sd, 
#                       colour = "C"), 
#                   width = 0.2, size = l_size) +  # Adjust size for error bars too
#     ggtitle(paste("Change in chemical C over time")) +  # Concatenating the title into a single string
#     xlab("Time") +
#     ylab("Amount") +
#     scale_colour_manual(values = c("A" = "blue", "B" = "red", "C" = "green")) +
#     scale_linetype_manual(values = c("Initial" = "dashed", "Calibrated" = "solid")) +
#     labs(colour = "Chemical")  # Set the color legend title
# # Arrange the plots in a grid
# grid.arrange(grobs = Fit_plot_v, nrow = 1, ncol = length(variable))