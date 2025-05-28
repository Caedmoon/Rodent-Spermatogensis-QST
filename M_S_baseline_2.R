#packages----
library("dplyr")
library("FME")
library("ggplot2")
library("deSolve")
library("gridExtra")

#Model Assumptions----
#bsaed on 3D computer simulation but does not account for movement. may need to consider?
#Post puberty
#34 day length (lectures)
#Levi (2022) found apoptosis was found in dividing germ cells such as spermatogonia and spermatocytes
#Nakata (2015) and Levi (2022) both have cell count per sminiferous tubule section
#Levi - Ki67 assay for proliferating cells
#Nakata (2015) counts different cell types - I will split them into the 4 categories instead of 12 subsections
#In Takashima (2011)-, the number of tubules in a control mouse was 1,321 in mice - will need to multiply by this for total volume
#Ray (2014) contains kinetic parameters identified from dynamic imaigng , cell kinetic studies, trying these first
#Age can have a large effect on production may need a modifier
#Elongated spermatid == mature spermatocyte
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
      #Spermatagonium Stem Cell - Type A
      d[1] <- k_mitosis * y["Sg_Sc"] * fdbc - k_l_Sg_D * y["Sg_Sc"]
      #Differentiatng spermatagonium - Type B
      d[2] <- k_mitosis * y["Sg_Sc"] * fdbc + k_Sg_D_mitosis * y["Sg_D"] - k_Sg_D_meiosis * y["Sg_D"] - k_l_Sg_D * y["Sg_D"]
      #P_tene
      d[3] <- k_Sg_D_meiosis * y["Sg_D"] - k_preleptotene * y["P_tene"] - k_l_preleptotene * y["P_tene"]
      #Leptotene
      d[4] <- k_preleptotene * y["P_tene"] - k_leptotene * y["L_tene"] - k_l_leptotene * y["L_tene"]
      #Zygotene
      d[5] <- k_leptotene * y["L_tene"] - k_zygotene * y["Z_tene"] - k_l_zygotene * y["Z_tene"]
      #Pachytene
      d[6] <- k_zygotene * y["Z_tene"] - k_pachytene * y["Pa_tene"] - k_l_pachytene * y["Pa_tene"]
      #Diplotene
      d[7] <- k_pachytene * y["Pa_tene"] - k_diplotene * y["D_tene"] - k_l_diplotene * y["D_tene"]
      #Secondary spermatocytes
      d[8] <- k_diplotene * y["D_tene"] - k_S_cyte * y["S_cyte"] - k_l_S_cyte * y["S_cyte"]
      #Round spermatids
      d[9] <- k_S_cyte * y["S_cyte"] - k_r_spermatid * y["R_tid"] - k_l_r_spermatid * y["R_tid"]
      #Elongated spermatids
      d[10] <- k_r_spermatid * y["R_tid"] - k_l_e_spermatid * y["E_tid"]

      
      
      # Return the rates of change
      return(list(d))
    })
  }
  
  # Initial values
  #Spermatagonium Stem Cell
  Sg_Sc_0 <- 8 #average from Nakata (2015)
  #Differentiating Spermatogonium
  Sg_D_0 <- 12 #average from Nakata (2015)
  #Preleptotene
  P_tene_0 <- 48
  #Leptotene
  L_tene_0 <- 0
  #Zygotene
  Z_tene_0 <- 0
  #Pachytene
  Pa_tene_0 <- 0
  #Diplotene
  D_tene_0 <- 0
  #Secondary spermatocytes
  S_cyte_0 <- 0 #meiotic cell count includes P1 and P2 (96.3 - P1 count) Nakata (2015)
  #Round Spermatid
  R_tid_0 <- 90 # round and elongating spermatid number Nakata (2015)
  #Elongated Spermatid
  E_tid_0 <- 90 # round and elongating spermatid number Nakata (2015)
  #Mature spermatozoa
  S_zoa_0 <- 108.8 #Final maturation stages, best representative Nakata (2015)
  initial <- c(Sg_Sc = Sg_Sc_0, Sg_D = Sg_D_0, P_tene = P_tene_0, L_tene = L_tene_0, Z_tene = Z_tene_0,
               Pa_tene = Pa_tene_0, D_tene = D_tene_0, S_cyte = S_cyte_0, R_tid = R_tid_0, E_tid = E_tid_0)
  
  #Fixed parameters
  fixed <- list()
  #all currently from Ray (2014)
  k_Sg_D_mitosis <- 1/205 * parms[["Time.MOD"]]
  k_preleptotene <- 1/44 * parms[["Time.MOD"]] #/h rate of differentiation of preleptotene
  k_leptotene <- 1/24 * parms[["Time.MOD"]] #/h rate of differentiation of leptotene
  k_zygotene <- 1/42 * parms[["Time.MOD"]] #/h rate of differentiation of zygotene
  k_pachytene <- 1/165 * parms[["Time.MOD"]] #/h rate of differentiation of leptotene
  k_r_spermatid <- 1/140 * parms[["Time.MOD"]] #/h rate of differentiation of round spermatid
  k_mitosis <- k_preleptotene * parms[["Time.MOD"]] #/h rate of spermatogonium stem cell division
  k_Sg_D_meiosis <- 1/88 * parms[["Time.MOD"]] #/h rate of differentiating spermatogonium division
  k_diplotene <- 1/24 * parms[["Time.MOD"]] #/h division time of diplotene
  k_S_cyte <- 1/22 * parms[["Time.MOD"]] #/h division time of secondary spermatocyte
  k_l_Sg_D <- 1/206 * parms[["Time.MOD"]] # /h lifespan-based rate of differentiating spermatogonium (206–206 h)
  k_l_preleptotene <- 1/((43 + 93)/2) * parms[["Time.MOD"]]   # /h lifespan-based rate of preleptotene (43–93 h)
  k_l_leptotene <- 1/((23 + 73)/2) * parms[["Time.MOD"]]       # /h lifespan-based rate of leptotene (23–73 h)
  k_l_zygotene <- 1/((41 + 91)/2) * parms[["Time.MOD"]]        # /h lifespan-based rate of zygotene (41–91 h)
  k_l_pachytene <- 1/((164 + 214)/2) * parms[["Time.MOD"]]     # /h lifespan-based rate of pachytene (164–214 h)
  k_l_diplotene <- 1/((21 + 23)/2) * parms[["Time.MOD"]]       # /h lifespan-based rate of diplotene (21–23 h)
  k_l_S_cyte <- 1/((22 + 25)/2) * parms[["Time.MOD"]]          # /h lifespan-based rate of secondary spermatocyte (22–25 h)
  k_l_r_spermatid <- 1/((139 + 189)/2) * parms[["Time.MOD"]]   # /h lifespan-based rate of round spermatid (139–189 h)
  k_l_e_spermatid <- 1/((200 + 230)/2) * parms[["Time.MOD"]]   # /h lifespan-based rate of elongated spermatid (200–230 h)
  # Time sequence for the simulation
  times <- seq(0, 100, by = 0.01)
  
  # Solve the ODE system
  output <- ode(y = initial, times = times, func = derivs, parms = parms, fixed = fixed)
  
  # Convert output to a data frame for easier visualization
  output_df <- as.data.frame(output)
  #Total per testis? singular
  #cell_cols <- setdiff(names(output_df), "time")
  #output_df[ , paste0(cell_cols, ".total")] <- output_df[ , cell_cols] * 1321
  return(output_df)
}

#Parameters to fit -----
parms <- list(
  Time.MOD = 1
)



#initial estimate -----
Initial_out <- Model(parms = parms)


test <- ggplot() +
  geom_line(data = Initial_out, aes(x = time, y = Sg_Sc, colour = "SSC")) +
  geom_line(data = Initial_out, aes(x = time, y = Sg_D, colour = "Differentiating spermatogonia")) +
  geom_line(data = Initial_out, aes(x = time, y = P_tene, colour = "Preleptotene")) +
  geom_line(data = Initial_out, aes(x = time, y = L_tene, colour = "Leptotene")) +
  geom_line(data = Initial_out, aes(x = time, y = Z_tene, colour = "Zygotene")) +
  geom_line(data = Initial_out, aes(x = time, y = Pa_tene, colour = "Pachytene")) +
  geom_line(data = Initial_out, aes(x = time, y = D_tene, colour = "Diplotene")) +
  geom_line(data = Initial_out, aes(x = time, y = S_cyte, colour = "Secondary spermatocytes")) +
  geom_line(data = Initial_out, aes(x = time, y = R_tid, colour = "Round spermatids")) +
  geom_line(data = Initial_out, aes(x = time, y = E_tid, colour = "Elongated spermatids")) +
  scale_colour_manual(values = c(
    "SSC" = "#1f77b4",
    "Differentiating spermatogonia" = "#17becf",
    "Preleptotene" = "#aec7e8",
    "Leptotene" = "#ff7f0e",
    "Zygotene" = "#ffbb78",
    "Pachytene" = "#2ca02c",
    "Diplotene" = "coral",
    "Secondary spermatocytes" = "#98df8a",
    "Round spermatids" = "#d62728",
    "Elongated spermatids" = "#9467bd"
  )) +
  xlab("Time (hr)") +
  ylab("Cell count (per tubule)") +
  ggtitle("Mouse Spermatogenesis Model Output") +
  theme_minimal()

print(test)

# 