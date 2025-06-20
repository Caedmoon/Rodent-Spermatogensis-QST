#Have to do data transformation
#Each species has its own time and cell 
#multiple timepoints with the same value - assuming constant levels 
obs_cell <- read.csv("data/nakata_t1.csv")


i <- 1
obs_v <- list()  # initialize list

for (i in 1:length(obs_cell$Cell)) {
  temp.name <- obs_cell$Cell[i]
  print(temp.name)
  
  # Extract row for current cell
  temp_row <- obs_cell[which(obs_cell$Cell == temp.name), ]
  
  # Generate 3 random timepoints
  timepoints <- runif(3, min = 0, max = 35)
  
  # Repeat the row 3 times and assign new time values
  temp_df <- temp_row[rep(1, 3), ]
  temp_df$time <- timepoints
  
  # Store in list with cell name
  obs_v[[i]] <- temp_df
  names(obs_v)[i] <- temp.name
}

## obs_v for cost

i <- 1
obs_v_cost <- obs_v
for (i in 1:length(obs_v_cost)) {
  new_name <- as.character(obs_v_cost[[i]][1, "Cell"])
  colnames(obs_v_cost[[i]])[colnames(obs_v_cost[[i]]) == "Count"] <- new_name
  obs_v_cost[[i]] <- obs_v_cost[[i]][-1]
  obs_v_cost[[i]] <- obs_v_cost[[i]][, c(3,1,2)]
}
