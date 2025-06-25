##############################################################
## Code adapted from: 
## Forrest, S.W., Pagendam, D., Bode, M., Drovandi, C., Potts, J.R., 
## Perry, J., Vanderduys, E., & Hoskins, A.J. (2024). 
## Predicting fine-scale distributions and emergent spatiotemporal patterns 
## from temporally dynamic step selection simulations. *Ecography*. 
## https://doi.org/10.1111/ecog.07421
##
## This adaptation includes extensions to incorporate habitat-specific rasters,
## memory-informed revisitation dynamics, and hourly coefficient-driven simulation.
##############################################################


##############################################
############# 1 Pair of harmonics ###########
##############################################

rm(list = ls())  # Clear the R environment

# Load required packages
library(sf)
library(tidyverse)
library(terra)

# Additional simulation-related packages
packages <- c("amt", "lubridate", "terra", "raster", "tictoc", "matrixStats", "Rfast")
walk(packages, require, character.only = T)

# Set working directory to location of harmonic coefficients
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair")

# Load coefficient estimates from model fitting output (for a single harmonic pair)
hourly_coefs <- read_csv("hours_temporal_coefficients_1hrm.csv")
head(hourly_coefs)

###########
#CODE ADDED
# Visualisation of how coefficients change with hour, particularly focusing on memory term
hourly_coefs_long <- hourly_coefs %>% 
  pivot_longer(cols = -c(hour, hour_transformed) , names_to = "covariate", values_to = "coef")

ggplot(hourly_coefs_long %>% filter(
  covariate %in% c(
    # "disturb_distance", 
    # "habitat_distance"#, 
    "memory"
  )), 
  aes(x = hour, y = coef, colour = covariate)) +
  geom_point() +
  geom_line() +
  labs(title = "Hourly coefficients", x = "Hour", y = "Coefficient")
###########

# Set working directory for SSF data (observed quoll movements and samples)
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")

# Load observed and available SSF points
all_ssf1 <-read.csv("Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")

# Filter for breeding season and select only true steps
all_ssf_breeding <- all_ssf1 %>% 
  filter(season == "breeding") %>% 
  arrange(t1_)
all_ssf_breeding_true <- all_ssf_breeding %>% 
  filter(case_ == "1") %>% 
  arrange(t1_)

#############
## Update the required landscape rasters below as needed (i.e., non-mining, dispersed, aggregated). It is coded for current mining at the moment.
#############

# Load NDVI and distance rasters for the current mining landscape
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/current mining")
ndvi <- rast("Current_habitat_classifications.tif")
habitat_dist <- rast("Current_habitat_distances.tif")
disturb_dist <- rast("Current_disturb_distances.tif")

ndvi  # Show NDVI raster info

# Resample distance rasters to match NDVI grid (resolution and extent)
habitat_dist_aligned <- terra::resample(habitat_dist, ndvi, method='bilinear') # or 'near' depending on your data
disturb_dist_aligned <- terra::resample(disturb_dist, ndvi, method='bilinear') # or 'near' depending on your data

## Log-transform the aligned distance rasters to reduce skew
constant = 0.1
habitat_dist_aligned_log_transformed <- log(habitat_dist_aligned+constant)
plot(habitat_dist_aligned_log_transformed)

constant = 0.1
disturb_dist_aligned_log_transformed <- log(disturb_dist_aligned+constant)
plot(disturb_dist_aligned_log_transformed)

#Split NDVI categorical raster into binary layers representing each habitat class

# Create binary maps for each habitat category
dense_veg_map <- ifel(ndvi == 1, 1, 0)
#plot(dense_veg_map)

grassland_map <- ifel(ndvi == 2, 1, 0)
#plot(grassland_map)

rocky_map <- ifel(ndvi == 3, 1, 0)
#plot(rocky_map)

other_disturbed_map <- ifel(ndvi == 4, 1, 0)
#plot(other_disturbed_map)

mine_pits_waste_dumps_map <- ifel(ndvi == 5, 1, 0)
#plot(mine_pits_waste_dumps_map)

# setting up for naive preds ----------------------------------------------

# combine resources into stack
resources <- c(dense_veg_map, grassland_map, rocky_map, other_disturbed_map, mine_pits_waste_dumps_map, habitat_dist_aligned_log_transformed, disturb_dist_aligned_log_transformed)

# get the extent of the resources
xmin <- ext(resources[[1]])[1]
xmax <- ext(resources[[1]])[2]
ymin <- ext(resources[[1]])[3]
ymax <- ext(resources[[1]])[4]

# set the origin to (0,0) so that simulation works in relative space
ext(resources) <- c(xmin - xmin, xmax - xmin, ymin - ymin, ymax - ymin)

#plot(resources)

hourly_coefs$true_hour<-hourly_coefs$hour  # Preserve original hour column

# Define the original hours in order (in decimal time)
original_hours <- c(18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)

# Create a numeric sequence to index those hours
new_sequence <- seq(1, 24, by=1)

# Map original hours to new sequence
hour_mapping <- setNames(new_sequence, original_hours)

# Apply the mapping to hourly_coefs
hourly_coefs$hour <- hour_mapping[as.character(hourly_coefs$hour)]

# Check updated hourly coefficients
hourly_coefs

# naive predictions -------------------------------------------------------

naive_pred_stack <- c()         # Initialise list for raster predictions
naive_pred_norm_stack <- c()    # (Optional) Normalised predictions if needed later

# Loop through each hour and generate logit-scale raster predictions
for(i in seq(1, 24, by = 1)) {
  
  # Ensure hour exists in coefficient table
  if(i %in% hourly_coefs$hour) {
    
    # Create prediction raster for each habitat/resource covariate
    dense_veg_lin <- resources[[1]]
    dense_veg_pred <- dense_veg_lin * hourly_coefs$ndvidense_veg[which(hourly_coefs$hour == i)]
    
    grassland_lin <- resources[[2]]
    grassland_pred <- grassland_lin * hourly_coefs$ndvigrassland[which(hourly_coefs$hour == i)]
    
    rocky_lin <- resources[[3]]
    rocky_pred <- rocky_lin * hourly_coefs$ndvirocky[which(hourly_coefs$hour == i)]
    
    other_disturbed_lin <- resources[[4]]
    other_disturbed_pred <- other_disturbed_lin * hourly_coefs$ndviother_disturbed[which(hourly_coefs$hour == i)]
    
    mine_pits_waste_dumps_lin <- resources[[5]]
    mine_pits_waste_dumps_pred <- mine_pits_waste_dumps_lin * hourly_coefs$ndvimine_pit_waste_dump[which(hourly_coefs$hour == i)]
    
    habitat_lin <- resources[[6]]
    habitat_pred <- habitat_lin * hourly_coefs$habitat_distance[which(hourly_coefs$hour == i)]
    
    disturb_lin <- resources[[7]]
    disturb_pred <- disturb_lin * hourly_coefs$disturb_distance[which(hourly_coefs$hour == i)]
    
    # Combine habitat predictions into one logit-scale raster
    naive_pred_logit <- dense_veg_pred +
      grassland_pred +
      rocky_pred +
      other_disturbed_pred +
      mine_pits_waste_dumps_pred +
      habitat_pred +
      disturb_pred
    
    # Append raster layer to stack
    naive_pred_stack <- c(naive_pred_stack, naive_pred_logit)
    
  }
}

# Convert list of rasters to SpatRaster stack
naive_predictions <- rast(naive_pred_stack)
#plot(naive_predictions)  # Optional visual check

# simulation function -----------------------------------------------------
# model_outputs<-readRDS("C:/Users/Mitch/OneDrive/Desktop/Simulation quoll paper/Step selection/Model outputs/model_twostep_2p_harms.RDS")
# summary(model_outputs)

# Define SSF simulation function with no memory (pure habitat-driven choice)
simulate_ssf_memory <- function(n_steps, n_ch, coef, xy0, resc) { #, memory_period, memory_delay, spatial_sd, hour, warmup, boundary
  # , exp_gamma
  
  tic()  # Start timing
  
  # Sample step lengths from gamma distribution (per hour)
  sl <- rgamma(n_steps * n_ch,
               shape = rep(coef$shape, each = n_ch),
               scale = rep(coef$scale, each = n_ch))
  
  # Turning angle parameters (kappa and mu depending on persistence)
  ta_coefs <- rep(coef$kappa, each = n_ch, length.out = n_steps)
  ta_coefs_mu <- ifelse(ta_coefs > 0, pi, 0)
  ta <- as.vector(mapply(Rfast::rvonmises, n = n_ch, m = ta_coefs_mu, k = abs(ta_coefs))) - pi
  
  # Setup sequences for simulation steps
  steps <- rep(1:n_steps, each = n_ch)
  
  ###########
  #CODE ADDED
  # Create hour values to be saved with the output
  hour <- rep(seq(1,24,by=1), length.out = n_steps)
  ###########
  
  # Initialise state variables for movement
  x_0 <- xy0[[1]]
  y_0 <- xy0[[2]]
  angle_0 <- runif(1, min = -pi, max = pi)
  
  x <- rep(NA, n_steps)
  y <- rep(NA, n_steps)
  step_length <- rep(NA, n_steps)
  angle <- rep(NA, n_steps)
  bearing <- rep(NA, n_steps)
  hab_p <- rep(NA, n_steps)
  mem_p <- rep(NA, n_steps)  # Not used here but left for consistency
  
  # Set first location and initial heading
  x[1] <- x_0
  y[1] <- y_0
  step_length[1] <- 0
  angle[1] <- angle_0
  bearing[1] <- 0
  hab_p[1] <- 0
  mem_p[1] <- 0
  
  for (i in 2:n_steps) {
    
    # Propose candidate steps from previous location
    x_prop <- (x[i - 1] + sl[steps == i] * cos(bearing[i - 1] + ta[steps == i]))
    y_prop <- (y[i - 1] + sl[steps == i] * sin(bearing[i - 1] + ta[steps == i]))
    sl_prop <- sl[steps == i]
    angle_prop <- ta[steps == i]
    bearing_prop <- atan2(y_prop - y[i - 1], x_prop - x[i - 1])
    
    ###########
    #CODE ADDED
    # Ensure correct raster index for current step (loops through naive prediction layers)
    resc_index <- ifelse(i %% terra::nlyr(naive_predictions) == 0, 
                         terra::nlyr(naive_predictions), 
                         i %% terra::nlyr(naive_predictions))
    ###########
    
    # Extract raster values at candidate endpoints
    p <- terra::extract(resc[[resc_index]], cbind(x_prop, y_prop))[,1]
    
    # Replace missing predictions with a low logit value
    p[is.na(p)] <- -200
    
    # Choose one location probabilistically based on habitat suitability
    w <- sample(n_ch, 1, prob = exp(p))
    
    # Update states
    x[i] <- x_prop[w]
    y[i] <- y_prop[w]
    step_length[i] <- sl_prop[w]
    angle[i] <- angle_prop[w]
    bearing[i] <- bearing_prop[w]
    hab_p[i] <- p[w]
    #mem_p[i] <- memory_p[w]  # Memory not used here
  }
  
  (toc())  # End timing
  
  # Return dataframe of steps and attributes
  return(data_frame(x = x, 
                    y = y, 
                    ###########
                    #CODE ADDED
                    hour = hour,
                    ###########
                    sl = step_length, 
                    angle = angle, 
                    bearing = bearing, 
                    hab_p = hab_p 
                    #mem_p = mem_p
  ))
}


# preparing to run function -----------------------------------------------

# Optional memory model parameters were here in older versions
# They are commented out since this is a memory-off version

# load(file = "C:/Users/Mitch/OneDrive/Desktop/Simulation quoll paper/Step selection/optim_space_time_KDE_param_list_mem1000_2024-02-05.RData")
# space_time_param_list_mem1000_ids <- space_time_KDE_param_list_mem1000[c(1,  2,  3,  5,  6,  8,  9, 11, 12, 13)]

# # spatial sd values
# spatial_sd_ids <- sapply(1:length(quoll_ids), function(i) {space_time_param_list_mem1000_ids[[i]]$par[1]})
# spatial_sd_ids

# # temporal decay values
# exp_gamma_ids <- sapply(1:length(quoll_ids), function(i) {space_time_param_list_mem1000_ids[[i]]$par[2]})
# exp_gamma_ids

# Previous version for extracting true individual-based starting locations
# start_x <- sapply(quoll_ids, function(id_val) {quoll_data_all %>% filter(y == 1 & id == id_val) %>% slice(1) %>% pull(x1) - xmin})
# start_y <- sapply(quoll_ids, function(id_val) {quoll_data_all %>% filter(y == 1 & id == id_val) %>% slice(1) %>% pull(y1) - ymin})
# start_quoll <- matrix(c(start_x, start_y), ncol = 2, nrow = 10)

# Alternative method for sampling uniformly from raster extent (commented)
# # for random starting locations
# start_unif <- cbind(runif(100,
#                           min = ext(naive_predictions[[1]])[1],
#                           max = ext(naive_predictions[[1]])[2]),
#                     runif(100,
#                           min = ext(naive_predictions[[1]])[3],
#                           max = ext(naive_predictions[[1]])[4]))

# Save for reuse (already done)
# save(start_unif,file="C:/Users/Mitch/OneDrive/Desktop/Simulation quoll paper/Step selection/quoll_starting_locations_landscape_centric.RData")

# Load pre-generated landscape-centric starting locations
load(file = "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/quoll_starting_locations_landscape_centric.RData")

# Adjustments for coordinate alignment were once done here, but commented
# start_unif[,1] <- start_unif[,1] - xmin
# start_unif[,2] <- start_unif[,2] - ymin

# Plot naive raster with starting locations
plot(naive_predictions[[1]])
points(x = start_unif[,1], y = start_unif[,2])

# Use first 50 starting coordinates
start_quoll <- matrix(c(x = start_unif[1:50,1], y = start_unif[1:50,2]), ncol = 2, nrow = 50)

# Check format
str(start_quoll)

# Prepare NDVI raster for ggplot if needed
ndvi_df <- as.data.frame(ndvi, xy = TRUE)
str(ndvi_df)
#str(sim_data_mem)

aes(x = x, y = y, fill = factor("Mining landscape single_large_version"))

#quoll_ids<-1:100

# Set simulation parameters
n_indvs <- 1     # One replicate per individual
n_steps <- 8760  # 12-month simulation (one step per hour)
n_ch <- 50       # Number of candidate steps per move
coef <- hourly_coefs
resc <- naive_predictions
# memory_period <- ...
# memory_delay <- ...
# spatial_sd <- ...
# exp_gamma <- ...

###########
#CODE ADDED
# For plotting, make sure you account for (xmin, ymin) if plotting true locations
plot(ndvi)
# Add points using original full coordinates
points(start_quoll[, 1] + xmin, start_quoll[, 2] + ymin, col = 'red')
###########

##Make sure to rename the memory column in excel
quoll_ids<-1:50  # Simulate 50 individuals

gc()       # Garbage collection to clear memory
dev.off()  # Clear plots

# Begin main simulation loop
for(k in quoll_ids) {
  
  #k=1
  stps_mem_list <- vector(mode = "list", length = n_indvs)  # Create list to hold paths
  
  tic()  # Start timing for this individual
  
  for(j in 1:n_indvs) {
    
    #j=1
    
    # Run simulation for this replicate
    stps_mem_list[[j]] <- simulate_ssf_memory(n_steps = n_steps,
                                              n_ch = n_ch,
                                              coef = hourly_coefs,
                                              xy0 = start_quoll[k,],
                                              resc = naive_predictions
                                              #memory_period = memory_period,
                                              #memory_delay = memory_delay,
                                              #spatial_sd = 689.4857
                                              # exp_gamma = exp_gamma_ids[k]
    )
  }
  
  toc()  # End timing
  
  # Convert to amt-format track
  animals_mem <- data_frame(
    id = paste0("id", k, 1:n_indvs),
    track = map(stps_mem_list, ~ amt::track(
      x = .$x, 
      y = .$y, 
      t = ymd_hm("2021-06-01 00:00") + hours(1:nrow(.)), 
      angle = .$angle,
      bearing = .$bearing,
      hab_p = .$hab_p#,
      #mem_p = .$mem_p
    )))
  
  # Flatten the nested tibble
  sim_data_mem <- unnest(animals_mem, cols = track)
  
  # Recode 1am to hour 24
  sim_data_mem$hour <- ifelse(hour(sim_data_mem$t_) == 1, 24, hour(sim_data_mem$t_))
  
  ###########
  #CODE ADDED
  # Adjust coordinates for export or plotting (back to real coordinates)
  plot(naive_predictions[[1]])
  points(sim_data_mem %>% pull(x_), 
         sim_data_mem %>% pull(y_))
  
  sim_data_mem_truexy <- sim_data_mem %>% mutate(x_ = x_ + xmin, y_ = y_ + ymin)
  ###########
  
  # Optional plotting block (commented out)
  # print(ggplot() +
  #         geom_raster(data = ndvi_df,
  #                     aes(x = x, y = y, fill = factor(Current_habitat_classifications)),
  #                     alpha = 0.5) +
  #         geom_point(data = sim_data_mem,
  #                    aes(x = x_, y = y_, colour = id),
  #                    size = 0.75, alpha = 0.1) +
  #         geom_path(data = sim_data_mem,
  #                   aes(x = x_, y = y_, colour = id),
  #                   alpha = 0.75) +
  #         scale_fill_brewer(palette = "Greens") +
  #         scale_color_viridis_d("Sim ID") +
  #         coord_equal() +
  #         theme_classic() +
  #         theme(legend.position = "none"))
  
  dev.off()  # Close any open graphical devices
  
  # Output path and filename
  base_path <- "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Simulation trajectories/Landscape-centric/Original landscape/"
  filename <- paste0( base_path,"1hrm", "id_",  k, "_EvM_mem", "_", n_indvs, "ind_", n_ch, "ch_", n_steps, "stp_", Sys.Date(), ".csv")
  print(filename)
  
  # Save results
  write_csv(sim_data_mem_truexy, filename)
  
  print(paste("Completed:", k, "of", 50))  # Progress log
  
  gc()  # Clean up memory
}

##################################################################
### Repeat for each landscape type with the respective layers ####
##################################################################