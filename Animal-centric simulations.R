##############################################################
## Code adapted from: 
## Forrest, S.W., Pagendam, D., Bode, M., Drovandi, C., Potts, J.R., 
## Perry, J., Vanderduys, E., & Hoskins, A.J. (2024). 
## Predicting fine-scale distributions and emergent spatiotemporal patterns 
## from temporally dynamic step selection simulations. *Ecography*. 
## https://doi.org/10.1111/ecog.07421
##############################################################



##############################################
############# 1 Pair of harmonics ############
##############################################

rm(list = ls())  # Clear environment

# Load spatial, tidyverse, and raster handling packages
library(sf)
library(tidyverse)
library(terra)

# Load additional required packages for movement simulation and optimisation
packages <- c("amt", "lubridate", "terra", "raster", "tictoc", "matrixStats", "Rfast")
walk(packages, require, character.only = T)

# Set working directory to where the harmonic coefficients are stored
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair")

##################################
## Change to different harmonics as required (i.e., 0hrm, 1hrm, 2hrm, 3hrm)
##################################

# Load coefficient estimates for the iSSF model for one pair of harmonics
hourly_coefs <- read_csv("hours_temporal_coefficients_1hrm.csv")
head(hourly_coefs)

###########
#CODE ADDED
# just having a look at some of the temporally dynamic coefs
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

# Set working directory for data
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")

# Load observed quoll movement data and habitat availability points
all_ssf1 <-read.csv("Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")

# Filter for breeding season and sort by timestamp
all_ssf_breeding <- all_ssf1 %>% 
  filter(season == "breeding") %>% 
  arrange(t1_)
all_ssf_breeding_true <- all_ssf_breeding %>% 
  arrange(t1_)

# Define 100 quoll IDs to simulate from
quoll_ids <- c(1:100)

#############
## Update the required landscape rasters below as needed (i.e., non-mining, dispersed, aggregated). It is coded for current mining at the moment.
#############

# Load habitat rasters for current mining landscape
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/current mining")
ndvi <- rast("Current_habitat_classifications.tif")
habitat_dist <- rast("Current_habitat_distances.tif")
disturb_dist <- rast("Current_disturb_distances.tif")

# Resample habitat and disturbance distance rasters to match NDVI resolution
habitat_dist_aligned <- terra::resample(habitat_dist, ndvi, method='bilinear') # or 'near' depending on your data
disturb_dist_aligned <- terra::resample(disturb_dist, ndvi, method='bilinear') # or 'near' depending on your data

##Log transform the distance rasters
constant = 0.1
habitat_dist_aligned_log_transformed <- log(habitat_dist_aligned+constant)
plot(habitat_dist_aligned_log_transformed)

constant = 0.1
disturb_dist_aligned_log_transformed <- log(disturb_dist_aligned+constant)
plot(disturb_dist_aligned_log_transformed)

#Split ndvi into habitats

# Create binary maps for each category
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

# set the origin to (0,0)
ext(resources) <- c(xmin - xmin, xmax - xmin, ymin - ymin, ymax - ymin)

#plot(resources)

hourly_coefs$true_hour<-hourly_coefs$hour

# Define the original hours in order
original_hours <- c(18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
# Create a sequence from 0 to 11.5 in steps of 0.5
new_sequence <- seq(1, 24, by=1)
# Create a mapping from original hours to new sequence
hour_mapping <- setNames(new_sequence, original_hours)
# Assuming hourly_coefs is your dataframe
# Map the original hour values to the new sequence
hourly_coefs$hour <- hour_mapping[as.character(hourly_coefs$hour)]
# Check the result
hourly_coefs

# naive predictions -------------------------------------------------------

naive_pred_stack <- c()
naive_pred_norm_stack <- c()

# Correct loop to iterate over specific half-hour intervals
for(i in seq(1, 24, by = 1)) {
  
  # Ensure 'i' matches an hour in 'hourly_coefs$hour'
  if(i %in% hourly_coefs$hour) {
    
    #i=0.5
    
    # dense_veg
    # Dense Vegetation Prediction
    dense_veg_lin <- resources[[1]]  # Assuming this is the binary map for dense vegetation
    dense_veg_pred <- dense_veg_lin * hourly_coefs$ndvidense_veg[which(hourly_coefs$hour == i)]
    
    # Grassland Prediction
    grassland_lin <- resources[[2]]  # Assuming this is the binary map for grassland
    grassland_pred <- grassland_lin * hourly_coefs$ndvigrassland[which(hourly_coefs$hour == i)]
    
    # Rocky Areas Prediction
    rocky_lin <- resources[[3]]  # Assuming this is the binary map for rocky areas
    rocky_pred <- rocky_lin * hourly_coefs$ndvirocky[which(hourly_coefs$hour == i)]
    
    # Other Disturbed Areas Prediction
    other_disturbed_lin <- resources[[4]]  # Assuming this is the binary map for other disturbed areas
    other_disturbed_pred <- other_disturbed_lin * hourly_coefs$ndviother_disturbed[which(hourly_coefs$hour == i)]
    
    # Mine Pits and Waste Dumps Prediction
    mine_pits_waste_dumps_lin <- resources[[5]]  # Assuming this is the binary map for mine pits and waste dumps
    mine_pits_waste_dumps_pred <- mine_pits_waste_dumps_lin * hourly_coefs$ndvimine_pit_waste_dump[which(hourly_coefs$hour == i)]
    
    # habitat distance
    habitat_lin <- resources[[6]] 
    habitat_pred <- habitat_lin * hourly_coefs$habitat_distance[which(hourly_coefs$hour == i)]
    
    # habitat distance
    disturb_lin <- resources[[7]] 
    disturb_pred <- disturb_lin * hourly_coefs$disturb_distance[which(hourly_coefs$hour == i)]
    
    # combine logit-scale predictions
    naive_pred_logit <- dense_veg_pred +
      grassland_pred +
      rocky_pred +
      other_disturbed_pred +
      mine_pits_waste_dumps_pred +
      habitat_pred +
      disturb_pred
    
    # add to list of raster layers
    naive_pred_stack <- c(naive_pred_stack, naive_pred_logit)
    
  }
}

# Combine stack into raster object
naive_predictions <- rast(naive_pred_stack)
#plot(naive_predictions)

# simulation function -----------------------------------------------------
# model_outputs<-readRDS("C:/Users/Mitch/OneDrive/Desktop/Simulation quoll paper/Step selection/Model outputs/model_twostep_2p_harms.RDS")
# summary(model_outputs)

# Define function to simulate SSF movement with memory component
simulate_ssf_memory <- function(n_steps, n_ch, coef, xy0, resc, memory_period, memory_delay, spatial_sd) { # hour, warmup, boundary
  # , exp_gamma
  
  tic()  # Start timer for benchmarking
  
  # Step length: gamma-distributed (shape and scale from coefficients)
  sl <- rgamma(n_steps * n_ch,
               shape = rep(coef$shape, each = n_ch),
               scale = rep(coef$scale, each = n_ch))
  
  # Turning angle concentration (kappa) and direction (mu)
  ta_coefs <- rep(coef$kappa, each = n_ch, length.out = n_steps)
  ta_coefs_mu <- ifelse(ta_coefs > 0, pi, 0)
  
  # Sample turning angles from von Mises distribution
  ta <- as.vector(mapply(Rfast::rvonmises, n = n_ch, m = ta_coefs_mu, k = abs(ta_coefs))) - pi
  
  # Step sequence and initial conditions
  steps <- rep(1:n_steps, each = n_ch)
  hour <- rep(seq(1,24,by=1), length.out = n_steps)
  x_0 <- xy0[[1]]
  y_0 <- xy0[[2]]
  angle_0 <- runif(1, min = -pi, max = pi)
  
  x <- rep(NA, n_steps)
  y <- rep(NA, n_steps)
  step_length <- rep(NA, n_steps)
  angle <- rep(NA, n_steps)
  bearing <- rep(NA, n_steps)
  hab_p <- rep(NA, n_steps)
  mem_p <- rep(NA, n_steps)
  
  x[1] <- x_0
  y[1] <- y_0
  step_length[1] <- 0
  angle[1] <- angle_0
  bearing[1] <- 0
  hab_p[1] <- 0
  mem_p[1] <- 0
  
  for (i in 2:n_steps) {
    
    x_prop <- (x[i - 1] + sl[steps == i] * cos(bearing[i - 1] + ta[steps == i]))
    y_prop <- (y[i - 1] + sl[steps == i] * sin(bearing[i - 1] + ta[steps == i]))
    sl_prop <- sl[steps == i]
    angle_prop <- ta[steps == i]
    bearing_prop <- atan2(y_prop - y[i - 1], x_prop - x[i - 1])
    
    # Determine which raster layer to use for naive predictions
    resc_index <- ifelse(i %% terra::nlyr(naive_predictions) == 0, 
                         terra::nlyr(naive_predictions), 
                         i %% terra::nlyr(naive_predictions))
    
    # Extract naive predictions at candidate locations
    p <- terra::extract(resc[[resc_index]], cbind(x_prop, y_prop))[,1]
    
    # Calculate memory density depending on trajectory length so far
    if(sum(!is.na(x)) < 2*memory_delay) {
      memory_density <- mapply(function(x_prop, y_prop) {
        recent_x <- x[1:(i-1)]
        recent_y <- y[1:(i-1)]
        diff_x <- x_prop - recent_x
        diff_y <- y_prop - recent_y
        return(logSumExp(dnorm(diff_x, mean = 0, sd = spatial_sd, log = TRUE) +
                           dnorm(diff_y, mean = 0, sd = spatial_sd, log = TRUE)))
      }, x_prop, y_prop)
      memory_p <- coef$memory[resc_index] * memory_density
      
    } else if(sum(!is.na(x)) <= memory_period) {
      memory_density <- mapply(function(x_prop, y_prop) {
        recent_x <- x[1:(i-memory_delay)]
        recent_y <- y[1:(i-memory_delay)]
        diff_x <- x_prop - recent_x
        diff_y <- y_prop - recent_y
        return(logSumExp(dnorm(diff_x, mean = 0, sd = spatial_sd, log = TRUE) +
                           dnorm(diff_y, mean = 0, sd = spatial_sd, log = TRUE)))
      }, x_prop, y_prop)
      memory_p <- coef$memory[resc_index] * memory_density
      
    } else if(sum(!is.na(x)) > memory_period) {
      memory_density <- mapply(function(x_prop, y_prop) {
        recent_x <- x[(i-memory_period):(i-memory_delay)]
        recent_y <- y[(i-memory_period):(i-memory_delay)]
        diff_x <- x_prop - recent_x
        diff_y <- y_prop - recent_y
        return(logSumExp(dnorm(diff_x, mean = 0, sd = spatial_sd, log = TRUE) +
                           dnorm(diff_y, mean = 0, sd = spatial_sd, log = TRUE)))
      }, x_prop, y_prop)
      memory_p <- coef$memory[resc_index] * memory_density
      
    } else {
      print("Out of specified memory conditions - check memory arguments")
    }
    
    p[is.na(p)] <- -200
    memory_p[is.na(memory_p)] <- -200
    
    w <- sample(n_ch, 1, prob = exp(memory_p + p))
    
    x[i] <- x_prop[w]
    y[i] <- y_prop[w]
    step_length[i] <- sl_prop[w]
    angle[i] <- angle_prop[w]
    bearing[i] <- bearing_prop[w]
    hab_p[i] <- p[w]
    mem_p[i] <- memory_p[w]
    
  }
  (toc())  # End timing
  
  return(data_frame(x = x, 
                    y = y, 
                    hour = hour,
                    sl = step_length, 
                    angle = angle, 
                    bearing = bearing, 
                    hab_p = hab_p, 
                    mem_p = mem_p))
}

# preparing to run function -----------------------------------------------

# Load pre-selected uniform starting locations (centred to raster origin)
load(file = "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/quoll_starting_locations_individual_centric.RData")

#####################################################
######### To run the harmonic validation simulations use the four quoll starting locations below. You will also have to adjust the number of steps and indvs below.
# 
# load(file = "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/quoll_starting_locations_for_comparing_harmonic_pairs.RData")
#####################################################


# Visualise starting points relative to full landscape
plot(ndvi)
points(x = start_unif[,1] + xmin, y = start_unif[,2] + ymin, col = 'red')

# Store as matrix with x and y
start_quoll <- matrix(c(x = start_unif[,1], y = start_unif[,2]), ncol = 2, nrow = 100)

# Confirm format
str(start_quoll)

# Convert NDVI raster to data frame for plotting
ndvi_df <- as.data.frame(ndvi, xy = TRUE)
str(ndvi_df)
#str(sim_data_mem)

aes(x = x, y = y, fill = factor("Mining_landscape_original"))

#############################################################################################
#### The below code can be adjusted to change the number of indvs and steps #################
#############################################################################################

n_indvs <- 5 # number of individual 'animals' per individual starting location
n_steps <- 480 # number of steps in the trajectory - 4880 is 6 months in hours + the warm-up period which should be discarded
n_ch <- 50 # number of steps to be selected from
coef <- hourly_coefs
resc <- naive_predictions
# most of the memory has declined to 0 by 500 locations
memory_period <- 100000
memory_delay <- 24
#spatial_sd <- spatial_sd_ids
#exp_gamma <- exp_gamma_ids

###########
#CODE ADDED
# this one is still on the original scale, so the locations won't 
# plot here unless we add the xmin and ymin values back on
# plot(ndvi)
# # Add the start_quoll points to the plot
# points(start_quoll[, 1], start_quoll[, 2], col = 'red')

# which would be like:
plot(ndvi)
# Add the start_quoll points to the plot
points(start_quoll[, 1] + xmin, start_quoll[, 2] + ymin, col = 'red')
###########

##Make sure to rename the memory column in excel

# Main simulation loop across all individuals
for(k in quoll_ids) {
  
  # List to store simulations for this quoll
  stps_mem_list <- vector(mode = "list", length = n_indvs)
  
  tic()  # Start timer for individual
  
  for(j in 1:n_indvs) {
    
    #j=1
    
    # Run SSF simulation with memory
    stps_mem_list[[j]] <- simulate_ssf_memory(n_steps = n_steps,
                                              n_ch = n_ch,
                                              coef = hourly_coefs,
                                              xy0 = start_quoll[k,],
                                              resc = naive_predictions,
                                              memory_period = memory_period,
                                              memory_delay = memory_delay,
                                              spatial_sd = 689.4857 #spatial_sd_ids[k],
                                              # exp_gamma = exp_gamma_ids[k]
    )
  }
  
  toc()  # End timer
  
  # Format list of trajectories into one `amt`-style dataframe
  animals_mem <- data_frame(
    id = paste0("id", k, 1:n_indvs),
    track = map(stps_mem_list, ~ amt::track(
      x = .$x, 
      y = .$y, 
      t = ymd_hm("2021-06-01 00:00") + hours(1:nrow(.)), 
      angle = .$angle,
      bearing = .$bearing,
      hab_p = .$hab_p,
      mem_p = .$mem_p)))
  
  # Flatten nested tibble
  sim_data_mem <- unnest(animals_mem, cols = track)
  
  # Fix hour = 1 to hour = 24
  sim_data_mem$hour <- ifelse(hour(sim_data_mem$t_) == 1, 24, hour(sim_data_mem$t_))
  
  ###########
  #CODE ADDED
  # in this line I was adding the xmin and ymin back on, so I'll add this line again below
  #sim_data_mem_truexy <- sim_data_mem %>% mutate(x_ = x_ + xmin, y_ = y_ + ymin)
  
  plot(naive_predictions[[1]])
  points(sim_data_mem %>% pull(x_),
         sim_data_mem %>% pull(y_))
  
  sim_data_mem_truexy <- sim_data_mem %>% mutate(x_ = x_ + xmin, y_ = y_ + ymin)
  ###########
  
  # Optional ggplot code for visualising paths (commented out)
  # print(ggplot() +
  #         geom_raster(data = ndvi_df,
  #                     aes(x = x, y = y, fill = factor(Habitat_reclassified)),  # Ensure discrete treatment
  #                     alpha = 0.5) +
  #         geom_point(data = sim_data_mem,
  #                    aes(x = x_, y = y_, colour = id),
  #                    size = 0.75, alpha = 0.1) +
  #         geom_path(data = sim_data_mem,
  #                   aes(x = x_, y = y_, colour = id),
  #                   alpha = 0.75) +
  #         scale_fill_brewer(palette = "Greens") +  # Applies to Habitat_reclassified
  #         scale_color_viridis_d("Sim ID") +  # Applies to id
  #         coord_equal() +
  #         theme_classic() +
  #         theme(legend.position = "none"))
  
  # Define path and filename for saving outputs
  base_path <- "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Simulation trajectories/Animal-centric/Original landscape/"
  filename <- paste0( base_path,"1hrm", "id_",  k, "_EvM_mem", memory_period, "_", memory_delay, "_", n_indvs, "ind_", n_ch, "ch_", n_steps, "stp_", Sys.Date(), ".csv")
  print(filename)
  
  # Save simulation output to CSV
  write_csv(sim_data_mem_truexy, filename)
  
  print(paste("Completed:", k, "of", 100))
  
}

##################################################################
### Repeat for each landscape type with the respective layers ####
##################################################################
