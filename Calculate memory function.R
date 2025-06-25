library(tidyverse)
library(lubridate)
packages <- c("amt", "terra", "tictoc", "matrixStats", "beepr", "ks")
walk(packages, require, character.only = T)

rm(list = ls())
#Import quoll data
#
#These data include the random steps and covariates. The random steps are included in this dataset so that the spatiotemporal memory density can be estimated at every used and random step, but they are not used to estimate the bandwidth or temporal decay parameters.

#{r}

setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")

##Read in CSV. This includes the data from step selection from Cowan et al., 2024, and includes the final result from the below code, but the below code steps through the method to calculate the memory function.

quoll_data_rand_steps <- read.csv("Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")

#quoll_data_rand_steps$t1_ <- as.POSIXct(quoll_data_rand_steps$t1_, format="%d/%m/%Y %H:%M", tz = "Australia/Perth")
#quoll_data_rand_steps$t2_ <- as.POSIXct(quoll_data_rand_steps$t2_, format="%d/%m/%Y %H:%M", tz = "Australia/Perth")
# quoll_data_rand_steps$t1_ <- dmy_hm(quoll_data_rand_steps$t1_)
# quoll_data_rand_steps$t2_ <- dmy_hm(quoll_data_rand_steps$t2_)
# quoll_data_rand_steps$t2_rounded<-quoll_data_rand_steps$t2_
#quoll_data_rand_steps$t2_rounded <- as.POSIXct(quoll_data_rand_steps$t2_rounded, format="%d/%m/%Y %H:%M")
quoll_data_rand_steps$t1_ <- ymd_hms(quoll_data_rand_steps$t1_, tz = "Australia/Perth")
quoll_data_rand_steps$t2_ <- ymd_hms(quoll_data_rand_steps$t2_, tz = "Australia/Perth")
quoll_data_rand_steps$t2_rounded <- quoll_data_rand_steps$t2_

quoll_data_rand_steps$y<-quoll_data_rand_steps$case_


quoll_ids <- unique(quoll_data_rand_steps$id)

 ggplot(quoll_data_rand_steps %>% filter(case_ == 0)) +
   geom_point(aes(x = x2_, y = y2_, colour = hour_t2), alpha = 0.25) +
   coord_equal() +
   scale_colour_viridis_c() +
   theme_classic()
 
 ggplot() +
   geom_point(data = quoll_data_rand_steps %>% filter(case_ == 0), 
              aes(x = x2_, y = y2_), colour = "black", size = 0.25, alpha = 0.25) +
   geom_point(data = quoll_data_rand_steps %>% filter(y == 1), 
              aes(x = x2_, y = y2_), colour = "red", size = 0.25, alpha = 0.25) +
   coord_equal() +
   theme_classic()


#Estimate the KDE kernel bandwidth
str(quoll_data_rand_steps)
# convert to track object to use the amt 

# keep a separate object for the original data
# quoll_data_rand_steps <- quoll_data_rand_steps |> mk_track(id = id, x1_, y1_, t1_, order_by_ts = T, all_cols = T, crs = 32751) |> arrange(id)
quoll_data_pres_track <- quoll_data_rand_steps |> filter(y == 1) %>% mk_track(id = id, x1_, y1_, t1_, order_by_ts = T, all_cols = T, crs = 32751) |> arrange(id)

quoll_data_pres_track


buffer <- 5000
res <- 25

h_ref <- vector(mode = "list", length = length(quoll_ids))
h_ref_vector <- c()
h_pi <- vector(mode = "list", length = length(quoll_ids))
h_pi_vector <- c()

tic()

for(i in 1:length(quoll_ids)) {
  
  # i = 1
  
  # subset by quoll id
  quoll_data_id <- quoll_data_pres_track |> filter(y == 1 & id == quoll_ids[i])
  
  # create template raster
  # create extent of the raster
  xmin <- min(quoll_data_id$x2_) - buffer
  xmax <- max(quoll_data_id$x2_) + buffer
  ymin <- min(quoll_data_id$y2_) - buffer
  ymax <- max(quoll_data_id$y2_) + buffer
  template_raster <- rast(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, res = res, crs = crs("epsg:32751"))
  
  
  # reference bandwidth
  tic()
  h_ref[[i]] <- hr_kde_ref(quoll_data_id)
  toc()
  print(h_ref[[i]])
  h_ref_vector[i] <- h_ref[[i]][1]
  
  
  quoll_hr_kde_ref <- hr_kde(quoll_data_id, 
                               h = h_ref[[i]],
                               trast = template_raster, levels = c(0.5, 0.75, 0.95))
  plot(quoll_hr_kde_ref$ud)
  
  # manual bandwidth
  # quoll_hr_kde_man <- hr_kde(quoll_data_id, 
  #                              h = c(400,400),
  #                              trast = template_raster, levels = c(0.5, 0.75, 0.95))
  # plot(quoll_hr_kde_man$ud)
  
  # plug-in bandwidth
  tic()
  h_pi[[i]] <- hr_kde_pi(quoll_data_id)
  toc()
  print(h_pi[[i]])
  h_pi_vector[i] <- h_pi[[i]][1]

  quoll_hr_kde_pi <- hr_kde(quoll_data_id,
                           h = h_pi[[i]],
                           trast = template_raster, levels = c(0.5, 0.75, 0.95))
  plot(quoll_hr_kde_pi$ud)

  # b_lscv <- hr_kde_lscv(quoll_data_id, trast = template_raster, which_min = "local")

  # least-squares cross validation bandwidth
  # quoll_hr_kde_lscv <- hr_kde(quoll_data_id,
  #                          h = hr_kde_lscv(quoll_data_id),
  #                          trast = template_raster, levels = c(0.5, 0.75, 0.95))
  # plot(quoll_hr_kde_lscv$ud)
  
}

toc()

hist(h_ref_vector, breaks = 10)
mean(h_ref_vector)
mean_kde_sd <- mean(h_ref_vector)
mean_kde_sd

# Estimating temporal decay parameter -------------------------------------
#Estimating the SD of the mixture density normal distributions (i.e. KDE) over each location


# SUBSETTING BY THE NUMBER OF HOURS

density_space_time_all <- function(locations_x, 
                                   locations_y,
                                   locations_time,
                                   # the params argument contains the sd that is shared in the x and y direction, and the temporal decay gamma parameter
                                   # params,
                                   spatial_sd,
                                   gamma_param,
                                   # the memory period should cover the duration that the memory decays to nearly 0. 
                                   memory_period, 
                                   # number of points to exclude
                                   memory_delay) {
  
  # create object for log_likelihood
  log_likelihood <- 0
  
  for(i in 1:(length(locations_x)-(i+memory_period))) {
    
    # extract the current location and time
    location_current_x <- locations_x[i+memory_period]
    location_current_y <- locations_y[i+memory_period]
    location_current_time <- locations_time[i+memory_period]
    
    # calculate the start and end times for the memory period (from the start of the memory period until the memory delay begins)
    memory_start_time <- location_current_time - as.difftime(memory_period, units = "hours")
    memory_end_time <- location_current_time - as.difftime(memory_delay, units = "hours")
    
    # subset the x, y, and time vectors by the locations BETWEEN the memory period and memory delay
    memory_period_x <- locations_x[locations_time >= memory_start_time & locations_time <= memory_end_time]
    memory_period_y <- locations_y[locations_time >= memory_start_time & locations_time <= memory_end_time]
    memory_period_time <- locations_time[locations_time >= memory_start_time & locations_time <= memory_end_time]
    num_locs <- length(memory_period_x)
    
    # spatial mixture density component
    # this subsets the locations from the start to the end (until the 24 hour delay) of the memory period for x and y coords and time
    # difference between locations in the x direction
    diff_x <- location_current_x - memory_period_x
    # difference between locations in the y direction
    diff_y <- location_current_y - memory_period_y
    # difference between locations in time
    diff_time <- as.numeric(difftime(location_current_time, memory_period_time, units = "hours"))
    
    # as all the parameters are on the log scale, then we can simply add them together to get the probability density for a given SD parameter and exponential decay parameter
    # we use the normal density function, where the density is defined by the distance (in the x or y direction, separately) and the SD parameter
    # the temporal decay component is already o
    log_joint_density <- dnorm(diff_x, mean = 0, sd = spatial_sd, log = TRUE) + 
      dnorm(diff_y, mean = 0, sd = spatial_sd, log = TRUE) + 
      (-gamma_param*diff_time)
    
    # we subtract the sum of the log temporal decay component, which is equivalent to dividing by the sum of the exponential components to normalise
    log_likelihood <- log_likelihood + logSumExp(log_joint_density) - log(num_locs) - logSumExp(-gamma_param*diff_time)
    
    # print(logSumExp(log_joint_density) - log(num_locs) - logSumExp(-gamma_param*diff_time))
    
  }
  return(-log_likelihood)
}






#Putting into a function to add the probability density at all used and random points
# EXCLUDING THE TEMPORAL DECAY COMPONENT

spatial_density_function <- function(data_input,
                                     id_val,
                                     memory_period,
                                     memory_delay,
                                     spatial_sd #,
                                     # temporal_decay_gamma
                                     ) {
  
  # for testing the function
  # data_input <- quoll_data_rand_steps %>% filter(id == quoll_ids[i])
  
  id_all_locations <- data_input %>% filter(id == id_val)
  id_used_locations <- data_input %>% filter(y == 1 & id == id_val)
  location_density <- c()
  
  for(i in 1:nrow(id_all_locations)){
    
    location_x <- id_all_locations[i,]$x2_
    location_y <- id_all_locations[i,]$y2_
    location_time <- id_all_locations[i,]$t2_
    
    memory_start_time <- location_time - as.difftime(memory_period, units = "hours")
    memory_end_time <- location_time - as.difftime(memory_delay, units = "hours")
    
    memory_x <- id_used_locations[id_used_locations$t2_ >= memory_start_time & id_used_locations$t2_ <= memory_end_time, ]$x1_
    memory_y <- id_used_locations[id_used_locations$t2_ >= memory_start_time & id_used_locations$t2_ <= memory_end_time, ]$y1_
    # memory_time <- id_used_locations[id_used_locations$t2_ >= memory_start_time & id_used_locations$t2_ <= memory_end_time, ]$t1_
    num_locs <- length(memory_x)
    
    diff_x <- location_x - memory_x
    diff_y <- location_y - memory_y
    # diff_time <- as.numeric(difftime(location_time, memory_time, units = "hours"))
    
    log_joint_density <- dnorm(diff_x, mean = 0, sd = spatial_sd, log = TRUE) + 
      dnorm(diff_y, mean = 0, sd = spatial_sd, log = TRUE) # + 
      # (-temporal_decay_gamma*diff_time)
    
    # normalise by the number of locations
    location_density[i] <- logSumExp(log_joint_density) - log(num_locs) # - logSumExp(-temporal_decay_gamma*diff_time)
    
  }
  
  return(location_density)
  
}



#Looping over individuals
spatial_sd_ids <- c()
exp_gamma_ids <- c()

for(i in 1:length(quoll_ids)) {
  spatial_sd_ids[i] <- h_ref_vector[i]
  # exp_gamma_ids[i] <- space_time_KDE_param_list_mem1000[[i]]$par
}

#memory_params <- data.frame("id" = quoll_ids, "kde_sd" = spatial_sd_ids, "temporal_decay" = exp_gamma_ids)
#write_csv(memory_params, paste0("memory_params_KDE_exp_decay_", Sys.Date(), ".csv"))

quoll_id_memory_spacetime1000_list <- vector(mode = "list", length = length(quoll_ids))

for(i in 1:length(quoll_ids)) {
  # tic()
  tic()
  quoll_id_memory_spacetime1000_list[[i]] <- quoll_data_rand_steps %>% filter(id == quoll_ids[i]) %>% 
    
    # note the 'spatial_temporal' in the name of the new column
    # mutate(kde_spatial_temporal_memory_density_log = spatial_temporal_density_function(., id_val = quoll_ids[i],
    #                                                                                    memory_period = 1000,
    #                                                                                    memory_delay = memory_delay,
    #                                                                                    spatial_sd = spatial_sd_ids[i],
    #                                                                                    temporal_decay_gamma = exp_gamma_ids[i])) %>% 
    
    # note the 'spatial' only in the name of the new column
    mutate(kde_ref_spatial_memory_density_log = spatial_density_function(., id_val = quoll_ids[i],
                                                                     memory_period = 1000,
                                                                     memory_delay = 24,
                                                                     spatial_sd = h_ref_vector[i])) %>% 
    
    mutate(kde_pi_spatial_memory_density_log = spatial_density_function(., id_val = quoll_ids[i],
                                                                     memory_period = 1000,
                                                                     memory_delay = 24,
                                                                     spatial_sd = h_pi_vector[i]))
  toc()
  
  # write_csv(quoll_id_memory_spacetime1000_list[[i]], paste0("st_memory_dfs/quoll_", quoll_ids[i], "_popn_EvM_covs_ST_KDEmemory1000_", Sys.Date(), ".csv"))
  
  # toc()
}

quoll_data_rand_steps_spatiotemporal_memory1000 <- bind_rows(quoll_id_memory_spacetime1000_list)
head(quoll_data_rand_steps_spatiotemporal_memory1000)

##Write the final CSV

#write_csv(quoll_data_rand_steps_spatiotemporal_memory1000, file = paste0("quoll_popn_GvM_covs_both_KDEmemory1000_10rs_", Sys.Date(), ".csv"))

# investigating the memory density

for(i in 1:length(quoll_ids)) {
  
  # WITHOUT TEMPORAL DECAY

  print(ggplot() +
    geom_point(data = quoll_data_rand_steps_spatiotemporal_memory1000 %>% filter(id == quoll_ids[i] & y == 0),
             aes(x = x2_, y = y2_, color = kde_ref_spatial_memory_density_log), size = 1) +
    scale_color_viridis_c() +
    geom_point(data = quoll_data_rand_steps_spatiotemporal_memory1000 %>% filter(id == quoll_ids[i] & y == 1),
               aes(x = x2_, y = y2_), colour = "red", size = 1) +
    theme_bw() +
      theme(legend.position = "bottom"))
  
}

#beep(sound = 2)

hist(quoll_data_rand_steps_spatiotemporal_memory1000$kde_ref_spatial_memory_density_log, breaks = 100)
hist(exp(quoll_data_rand_steps_spatiotemporal_memory1000$kde_ref_spatial_memory_density_log), breaks = 100)
# the first locations will be NA - before the memory delay kicks in, which will be multiplied by the number of random steps
# sum(is.na(quoll_data_rand_steps_spatiotemporal_memory1000$spatiotemporal_memory_density_log))

# quoll_data_rand_steps_spatiotemporal_memory1000 |> ggplot() +
#   geom_point(aes(x = t_, y = spatiotemporal_memory_density_log, color = id)) 



# for(i in 1:length(spatial_sd_ids)) plot(1:memory_period, exp(-exp_gamma_ids[i]*(memory_period:1)), type = "l", ylim = c(0, 1))


#Plotting the optimised spatial sd parameter

range <- 3500

# for(i in 1:length(spatial_sd_ids)) plot(-range:range, dnorm(-range:range, mean = 0, sd = spatial_sd_ids[i]), type = "l")

spatial_sds <- data.frame(-range:range, sapply(spatial_sd_ids, function(sd) dnorm(-range:range, mean = 0, sd = sd)))
colnames(spatial_sds) <- c("x", quoll_ids)

spatial_sds_long <- pivot_longer(spatial_sds, cols = !1)

ggplot(spatial_sds_long |> filter(name != "2154")) +
  geom_line(aes(x = -x, y = value, colour = name), size = 1) +
  scale_x_continuous("Metres from location") +
  scale_y_continuous("Density") +
  scale_colour_viridis_d("ID") +
  ggtitle("Estimated spatial SD function") +
  theme_classic() +
  theme(legend.position = "none")

#ggsave(paste0("plots/memory_process/spatial_KDE_sd_byID_", Sys.Date(), ".png"), width=150, height=90, units="mm", dpi = 300)


