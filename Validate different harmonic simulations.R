###Creates summaries of simulations from models with 0-4 pairs of harmonics

setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")
rm(list = ls())

library(tidyverse)
packages <- c("amt", "lubridate", "terra", "tictoc", "beepr", "ks",
              "adehabitatLT", "adehabitatHR", "ggpubr", "patchwork","dplyr")
walk(packages, require, character.only = T)

#quoll_ids <- unique(all_ssf_breeding_true$id)
quoll_ids <- c(1:4)

# unscaled rasters
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/current mining")
ndvi <- rast("Current_habitat_classifications.tif")
habitat_dist <- rast("Current_habitat_distances.tif")
disturb_dist <- rast("Current_disturb_distances.tif")
# Resample habitat_dist and disturb_dist to match ndvi's extent and resolution
habitat_dist_aligned <- terra::resample(habitat_dist, ndvi, method='bilinear') # or 'near' depending on your data
disturb_dist_aligned <- terra::resample(disturb_dist, ndvi, method='bilinear') # or 'near' depending on your data
##Log transform the distance rasters
constant = 0.1
habitat_dist_aligned_log_transformed <- log(habitat_dist_aligned+constant)
plot(habitat_dist_aligned_log_transformed)
constant = 0.1
disturb_dist_aligned_log_transformed <- log(disturb_dist_aligned+constant)
plot(disturb_dist_aligned_log_transformed)

xmin <- ext(ndvi[[1]])[1]
xmax <- ext(ndvi[[1]])[2]
ymin <- ext(ndvi[[1]])[3]
ymax <- ext(ndvi[[1]])[4]


original_hours <- c(18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
# Create a sequence from 0 to 11.5 in steps of 0.5
sequence <- seq(1, 24, by=1)
# Create a mapping from original hours to new sequence
hour_mapping <- setNames(original_hours,sequence)


############################################
########### 0 pair validation ##############
############################################
# read in multiple csv files with similar filenames and bind them together
sim_data_full_list <- 
  list.files("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/",  pattern = "*.csv", full.names = T)

sim_data_all <- grep("0hrm", sim_data_full_list, value = T) %>% 
  map_dfr(read_csv)%>%
  mutate(id = as.factor(id))

sim_data_all %>% group_by(id) %>% slice(1)

# Read in simulation summaries
#sim_data_all <-
#  read_csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/hrm0_all.csv") %>%
#  mutate(id = as.factor(id))
###Only do this once
##### Only do this if it wasn't done in the simulation save #
#sim_data_all <- sim_data_all %>% mutate(x_ = x_ + xmin, y_ = y_ + ymin) #####Only do this if it wasn't done in the simulation save
#############################################################

sim_data_all <- sim_data_all %>% ##Changes 0's to 1's.
  mutate(hour = ifelse(hour == 24, 1, hour))
sim_data_all <- sim_data_all %>% ##Changes 0's to 1's.
  mutate(hour = ifelse(hour == 0, 24, hour))
sim_data_all$true_hour <- hour_mapping[as.character(sim_data_all$hour)]
#Update time back to 'real' time 6pm to 6am
sim_data_all <- sim_data_all %>%
  mutate(new_hour = sprintf("%02d:%02d:00", true_hour %/% 1, (true_hour * 60) %% 60))
sim_data_all <- sim_data_all %>%
  mutate(
    date_component = as.Date(t_), # Extract date component
    time_component = format(strptime(new_hour, format = "%H:%M:%S"), "%T"), # Format new_hour as time
    t_ = paste(date_component, time_component) # Combine date and new time
  ) 
sim_data_all <- sim_data_all %>%
  mutate(t_ = lubridate::with_tz(sim_data_all$t_, tzone = "Australia/Perth"))
# Now remove the temporary columns by directly unsetting them
sim_data_all$date_component <- NULL
sim_data_all$time_component <- NULL
sim_data_all$new_hour <- NULL
sim_data_all$true_hour <- NULL
sim_data_all$hour <- NULL


sim_data_all %>% group_by(id) %>% slice(1)

sim_data_all$traj<-sim_data_all$id
# create vector of unique ids for subsetting
sim_data_traj_ids <- unique(sim_data_all$traj)

# convert to track object for step lengths and turning angles etc
sim_data_all <- sim_data_all %>% mk_track(id = id, x_, y_, t_, all_cols = T, crs = 32751)
sim_data_all_nested <- sim_data_all %>% arrange(traj) %>% nest(data = -"traj")

plot(sim_data_all)

# extract covariate values at the end of the all the simulated steps
sim_data_all_nested_steps <- sim_data_all_nested %>%
  mutate(steps = map(data, function(x)
    x %>% steps(keep_cols = "end") %>%
      extract_covariates(ndvi,
                         where = "end") %>%
      mutate(ndvi = factor(Current_habitat_classifications, levels = 1:5, labels = c("dense_veg","grassland", "rocky", "other_disturbed","mine_pit_waste_dump"))) %>%
      extract_covariates(habitat_dist_aligned_log_transformed,
                         where = "end") %>%
      mutate(habitat_distance_end = layer) %>%
      extract_covariates(disturb_dist_aligned_log_transformed,
                         where = "end") %>%
      mutate(disturb_distance_end = layer)))

# sim_data_all_nested_steps$steps[[1]]
sim_data_all_steps <- sim_data_all_nested_steps %>%
  amt::select(traj, steps) %>%
  amt::unnest(cols = steps)

sim_data_all_track <- sim_data_all_steps %>%
  mk_track(id = id, x1_, y1_, t1_, order_by_ts = T, all_cols = T, crs = 32751) %>%
  arrange(traj)
sim_ids <- unique(sim_data_all_steps$id)

#Plot
ggplot() +
  geom_path(data = sim_data_all_steps,
            aes(x = x2_, y = y2_, colour = traj),
            alpha = 0.1) +
  scale_color_viridis_d("id") +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "none")


for(i in 1:4) { # length(sim_data_ids)
  print(ggplot() +
          geom_point(data = sim_data_all_steps %>%
                       filter(traj == sim_data_traj_ids[i]),
                     aes(x = x2_, y = y2_, colour = date(t2_)),
                     size = 0.5, alpha = 0.5) +
          geom_path(data = sim_data_all_steps %>%
                      filter(traj == sim_data_traj_ids[i]),
                    aes(x = x2_, y = y2_, colour = date(t2_)),
                    alpha = 0.5) +
          scale_color_viridis_c("Time", trans = "date") +
          coord_equal() +
          theme_bw())
}


#Now for observed data
all_ssf1 <-read.csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")
all_ssf_breeding <- all_ssf1 %>% 
  filter(season == "breeding") %>% 
  arrange(t1_)
quoll_data_all <- all_ssf_breeding %>% 
  filter(case_ == "1") %>% 
  arrange(t1_)

quoll_data_all$t1_ <- as.POSIXct(quoll_data_all$t1_, format="%Y-%m-%dT%H:%M:%SZ", tz="Australia/Perth")
quoll_data_all$t2_ <- as.POSIXct(quoll_data_all$t2_, format="%Y-%m-%dT%H:%M:%SZ", tz="Australia/Perth")

quoll_data_all <- quoll_data_all %>%
  mutate(t1_ = lubridate::with_tz(quoll_data_all$t1_, tzone = "Australia/Perth"),
         t2_ = lubridate::with_tz(quoll_data_all$t2_, tzone = "Australia/Perth"))

quoll_data_all <- quoll_data_all %>%
  mutate(id_num = as.numeric(factor(id)),
         step_id = step_id_,
         x1 = x1_, x2 = x2_,
         y1 = y1_, y2 = y2_,
         t1 = t1_,
         t1_rounded = round_date(t1_, "hour"),
         hour_t1 = hour(t1_rounded),
         t2 = t2_,
         t2_rounded = round_date(t2_, "hour"),
         hour_t2 = hour(t2_rounded),
         hour = ifelse(hour_t2 == 0, 24, hour_t2),
         yday = yday(t1_),
         year = year(t1_),
         month = month(t1_),
         sl = sl_,
         log_sl = log(sl_),
         ta = ta_,
         cos_ta = cos(ta_),
         spatiotemporal_memory_density_log = kde_ref_spatial_memory_density_log,
         spatiotemporal_memory_density = exp(kde_ref_spatial_memory_density_log))


# ensure that we have just the individuals that were used for model fitting
quoll_year_ids <- unique(quoll_data_all$id)
quoll_data_all <- quoll_data_all %>% filter(id %in% quoll_year_ids)
quoll_ids <- unique(quoll_data_all$id)
# check the initial location and time for each individual
quoll_data_all %>% group_by(id) %>% slice(1)

# make a dataframe of only the presence locations
quoll_data_all_pres <- quoll_data_all %>% filter(y == 1)
# convert to track object for step lengths and turning angles etc
quoll_data_all_pres <- quoll_data_all_pres %>%
  mk_track(id = id, x1_, y1_, t1_, order_by_ts = T, all_cols = T, crs = 32751) %>%
  arrange(id)
#Subset the quoll data to the period that simulations were generated for
quoll_data_nested <- quoll_data_all_pres %>% arrange(id) %>% nest(data = -"id")
#Plot timeline of GPS data
quoll_data_all_pres %>% ggplot(aes(x = t1, y = factor(id), colour = factor(id))) +
  geom_point(alpha = 0.1) +
  scale_y_discrete("quoll ID") +
  scale_x_datetime("Date") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = "none")

for(i in 1:length(unique(quoll_data_all_pres$id))) {
  print(ggplot() +
          geom_point(data = quoll_data_all_pres %>% filter(id == unique(quoll_data_all_pres$id)[i]),
                     aes(x = x2_, y = y2_, colour = date(t2_)),
                     size = 0.5, alpha = 0.5) +
          geom_path(data = quoll_data_all_pres %>% filter(id == unique(quoll_data_all_pres$id)[i]),
                    aes(x = x2_, y = y2_, colour = date(t2_)),
                    alpha = 0.5) +
          scale_color_viridis_c("Time", trans = "date") +
          coord_equal() +
          theme_bw())
}

quoll_data <- quoll_data_all_pres


# points for each individual quoll
# combine list elements into data frame
sim_data_btime <- sim_data_all_track
sim_data_traj_ids <- unique(sim_data_btime$traj)
sim_data_btime$hour<-hour(sim_data_btime$t_)
str(sim_data_btime)
str(quoll_ids)

ndvi_xy <- as.data.frame(ndvi, xy = TRUE)
ndvi_xy$ndvi_discrete <- as.factor(ndvi_xy$Current_habitat_classifications)
n_sims <- 4
buffer <- 2500

for(i in 1:length(quoll_ids)) {
  # Find initial location for quoll data
  quoll_id_initial_x <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    slice(1) %>%
    pull(x_)
  quoll_id_initial_y <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    slice(1) %>%
    pull(y_)
  quoll_id_initial_df <- data.frame("X" = quoll_id_initial_x, "y" = quoll_id_initial_y)
  
  # Set the extent of the plot for quoll data
  extent_quoll <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    summarise(min_x = min(x_), min_y = min(y_), max_x = max(x_), max_y = max(y_))
  
  # Since we want to plot all sim_data_btime ids, no need to filter by quoll_ids for sim_data_btime
  extent_sim <- sim_data_btime %>%
    summarise(min_x = min(x_), min_y = min(y_), max_x = max(x_), max_y = max(y_))
  
  # Combine extents and add a buffer
  combined_min_x <- min(extent_sim$min_x, extent_quoll$min_x) - buffer
  combined_max_x <- max(extent_sim$max_x, extent_quoll$max_x) + buffer
  combined_min_y <- min(extent_sim$min_y, extent_quoll$min_y) - buffer
  combined_max_y <- max(extent_sim$max_y, extent_quoll$max_y) + buffer
  
  # Plot
  # plot <-  ggplot() +
  #     geom_raster(data = ndvi_xy,
  #                 aes(x = x, y = y, fill = ndvi_discrete),
  #                 alpha = 0.5) +
  #     scale_fill_brewer("NDVI", palette = "Greens", guide = guide_legend(reverse = TRUE)) +
  #     geom_path(data = sim_data_btime,
  #               aes(x = x2_, y = y2_, colour = traj),
  #               alpha = 0.75,
  #               linewidth = 0.25) +
  #     geom_point(data = sim_data_btime,
  #                aes(x = x2_, y = y2_, colour = traj),
  #                alpha = 0.75,
  #                size = 0.01) +
  #     geom_path(data = quoll_data %>% filter(id == quoll_ids[i]),
  #               aes(x = x2_, y = y2_),
  #               colour = "red",
  #               alpha = 0.75,
  #               linewidth = 0.25) +
  #     geom_point(data = quoll_data %>% filter(id == quoll_ids[i]),
  #                aes(x = x2_, y = y2_),
  #                colour = "red",
  #                alpha = 0.75,
  #                size = 0.01) +
  #     geom_point(data = quoll_id_initial_df, aes(x = X, y = y),
  #                colour = "blue",
  #                alpha = 1,
  #                shape = 4) +
  #     scale_color_viridis_d(guide = "none") +
  #     scale_x_continuous("Easting (m)", limits = c(combined_min_x, combined_max_x)) +
  #     scale_y_continuous("Northing (m)", limits = c(combined_min_y, combined_max_y)) +
  #     ggtitle(paste("Quoll ID:", quoll_ids[i])) +
  #     coord_equal() +
  #     theme_classic() +
  #     theme(legend.position = "right")
  #   
  #   # It's a good idea to explicitly print the plot if you're running this in a script
  #   print(plot)
  # 
}

head(quoll_data)

library(dplyr)
library(tidyr)

quoll_hourly_habitat <- quoll_data %>% 
  group_by(hour, id, NDVI_class_end) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(hour, id) %>%
  mutate(proportion = count / sum(count)) %>%
  pivot_wider(names_from = NDVI_class_end, 
              values_from = proportion, 
              values_fill = list(proportion = 0)) %>%
  left_join(
    quoll_data %>% 
      group_by(hour, id) %>%
      summarise(
        step_length_mean = mean(sl_, na.rm = TRUE),
        step_length_median = median(sl_, na.rm = TRUE),
        step_length_sd = sd(sl_, na.rm = TRUE),
        habitat_dist_mean = mean(habitat_distance_end, na.rm = TRUE),
        habitat_dist_median = median(habitat_distance_end, na.rm = TRUE),
        habitat_dist_sd = sd(habitat_distance_end, na.rm = TRUE),
        disturb_dist_mean = mean(disturb_distance_end, na.rm = TRUE),
        disturb_dist_median = median(disturb_distance_end, na.rm = TRUE),
        disturb_dist_sd = sd(disturb_distance_end, na.rm = TRUE),
        .groups = "drop")
  )



quoll_hourly_habitat <- data.frame("data" = "quoll", quoll_hourly_habitat) %>%
  mutate(id = as.factor(id))
# write.csv(quoll_hourly_habitat,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/quoll_summaries_hourly_habitat_", Sys.Date(), ".csv"))

quoll_hourly_habitat_long <- quoll_hourly_habitat %>%
  pivot_longer(cols = !c(data, hour, id), names_to = "variable", values_to = "value")

##Simulated summary
sim_hourly_habitat <- sim_data_btime %>%
  group_by(hour, id, ndvi) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(hour, id) %>%
  mutate(proportion = count / sum(count)) %>%
  pivot_wider(names_from = ndvi, 
              values_from = proportion, 
              values_fill = list(proportion = 0)) %>%
  left_join(
    sim_data_btime %>% 
      group_by(hour, id) %>%
      summarise(
        step_length_mean = mean(sl_, na.rm = TRUE),
        step_length_median = median(sl_, na.rm = TRUE),
        step_length_sd = sd(sl_, na.rm = TRUE),
        habitat_dist_mean = mean(habitat_distance_end, na.rm = TRUE),
        habitat_dist_median = median(habitat_distance_end, na.rm = TRUE),
        habitat_dist_sd = sd(habitat_distance_end, na.rm = TRUE),
        disturb_dist_mean = mean(disturb_distance_end, na.rm = TRUE),
        disturb_dist_median = median(disturb_distance_end, na.rm = TRUE),
        disturb_dist_sd = sd(disturb_distance_end, na.rm = TRUE),
        .groups = "drop")
  )


sim_hourly_habitat <- data.frame("data" = "0p", sim_hourly_habitat) %>%
  rename(id = id) %>% mutate(id = as.factor(id))
# write.csv(sim_hourly_habitat,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_0p_memALL_summaries_hourly_habitat_", Sys.Date(), ".csv"))

sim_hourly_habitat_long <- sim_hourly_habitat %>%
  pivot_longer(cols = !c(data, hour, id), names_to = "variable", values_to = "value")

# combine the dataframe
hourly_habitat_long <- bind_rows(quoll_hourly_habitat_long, sim_hourly_habitat_long)

#Plot the quoll vs. sim
for(i in 1:length(unique(hourly_habitat_long$variable))) {
  print(ggplot(data = hourly_habitat_long %>%
                 filter(variable == unique(variable)[i]),
               aes(x = factor(hour), y = value, colour = data)) +
          geom_boxplot() +
          ggtitle(unique(hourly_habitat_long$variable)[i]) +
          theme_classic())
}



quoll_mean_sl <- mean(quoll_data$sl_)
quoll_data_ltraj <- amt::as_ltraj(quoll_data, id = quoll_data$id)
quoll_RT <- residenceTime(quoll_data_ltraj,
                          radius = quoll_mean_sl,
                          maxt = 12, units = "hours", addinfo = FALSE)


###Get proportions of ndvi data
# Rocky habitat proportion by ID
rocky_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    rocky_count = sum(NDVI_class_end == "rocky", na.rm = TRUE),
    total_count = n(),
    rocky_prop = rocky_count / total_count,
    .groups = 'drop'
  )
# Grassland habitat proportion by ID
grassland_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    grassland_count = sum(NDVI_class_end == "grassland", na.rm = TRUE),
    total_count = n(),
    grassland_prop = grassland_count / total_count,
    .groups = 'drop'
  )
# Dense vegetation habitat proportion by ID
dense_veg_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    dense_veg_count = sum(NDVI_class_end == "dense_veg", na.rm = TRUE),
    total_count = n(),
    dense_veg_prop = dense_veg_count / total_count,
    .groups = 'drop'
  )
# Mine pit, waste dump habitat proportion by ID
mine_pit_waste_dump_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    mine_pit_waste_dump_count = sum(NDVI_class_end == "mine_pit_waste_dump", na.rm = TRUE),
    total_count = n(),
    mine_pit_waste_dump_prop = mine_pit_waste_dump_count / total_count,
    .groups = 'drop'
  )
# Other disturbed habitat proportion by ID
other_disturbed_prop_by_id <- quoll_data %>%  group_by(id) %>%  summarise(other_disturbed_count = sum(NDVI_class_end == "other_disturbed", na.rm = TRUE), total_count = n(), other_disturbed_prop = other_disturbed_count / total_count, .groups = 'drop')

####Calculate summary statistics for quoll
#Here we’ve created a loop that contains the summary statistics. We subset each individual animal’s trajectory
#and then below each simulated individual’s trajectory and calculate values for all of the summary statistics.
buffer <- 10000
res <- 100
# setup empty objects to store the results
id <- c()
step_length_median <- c()
step_length_mean <- c()
step_length_mean <- c()
step_length_median <- c()
step_length_sd <- c()
rocky_prop<-c()
grassland_prop<-c()
dense_veg_prop<-c()
mine_pit_waste_dump_prop<-c()
other_disturbed_prop<-c()
habitat_dist_mean <- c()
habitat_dist_median <- c()
habitat_dist_sd <- c()
disturb_dist_mean <- c()
disturb_dist_median <- c()
disturb_dist_sd <- c()
gamma_shape <- c()
gamma_scale <- c()
vm_kappa <- c()
straightness <- c()
tot_dist <- c()
intensity_use <- c()
sinuosity <- c()
tac <- c()
#residence_time <- c() # residence time in hours
hr_area_50 <- c()
hr_area_75 <- c()
hr_area_95 <- c()

for(k in 1:length(quoll_ids)) {
  quoll_data_id <- quoll_data %>% filter(id == quoll_ids[k])
  xmin <- min(quoll_data_id$x2_) - buffer
  xmax <- max(quoll_data_id$x2_) + buffer
  ymin <- min(quoll_data_id$y2_) - buffer
  ymax <- max(quoll_data_id$y2_) + buffer
  template_raster <- rast(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax,
                          res = res, crs = crs("epsg:32751"))
  id[k] <- quoll_data_id$id[1]
  step_length_median[k] <- median(quoll_data_id$sl_)
  step_length_mean[k] <- mean(quoll_data_id$sl_)
  rocky_prop<-rocky_prop_by_id$rocky_prop
  grassland_prop<-grassland_prop_by_id$grassland_prop
  dense_veg_prop<-dense_veg_prop_by_id$dense_veg_prop
  mine_pit_waste_dump_prop<-mine_pit_waste_dump_prop_by_id$mine_pit_waste_dump_prop
  other_disturbed_prop<-other_disturbed_prop_by_id$other_disturbed_prop
  habitat_dist_mean<-mean(quoll_data_id$habitat_distance_end)
  habitat_dist_median<-median(quoll_data_id$habitat_distance_end)
  habitat_dist_sd<-sd(quoll_data_id$habitat_distance_end)
  disturb_dist_mean<-mean(quoll_data_id$disturb_distance_end)
  disturb_dist_median<-median(quoll_data_id$disturb_distance_end)
  disturb_dist_sd<-sd(quoll_data_id$disturb_distance_end)
  gamma_fit <- fit_distr(quoll_data_id$sl_, "gamma")
  gamma_shape[k] <- gamma_fit$params$shape
  gamma_scale[k] <- gamma_fit$params$scale
  vM_fit <- fit_distr(quoll_data_id$sl_, "vonmises")
  vm_kappa[k] <- vM_fit$params$kappa
  straightness[k] <- amt::straightness(quoll_data_id)
  tot_dist[k] <- amt::tot_dist(quoll_data_id)
  intensity_use[k] <- amt::intensity_use(quoll_data_id)
  sinuosity[k] <- amt::sinuosity(quoll_data_id)
  tac[k] <- amt::tac(quoll_data_id)
  #residence_time[k] <- mean(quoll_RT[[k]][,2], na.rm = TRUE)/60/60
  quoll_hr_kde <- hr_kde(quoll_data_id, trast = template_raster,
                         levels = c(0.5, 0.75, 0.95))
  quoll_hr_kde_area <- hr_area(quoll_hr_kde)
  hr_area_50[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.5),]$area/1e6
  hr_area_75[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.75),]$area/1e6
  hr_area_95[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.95),]$area/1e6
}

# create a data frame that has traj, id, and all the summaries
quoll_summary_df <- data.frame(
  traj = "obs",
  id = as.numeric(id),
  data = "obs",
  sim = "obs",
  step_length_median = step_length_median,
  step_length_mean = step_length_mean,
  rocky_prop=rocky_prop,
  grassland_prop=grassland_prop,
  dense_veg_prop=dense_veg_prop,
  mine_pit_waste_dump_prop=mine_pit_waste_dump_prop,
  other_disturbed_prop=other_disturbed_prop,
  habitat_dist_mean=habitat_dist_mean,
  habitat_dist_median=habitat_dist_median,
  habitat_dist_sd=habitat_dist_sd,
  disturb_dist_mean=disturb_dist_mean,
  disturb_dist_median=disturb_dist_median,
  disturb_dist_sd = disturb_dist_sd,
  gamma_shape = gamma_shape,
  gamma_scale = gamma_scale,
  vm_kappa = vm_kappa,
  straightness = straightness,
  tot_dist = tot_dist,
  intensity_use = intensity_use,
  sinuosity = sinuosity,
  tac = tac,
  #residence_time = residence_time,
  hr_area_50 = hr_area_50,
  hr_area_75 = hr_area_75,
  hr_area_95 = hr_area_95)

# write_csv(quoll_summary_df,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/quoll_summary_statistics_df_", Sys.Date(), ".csv"))


###Now for simulated data
sim_mean_sl <- mean(sim_data_btime$sl_)
sim_data_ltraj <- amt::as_ltraj(sim_data_btime, id = sim_data_btime$traj)
sim_RT <- residenceTime(sim_data_ltraj, radius = sim_mean_sl, maxt = 12,
                        units = "hours", addinfo = FALSE)


###Get proportions of ndvi data
# Rocky habitat proportion by ID
rocky_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    rocky_count = sum(ndvi == "rocky", na.rm = TRUE),
    total_count = n(),
    rocky_prop = rocky_count / total_count,
    .groups = 'drop'
  )
# Grassland habitat proportion by ID
grassland_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    grassland_count = sum(ndvi == "grassland", na.rm = TRUE),
    total_count = n(),
    grassland_prop = grassland_count / total_count,
    .groups = 'drop'
  )
# Dense vegetation habitat proportion by ID
dense_veg_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    dense_veg_count = sum(ndvi == "dense_veg", na.rm = TRUE),
    total_count = n(),
    dense_veg_prop = dense_veg_count / total_count,
    .groups = 'drop'
  )
# Mine pit, waste dump habitat proportion by ID
mine_pit_waste_dump_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    mine_pit_waste_dump_count = sum(ndvi == "mine_pit_waste_dump", na.rm = TRUE),
    total_count = n(),
    mine_pit_waste_dump_prop = mine_pit_waste_dump_count / total_count,
    .groups = 'drop'
  )
# Other disturbed habitat proportion by ID
other_disturbed_prop_by_id <- sim_data_btime %>%  group_by(id) %>%  summarise(other_disturbed_count = sum(ndvi == "other_disturbed", na.rm = TRUE), total_count = n(), other_disturbed_prop = other_disturbed_count / total_count, .groups = 'drop')

####Calculate summary statistics for quoll
#Here we’ve created a loop that contains the summary statistics. We subset each individual animal’s trajectory
#and then below each simulated individual’s trajectory and calculate values for all of the summary statistics.
buffer <- 10000
res <- 100
# setup empty objects to store the results
id <- c()
step_length_median <- c()
step_length_mean <- c()
step_length_mean <- c()
step_length_median <- c()
step_length_sd <- c()
rocky_prop<-c()
grassland_prop<-c()
dense_veg_prop<-c()
mine_pit_waste_dump_prop<-c()
other_disturbed_prop<-c()
habitat_dist_mean <- c()
habitat_dist_median <- c()
habitat_dist_sd <- c()
disturb_dist_mean <- c()
disturb_dist_median <- c()
disturb_dist_sd <- c()
gamma_shape <- c()
gamma_scale <- c()
vm_kappa <- c()
straightness <- c()
tot_dist <- c()
intensity_use <- c()
sinuosity <- c()
tac <- c()
#residence_time <- c() # residence time in hours
hr_area_50 <- c()
hr_area_75 <- c()
hr_area_95 <- c()

for(k in 1:length(sim_data_traj_ids)) {
  sim_data_traj <- sim_data_btime %>% filter(id == sim_data_traj_ids[k])
  xmin <- min(sim_data_traj$x2_) - buffer
  xmax <- max(sim_data_traj$x2_) + buffer
  ymin <- min(sim_data_traj$y2_) - buffer
  ymax <- max(sim_data_traj$y2_) + buffer
  template_raster <- rast(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax,
                          res = res, crs = crs("epsg:32751"))
  id[k] <- sim_data_traj$id[1]
  step_length_median[k] <- median(sim_data_traj$sl_)
  step_length_mean[k] <- mean(sim_data_traj$sl_)
  rocky_prop<-rocky_prop_by_id$rocky_prop
  grassland_prop<-grassland_prop_by_id$grassland_prop
  dense_veg_prop<-dense_veg_prop_by_id$dense_veg_prop
  mine_pit_waste_dump_prop<-mine_pit_waste_dump_prop_by_id$mine_pit_waste_dump_prop
  other_disturbed_prop<-other_disturbed_prop_by_id$other_disturbed_prop
  habitat_dist_mean<-mean(sim_data_traj$habitat_distance_end)
  habitat_dist_median<-median(sim_data_traj$habitat_distance_end)
  habitat_dist_sd<-sd(sim_data_traj$habitat_distance_end)
  disturb_dist_mean<-mean(sim_data_traj$disturb_distance_end)
  disturb_dist_median<-median(sim_data_traj$disturb_distance_end)
  disturb_dist_sd<-sd(sim_data_traj$disturb_distance_end)
  gamma_fit <- fit_distr(sim_data_traj$sl_, "gamma")
  gamma_shape[k] <- gamma_fit$params$shape
  gamma_scale[k] <- gamma_fit$params$scale
  vM_fit <- fit_distr(sim_data_traj$sl_, "vonmises")
  vm_kappa[k] <- vM_fit$params$kappa
  straightness[k] <- amt::straightness(sim_data_traj)
  tot_dist[k] <- amt::tot_dist(sim_data_traj)
  intensity_use[k] <- amt::intensity_use(sim_data_traj)
  sinuosity[k] <- amt::sinuosity(sim_data_traj)
  tac[k] <- amt::tac(sim_data_traj)
  #residence_time[k] <- mean(quoll_RT[[k]][,2], na.rm = TRUE)/60/60
  sim_hr_kde <- hr_kde(sim_data_traj, trast = template_raster,
                       levels = c(0.5, 0.75, 0.95))
  sim_hr_kde_area <- hr_area(sim_hr_kde)
  hr_area_50[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.5),]$area/1e6
  hr_area_75[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.75),]$area/1e6
  hr_area_95[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.95),]$area/1e6
}

# create a data frame that has traj, id, and all the summaries
sim_summary_df <- data.frame(
  traj = "obs",
  id = as.numeric(id),
  data = "obs",
  sim = "obs",
  step_length_median = step_length_median,
  step_length_mean = step_length_mean,
  rocky_prop=rocky_prop,
  grassland_prop=grassland_prop,
  dense_veg_prop=dense_veg_prop,
  mine_pit_waste_dump_prop=mine_pit_waste_dump_prop,
  other_disturbed_prop=other_disturbed_prop,
  habitat_dist_mean=habitat_dist_mean,
  habitat_dist_median=habitat_dist_median,
  habitat_dist_sd=habitat_dist_sd,
  disturb_dist_mean=disturb_dist_mean,
  disturb_dist_median=disturb_dist_median,
  disturb_dist_sd = disturb_dist_sd,
  gamma_shape = gamma_shape,
  gamma_scale = gamma_scale,
  vm_kappa = vm_kappa,
  straightness = straightness,
  tot_dist = tot_dist,
  intensity_use = intensity_use,
  sinuosity = sinuosity,
  tac = tac,
  #residence_time = residence_time,
  hr_area_50 = hr_area_50,
  hr_area_75 = hr_area_75,
  hr_area_95 = hr_area_95)


# write_csv(sim_summary_df,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_0p_memALL_daily_summary_statistics_df_", Sys.Date(), ".csv"))




#


##########1 pair of harmonics

rm(list = ls())

library(tidyverse)
packages <- c("amt", "lubridate", "terra", "tictoc", "beepr", "ks",
              "adehabitatLT", "adehabitatHR", "ggpubr", "patchwork","dplyr")
walk(packages, require, character.only = T)

setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")

#quoll_ids <- unique(all_ssf_breeding_true$id)
quoll_ids <- c(1:4)

# unscaled rasters
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/current mining")
ndvi <- rast("Current_habitat_classifications.tif")
habitat_dist <- rast("Current_habitat_distances.tif")
disturb_dist <- rast("Current_disturb_distances.tif")
# Resample habitat_dist and disturb_dist to match ndvi's extent and resolution
habitat_dist_aligned <- terra::resample(habitat_dist, ndvi, method='bilinear') # or 'near' depending on your data
disturb_dist_aligned <- terra::resample(disturb_dist, ndvi, method='bilinear') # or 'near' depending on your data
##Log transform the distance rasters
constant = 0.1
habitat_dist_aligned_log_transformed <- log(habitat_dist_aligned+constant)
plot(habitat_dist_aligned_log_transformed)
constant = 0.1
disturb_dist_aligned_log_transformed <- log(disturb_dist_aligned+constant)
plot(disturb_dist_aligned_log_transformed)

xmin <- ext(ndvi[[1]])[1]
xmax <- ext(ndvi[[1]])[2]
ymin <- ext(ndvi[[1]])[3]
ymax <- ext(ndvi[[1]])[4]


original_hours <- c(18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
# Create a sequence from 0 to 11.5 in steps of 0.5
sequence <- seq(1, 24, by=1)
# Create a mapping from original hours to new sequence
hour_mapping <- setNames(original_hours,sequence)


############################################
########### 1 pair validation ##############
############################################
sim_data_full_list <- 
  list.files("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/",  pattern = "*.csv", full.names = T)

sim_data_all <- grep("1hrm", sim_data_full_list, value = T) %>% 
  map_dfr(read_csv)%>%
  mutate(id = as.factor(id))

sim_data_all %>% group_by(id) %>% slice(1)

###Only do this once
##### Only do this if it wasn't done in the simulation save #
#sim_data_all <- sim_data_all %>% mutate(x_ = x_ + xmin, y_ = y_ + ymin) #####Only do this if it wasn't done in the simulation save
#############################################################

sim_data_all <- sim_data_all %>% ##Changes 0's to 1's.
  mutate(hour = ifelse(hour == 24, 1, hour))
sim_data_all <- sim_data_all %>% ##Changes 0's to 1's.
  mutate(hour = ifelse(hour == 0, 24, hour))
sim_data_all$true_hour <- hour_mapping[as.character(sim_data_all$hour)]
#Update time back to 'real' time 6pm to 6am
sim_data_all <- sim_data_all %>%
  mutate(new_hour = sprintf("%02d:%02d:00", true_hour %/% 1, (true_hour * 60) %% 60))
sim_data_all <- sim_data_all %>%
  mutate(
    date_component = as.Date(t_), # Extract date component
    time_component = format(strptime(new_hour, format = "%H:%M:%S"), "%T"), # Format new_hour as time
    t_ = paste(date_component, time_component) # Combine date and new time
  ) 
sim_data_all <- sim_data_all %>%
  mutate(t_ = lubridate::with_tz(sim_data_all$t_, tzone = "Australia/Perth"))
# Now remove the temporary columns by directly unsetting them
sim_data_all$date_component <- NULL
sim_data_all$time_component <- NULL
sim_data_all$new_hour <- NULL
sim_data_all$true_hour <- NULL
sim_data_all$hour <- NULL


sim_data_all %>% group_by(id) %>% slice(1)

sim_data_all$traj<-sim_data_all$id
# create vector of unique ids for subsetting
sim_data_traj_ids <- unique(sim_data_all$traj)

# convert to track object for step lengths and turning angles etc
sim_data_all <- sim_data_all %>% mk_track(id = id, x_, y_, t_, all_cols = T, crs = 32751)
sim_data_all_nested <- sim_data_all %>% arrange(traj) %>% nest(data = -"traj")

plot(sim_data_all)

# extract covariate values at the end of the all the simulated steps
sim_data_all_nested_steps <- sim_data_all_nested %>%
  mutate(steps = map(data, function(x)
    x %>% steps(keep_cols = "end") %>%
      extract_covariates(ndvi,
                         where = "end") %>%
      mutate(ndvi = factor(Current_habitat_classifications, levels = 1:5, labels = c("dense_veg","grassland", "rocky", "other_disturbed","mine_pit_waste_dump"))) %>%
      extract_covariates(habitat_dist_aligned_log_transformed,
                         where = "end") %>%
      mutate(habitat_distance_end = layer) %>%
      extract_covariates(disturb_dist_aligned_log_transformed,
                         where = "end") %>%
      mutate(disturb_distance_end = layer)))

# sim_data_all_nested_steps$steps[[1]]
sim_data_all_steps <- sim_data_all_nested_steps %>%
  amt::select(traj, steps) %>%
  amt::unnest(cols = steps)

sim_data_all_track <- sim_data_all_steps %>%
  mk_track(id = id, x1_, y1_, t1_, order_by_ts = T, all_cols = T, crs = 32751) %>%
  arrange(traj)
sim_ids <- unique(sim_data_all_steps$id)

#Plot
ggplot() +
  geom_path(data = sim_data_all_steps,
            aes(x = x2_, y = y2_, colour = traj),
            alpha = 0.1) +
  scale_color_viridis_d("id") +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "none")


for(i in 1:4) { # length(sim_data_ids)
  print(ggplot() +
          geom_point(data = sim_data_all_steps %>%
                       filter(traj == sim_data_traj_ids[i]),
                     aes(x = x2_, y = y2_, colour = date(t2_)),
                     size = 0.5, alpha = 0.5) +
          geom_path(data = sim_data_all_steps %>%
                      filter(traj == sim_data_traj_ids[i]),
                    aes(x = x2_, y = y2_, colour = date(t2_)),
                    alpha = 0.5) +
          scale_color_viridis_c("Time", trans = "date") +
          coord_equal() +
          theme_bw())
}


#Now for observed data
all_ssf1 <-read.csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")
all_ssf_breeding <- all_ssf1 %>% 
  filter(season == "breeding") %>% 
  arrange(t1_)
quoll_data_all <- all_ssf_breeding %>% 
  filter(case_ == "1") %>% 
  arrange(t1_)

quoll_data_all$t1_ <- as.POSIXct(quoll_data_all$t1_, format="%Y-%m-%dT%H:%M:%SZ", tz="Australia/Perth")
quoll_data_all$t2_ <- as.POSIXct(quoll_data_all$t2_, format="%Y-%m-%dT%H:%M:%SZ", tz="Australia/Perth")

quoll_data_all <- quoll_data_all %>%
  mutate(t1_ = lubridate::with_tz(quoll_data_all$t1_, tzone = "Australia/Perth"),
         t2_ = lubridate::with_tz(quoll_data_all$t2_, tzone = "Australia/Perth"))

quoll_data_all <- quoll_data_all %>%
  mutate(id_num = as.numeric(factor(id)),
         step_id = step_id_,
         x1 = x1_, x2 = x2_,
         y1 = y1_, y2 = y2_,
         t1 = t1_,
         t1_rounded = round_date(t1_, "hour"),
         hour_t1 = hour(t1_rounded),
         t2 = t2_,
         t2_rounded = round_date(t2_, "hour"),
         hour_t2 = hour(t2_rounded),
         hour = ifelse(hour_t2 == 0, 24, hour_t2),
         yday = yday(t1_),
         year = year(t1_),
         month = month(t1_),
         sl = sl_,
         log_sl = log(sl_),
         ta = ta_,
         cos_ta = cos(ta_),
         spatiotemporal_memory_density_log = kde_ref_spatial_memory_density_log,
         spatiotemporal_memory_density = exp(kde_ref_spatial_memory_density_log))


# ensure that we have just the individuals that were used for model fitting
quoll_year_ids <- unique(quoll_data_all$id)
quoll_data_all <- quoll_data_all %>% filter(id %in% quoll_year_ids)
quoll_ids <- unique(quoll_data_all$id)
# check the initial location and time for each individual
quoll_data_all %>% group_by(id) %>% slice(1)

# make a dataframe of only the presence locations
quoll_data_all_pres <- quoll_data_all %>% filter(y == 1)
# convert to track object for step lengths and turning angles etc
quoll_data_all_pres <- quoll_data_all_pres %>%
  mk_track(id = id, x1_, y1_, t1_, order_by_ts = T, all_cols = T, crs = 32751) %>%
  arrange(id)
#Subset the quoll data to the period that simulations were generated for
quoll_data_nested <- quoll_data_all_pres %>% arrange(id) %>% nest(data = -"id")
#Plot timeline of GPS data
quoll_data_all_pres %>% ggplot(aes(x = t1, y = factor(id), colour = factor(id))) +
  geom_point(alpha = 0.1) +
  scale_y_discrete("quoll ID") +
  scale_x_datetime("Date") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = "none")

for(i in 1:length(unique(quoll_data_all_pres$id))) {
  print(ggplot() +
          geom_point(data = quoll_data_all_pres %>% filter(id == unique(quoll_data_all_pres$id)[i]),
                     aes(x = x2_, y = y2_, colour = date(t2_)),
                     size = 0.5, alpha = 0.5) +
          geom_path(data = quoll_data_all_pres %>% filter(id == unique(quoll_data_all_pres$id)[i]),
                    aes(x = x2_, y = y2_, colour = date(t2_)),
                    alpha = 0.5) +
          scale_color_viridis_c("Time", trans = "date") +
          coord_equal() +
          theme_bw())
}

quoll_data <- quoll_data_all_pres


# points for each individual quoll
# combine list elements into data frame
sim_data_btime <- sim_data_all_track
sim_data_traj_ids <- unique(sim_data_btime$traj)
sim_data_btime$hour<-hour(sim_data_btime$t_)
str(sim_data_btime)
str(quoll_ids)

ndvi_xy <- as.data.frame(ndvi, xy = TRUE)
ndvi_xy$ndvi_discrete <- as.factor(ndvi_xy$Current_habitat_classifications)
n_sims <- 4
buffer <- 2500

for(i in 1:length(quoll_ids)) {
  # Find initial location for quoll data
  quoll_id_initial_x <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    slice(1) %>%
    pull(x_)
  quoll_id_initial_y <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    slice(1) %>%
    pull(y_)
  quoll_id_initial_df <- data.frame("X" = quoll_id_initial_x, "y" = quoll_id_initial_y)
  
  # Set the extent of the plot for quoll data
  extent_quoll <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    summarise(min_x = min(x_), min_y = min(y_), max_x = max(x_), max_y = max(y_))
  
  # Since we want to plot all sim_data_btime ids, no need to filter by quoll_ids for sim_data_btime
  extent_sim <- sim_data_btime %>%
    summarise(min_x = min(x_), min_y = min(y_), max_x = max(x_), max_y = max(y_))
  
  # Combine extents and add a buffer
  combined_min_x <- min(extent_sim$min_x, extent_quoll$min_x) - buffer
  combined_max_x <- max(extent_sim$max_x, extent_quoll$max_x) + buffer
  combined_min_y <- min(extent_sim$min_y, extent_quoll$min_y) - buffer
  combined_max_y <- max(extent_sim$max_y, extent_quoll$max_y) + buffer
  
  # # Plot
  # plot <-  ggplot() +
  #   geom_raster(data = ndvi_xy,
  #               aes(x = x, y = y, fill = ndvi_discrete),
  #               alpha = 0.5) +
  #   scale_fill_brewer("NDVI", palette = "Greens", guide = guide_legend(reverse = TRUE)) +
  #   geom_path(data = sim_data_btime,
  #             aes(x = x2_, y = y2_, colour = traj),
  #             alpha = 0.75,
  #             linewidth = 0.25) +
  #   geom_point(data = sim_data_btime,
  #              aes(x = x2_, y = y2_, colour = traj),
  #              alpha = 0.75,
  #              size = 0.01) +
  #   geom_path(data = quoll_data %>% filter(id == quoll_ids[i]),
  #             aes(x = x2_, y = y2_),
  #             colour = "red",
  #             alpha = 0.75,
  #             linewidth = 0.25) +
  #   geom_point(data = quoll_data %>% filter(id == quoll_ids[i]),
  #              aes(x = x2_, y = y2_),
  #              colour = "red",
  #              alpha = 0.75,
  #              size = 0.01) +
  #   geom_point(data = quoll_id_initial_df, aes(x = X, y = y),
  #              colour = "blue",
  #              alpha = 1,
  #              shape = 4) +
  #   scale_color_viridis_d(guide = "none") +
  #   scale_x_continuous("Easting (m)", limits = c(combined_min_x, combined_max_x)) +
  #   scale_y_continuous("Northing (m)", limits = c(combined_min_y, combined_max_y)) +
  #   ggtitle(paste("Quoll ID:", quoll_ids[i])) +
  #   coord_equal() +
  #   theme_classic() +
  #   theme(legend.position = "right")
  # 
  # # It's a good idea to explicitly print the plot if you're running this in a script
  # print(plot)
}

head(quoll_data)

library(dplyr)
library(tidyr)

quoll_hourly_habitat <- quoll_data %>% 
  group_by(hour, id, NDVI_class_end) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(hour, id) %>%
  mutate(proportion = count / sum(count)) %>%
  pivot_wider(names_from = NDVI_class_end, 
              values_from = proportion, 
              values_fill = list(proportion = 0)) %>%
  left_join(
    quoll_data %>% 
      group_by(hour, id) %>%
      summarise(
        step_length_mean = mean(sl_, na.rm = TRUE),
        step_length_median = median(sl_, na.rm = TRUE),
        step_length_sd = sd(sl_, na.rm = TRUE),
        habitat_dist_mean = mean(habitat_distance_end, na.rm = TRUE),
        habitat_dist_median = median(habitat_distance_end, na.rm = TRUE),
        habitat_dist_sd = sd(habitat_distance_end, na.rm = TRUE),
        disturb_dist_mean = mean(disturb_distance_end, na.rm = TRUE),
        disturb_dist_median = median(disturb_distance_end, na.rm = TRUE),
        disturb_dist_sd = sd(disturb_distance_end, na.rm = TRUE),
        .groups = "drop")
  )



quoll_hourly_habitat <- data.frame("data" = "quoll", quoll_hourly_habitat) %>%
  mutate(id = as.factor(id))
# write.csv(quoll_hourly_habitat,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/quoll_summaries_hourly_habitat_", Sys.Date(), ".csv"))

quoll_hourly_habitat_long <- quoll_hourly_habitat %>%
  pivot_longer(cols = !c(data, hour, id), names_to = "variable", values_to = "value")

##Simulated summary
sim_hourly_habitat <- sim_data_btime %>%
  group_by(hour, id, ndvi) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(hour, id) %>%
  mutate(proportion = count / sum(count)) %>%
  pivot_wider(names_from = ndvi, 
              values_from = proportion, 
              values_fill = list(proportion = 0)) %>%
  left_join(
    sim_data_btime %>% 
      group_by(hour, id) %>%
      summarise(
        step_length_mean = mean(sl_, na.rm = TRUE),
        step_length_median = median(sl_, na.rm = TRUE),
        step_length_sd = sd(sl_, na.rm = TRUE),
        habitat_dist_mean = mean(habitat_distance_end, na.rm = TRUE),
        habitat_dist_median = median(habitat_distance_end, na.rm = TRUE),
        habitat_dist_sd = sd(habitat_distance_end, na.rm = TRUE),
        disturb_dist_mean = mean(disturb_distance_end, na.rm = TRUE),
        disturb_dist_median = median(disturb_distance_end, na.rm = TRUE),
        disturb_dist_sd = sd(disturb_distance_end, na.rm = TRUE),
        .groups = "drop")
  )


sim_hourly_habitat <- data.frame("data" = "1p", sim_hourly_habitat) %>%
  rename(id = id) %>% mutate(id = as.factor(id))
# write.csv(sim_hourly_habitat,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_1p_memALL_summaries_hourly_habitat_", Sys.Date(), ".csv"))

sim_hourly_habitat_long <- sim_hourly_habitat %>%
  pivot_longer(cols = !c(data, hour, id), names_to = "variable", values_to = "value")

# combine the dataframe
hourly_habitat_long <- bind_rows(quoll_hourly_habitat_long, sim_hourly_habitat_long)

#Plot the quoll vs. sim
for(i in 1:length(unique(hourly_habitat_long$variable))) {
  print(ggplot(data = hourly_habitat_long %>%
                 filter(variable == unique(variable)[i]),
               aes(x = factor(hour), y = value, colour = data)) +
          geom_boxplot() +
          ggtitle(unique(hourly_habitat_long$variable)[i]) +
          theme_classic())
}



quoll_mean_sl <- mean(quoll_data$sl_)
quoll_data_ltraj <- amt::as_ltraj(quoll_data, id = quoll_data$id)
quoll_RT <- residenceTime(quoll_data_ltraj,
                          radius = quoll_mean_sl,
                          maxt = 12, units = "hours", addinfo = FALSE)


###Get proportions of ndvi data
# Rocky habitat proportion by ID
rocky_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    rocky_count = sum(NDVI_class_end == "rocky", na.rm = TRUE),
    total_count = n(),
    rocky_prop = rocky_count / total_count,
    .groups = 'drop'
  )
# Grassland habitat proportion by ID
grassland_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    grassland_count = sum(NDVI_class_end == "grassland", na.rm = TRUE),
    total_count = n(),
    grassland_prop = grassland_count / total_count,
    .groups = 'drop'
  )
# Dense vegetation habitat proportion by ID
dense_veg_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    dense_veg_count = sum(NDVI_class_end == "dense_veg", na.rm = TRUE),
    total_count = n(),
    dense_veg_prop = dense_veg_count / total_count,
    .groups = 'drop'
  )
# Mine pit, waste dump habitat proportion by ID
mine_pit_waste_dump_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    mine_pit_waste_dump_count = sum(NDVI_class_end == "mine_pit_waste_dump", na.rm = TRUE),
    total_count = n(),
    mine_pit_waste_dump_prop = mine_pit_waste_dump_count / total_count,
    .groups = 'drop'
  )
# Other disturbed habitat proportion by ID
other_disturbed_prop_by_id <- quoll_data %>%  group_by(id) %>%  summarise(other_disturbed_count = sum(NDVI_class_end == "other_disturbed", na.rm = TRUE), total_count = n(), other_disturbed_prop = other_disturbed_count / total_count, .groups = 'drop')

####Calculate summary statistics for quoll
#Here we’ve created a loop that contains the summary statistics. We subset each individual animal’s trajectory
#and then below each simulated individual’s trajectory and calculate values for all of the summary statistics.
buffer <- 10000
res <- 100
# setup empty objects to store the results
id <- c()
step_length_median <- c()
step_length_mean <- c()
step_length_mean <- c()
step_length_median <- c()
step_length_sd <- c()
rocky_prop<-c()
grassland_prop<-c()
dense_veg_prop<-c()
mine_pit_waste_dump_prop<-c()
other_disturbed_prop<-c()
habitat_dist_mean <- c()
habitat_dist_median <- c()
habitat_dist_sd <- c()
disturb_dist_mean <- c()
disturb_dist_median <- c()
disturb_dist_sd <- c()
gamma_shape <- c()
gamma_scale <- c()
vm_kappa <- c()
straightness <- c()
tot_dist <- c()
intensity_use <- c()
sinuosity <- c()
tac <- c()
#residence_time <- c() # residence time in hours
hr_area_50 <- c()
hr_area_75 <- c()
hr_area_95 <- c()

for(k in 1:length(quoll_ids)) {
  quoll_data_id <- quoll_data %>% filter(id == quoll_ids[k])
  xmin <- min(quoll_data_id$x2_) - buffer
  xmax <- max(quoll_data_id$x2_) + buffer
  ymin <- min(quoll_data_id$y2_) - buffer
  ymax <- max(quoll_data_id$y2_) + buffer
  template_raster <- rast(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax,
                          res = res, crs = crs("epsg:32751"))
  id[k] <- quoll_data_id$id[1]
  step_length_median[k] <- median(quoll_data_id$sl_)
  step_length_mean[k] <- mean(quoll_data_id$sl_)
  rocky_prop<-rocky_prop_by_id$rocky_prop
  grassland_prop<-grassland_prop_by_id$grassland_prop
  dense_veg_prop<-dense_veg_prop_by_id$dense_veg_prop
  mine_pit_waste_dump_prop<-mine_pit_waste_dump_prop_by_id$mine_pit_waste_dump_prop
  other_disturbed_prop<-other_disturbed_prop_by_id$other_disturbed_prop
  habitat_dist_mean<-mean(quoll_data_id$habitat_distance_end)
  habitat_dist_median<-median(quoll_data_id$habitat_distance_end)
  habitat_dist_sd<-sd(quoll_data_id$habitat_distance_end)
  disturb_dist_mean<-mean(quoll_data_id$disturb_distance_end)
  disturb_dist_median<-median(quoll_data_id$disturb_distance_end)
  disturb_dist_sd<-sd(quoll_data_id$disturb_distance_end)
  gamma_fit <- fit_distr(quoll_data_id$sl_, "gamma")
  gamma_shape[k] <- gamma_fit$params$shape
  gamma_scale[k] <- gamma_fit$params$scale
  vM_fit <- fit_distr(quoll_data_id$sl_, "vonmises")
  vm_kappa[k] <- vM_fit$params$kappa
  straightness[k] <- amt::straightness(quoll_data_id)
  tot_dist[k] <- amt::tot_dist(quoll_data_id)
  intensity_use[k] <- amt::intensity_use(quoll_data_id)
  sinuosity[k] <- amt::sinuosity(quoll_data_id)
  tac[k] <- amt::tac(quoll_data_id)
  #residence_time[k] <- mean(quoll_RT[[k]][,2], na.rm = TRUE)/60/60
  quoll_hr_kde <- hr_kde(quoll_data_id, trast = template_raster,
                         levels = c(0.5, 0.75, 0.95))
  quoll_hr_kde_area <- hr_area(quoll_hr_kde)
  hr_area_50[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.5),]$area/1e6
  hr_area_75[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.75),]$area/1e6
  hr_area_95[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.95),]$area/1e6
}

# create a data frame that has traj, id, and all the summaries
quoll_summary_df <- data.frame(
  traj = "obs",
  id = as.numeric(id),
  data = "obs",
  sim = "obs",
  step_length_median = step_length_median,
  step_length_mean = step_length_mean,
  rocky_prop=rocky_prop,
  grassland_prop=grassland_prop,
  dense_veg_prop=dense_veg_prop,
  mine_pit_waste_dump_prop=mine_pit_waste_dump_prop,
  other_disturbed_prop=other_disturbed_prop,
  habitat_dist_mean=habitat_dist_mean,
  habitat_dist_median=habitat_dist_median,
  habitat_dist_sd=habitat_dist_sd,
  disturb_dist_mean=disturb_dist_mean,
  disturb_dist_median=disturb_dist_median,
  disturb_dist_sd = disturb_dist_sd,
  gamma_shape = gamma_shape,
  gamma_scale = gamma_scale,
  vm_kappa = vm_kappa,
  straightness = straightness,
  tot_dist = tot_dist,
  intensity_use = intensity_use,
  sinuosity = sinuosity,
  tac = tac,
  #residence_time = residence_time,
  hr_area_50 = hr_area_50,
  hr_area_75 = hr_area_75,
  hr_area_95 = hr_area_95)

# write_csv(quoll_summary_df,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/quoll_summary_statistics_df_", Sys.Date(), ".csv"))


###Now for simulated data
sim_mean_sl <- mean(sim_data_btime$sl_)
sim_data_ltraj <- amt::as_ltraj(sim_data_btime, id = sim_data_btime$traj)
sim_RT <- residenceTime(sim_data_ltraj, radius = sim_mean_sl, maxt = 12,
                        units = "hours", addinfo = FALSE)


###Get proportions of ndvi data
# Rocky habitat proportion by ID
rocky_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    rocky_count = sum(ndvi == "rocky", na.rm = TRUE),
    total_count = n(),
    rocky_prop = rocky_count / total_count,
    .groups = 'drop'
  )
# Grassland habitat proportion by ID
grassland_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    grassland_count = sum(ndvi == "grassland", na.rm = TRUE),
    total_count = n(),
    grassland_prop = grassland_count / total_count,
    .groups = 'drop'
  )
# Dense vegetation habitat proportion by ID
dense_veg_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    dense_veg_count = sum(ndvi == "dense_veg", na.rm = TRUE),
    total_count = n(),
    dense_veg_prop = dense_veg_count / total_count,
    .groups = 'drop'
  )
# Mine pit, waste dump habitat proportion by ID
mine_pit_waste_dump_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    mine_pit_waste_dump_count = sum(ndvi == "mine_pit_waste_dump", na.rm = TRUE),
    total_count = n(),
    mine_pit_waste_dump_prop = mine_pit_waste_dump_count / total_count,
    .groups = 'drop'
  )
# Other disturbed habitat proportion by ID
other_disturbed_prop_by_id <- sim_data_btime %>%  group_by(id) %>%  summarise(other_disturbed_count = sum(ndvi == "other_disturbed", na.rm = TRUE), total_count = n(), other_disturbed_prop = other_disturbed_count / total_count, .groups = 'drop')

####Calculate summary statistics for quoll
#Here we’ve created a loop that contains the summary statistics. We subset each individual animal’s trajectory
#and then below each simulated individual’s trajectory and calculate values for all of the summary statistics.
buffer <- 10000
res <- 100
# setup empty objects to store the results
id <- c()
step_length_median <- c()
step_length_mean <- c()
step_length_mean <- c()
step_length_median <- c()
step_length_sd <- c()
rocky_prop<-c()
grassland_prop<-c()
dense_veg_prop<-c()
mine_pit_waste_dump_prop<-c()
other_disturbed_prop<-c()
habitat_dist_mean <- c()
habitat_dist_median <- c()
habitat_dist_sd <- c()
disturb_dist_mean <- c()
disturb_dist_median <- c()
disturb_dist_sd <- c()
gamma_shape <- c()
gamma_scale <- c()
vm_kappa <- c()
straightness <- c()
tot_dist <- c()
intensity_use <- c()
sinuosity <- c()
tac <- c()
#residence_time <- c() # residence time in hours
hr_area_50 <- c()
hr_area_75 <- c()
hr_area_95 <- c()

for(k in 1:length(sim_data_traj_ids)) {
  sim_data_traj <- sim_data_btime %>% filter(id == sim_data_traj_ids[k])
  xmin <- min(sim_data_traj$x2_) - buffer
  xmax <- max(sim_data_traj$x2_) + buffer
  ymin <- min(sim_data_traj$y2_) - buffer
  ymax <- max(sim_data_traj$y2_) + buffer
  template_raster <- rast(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax,
                          res = res, crs = crs("epsg:32751"))
  id[k] <- sim_data_traj$id[1]
  step_length_median[k] <- median(sim_data_traj$sl_)
  step_length_mean[k] <- mean(sim_data_traj$sl_)
  rocky_prop<-rocky_prop_by_id$rocky_prop
  grassland_prop<-grassland_prop_by_id$grassland_prop
  dense_veg_prop<-dense_veg_prop_by_id$dense_veg_prop
  mine_pit_waste_dump_prop<-mine_pit_waste_dump_prop_by_id$mine_pit_waste_dump_prop
  other_disturbed_prop<-other_disturbed_prop_by_id$other_disturbed_prop
  habitat_dist_mean<-mean(sim_data_traj$habitat_distance_end)
  habitat_dist_median<-median(sim_data_traj$habitat_distance_end)
  habitat_dist_sd<-sd(sim_data_traj$habitat_distance_end)
  disturb_dist_mean<-mean(sim_data_traj$disturb_distance_end)
  disturb_dist_median<-median(sim_data_traj$disturb_distance_end)
  disturb_dist_sd<-sd(sim_data_traj$disturb_distance_end)
  gamma_fit <- fit_distr(sim_data_traj$sl_, "gamma")
  gamma_shape[k] <- gamma_fit$params$shape
  gamma_scale[k] <- gamma_fit$params$scale
  vM_fit <- fit_distr(sim_data_traj$sl_, "vonmises")
  vm_kappa[k] <- vM_fit$params$kappa
  straightness[k] <- amt::straightness(sim_data_traj)
  tot_dist[k] <- amt::tot_dist(sim_data_traj)
  intensity_use[k] <- amt::intensity_use(sim_data_traj)
  sinuosity[k] <- amt::sinuosity(sim_data_traj)
  tac[k] <- amt::tac(sim_data_traj)
  #residence_time[k] <- mean(quoll_RT[[k]][,2], na.rm = TRUE)/60/60
  sim_hr_kde <- hr_kde(sim_data_traj, trast = template_raster,
                       levels = c(0.5, 0.75, 0.95))
  sim_hr_kde_area <- hr_area(sim_hr_kde)
  hr_area_50[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.5),]$area/1e6
  hr_area_75[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.75),]$area/1e6
  hr_area_95[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.95),]$area/1e6
}

# create a data frame that has traj, id, and all the summaries
sim_summary_df <- data.frame(
  traj = "obs",
  id = as.numeric(id),
  data = "obs",
  sim = "obs",
  step_length_median = step_length_median,
  step_length_mean = step_length_mean,
  rocky_prop=rocky_prop,
  grassland_prop=grassland_prop,
  dense_veg_prop=dense_veg_prop,
  mine_pit_waste_dump_prop=mine_pit_waste_dump_prop,
  other_disturbed_prop=other_disturbed_prop,
  habitat_dist_mean=habitat_dist_mean,
  habitat_dist_median=habitat_dist_median,
  habitat_dist_sd=habitat_dist_sd,
  disturb_dist_mean=disturb_dist_mean,
  disturb_dist_median=disturb_dist_median,
  disturb_dist_sd = disturb_dist_sd,
  gamma_shape = gamma_shape,
  gamma_scale = gamma_scale,
  vm_kappa = vm_kappa,
  straightness = straightness,
  tot_dist = tot_dist,
  intensity_use = intensity_use,
  sinuosity = sinuosity,
  tac = tac,
  #residence_time = residence_time,
  hr_area_50 = hr_area_50,
  hr_area_75 = hr_area_75,
  hr_area_95 = hr_area_95)


# write_csv(sim_summary_df,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_1p_memALL_daily_summary_statistics_df_", Sys.Date(), ".csv"))



#





##########2 pair of harmonics

setwd("C:/Users/Mitch/OneDrive/Desktop/Simulation quoll paper/Step selection")
rm(list = ls())

library(tidyverse)
packages <- c("amt", "lubridate", "terra", "tictoc", "beepr", "ks",
              "adehabitatLT", "adehabitatHR", "ggpubr", "patchwork","dplyr")
walk(packages, require, character.only = T)

setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")

#quoll_ids <- unique(all_ssf_breeding_true$id)
quoll_ids <- c(1:4)

# unscaled rasters
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/current mining")
ndvi <- rast("Current_habitat_classifications.tif")
habitat_dist <- rast("Current_habitat_distances.tif")
disturb_dist <- rast("Current_disturb_distances.tif")
# Resample habitat_dist and disturb_dist to match ndvi's extent and resolution
habitat_dist_aligned <- terra::resample(habitat_dist, ndvi, method='bilinear') # or 'near' depending on your data
disturb_dist_aligned <- terra::resample(disturb_dist, ndvi, method='bilinear') # or 'near' depending on your data
##Log transform the distance rasters
constant = 0.1
habitat_dist_aligned_log_transformed <- log(habitat_dist_aligned+constant)
plot(habitat_dist_aligned_log_transformed)
constant = 0.1
disturb_dist_aligned_log_transformed <- log(disturb_dist_aligned+constant)
plot(disturb_dist_aligned_log_transformed)

xmin <- ext(ndvi[[1]])[1]
xmax <- ext(ndvi[[1]])[2]
ymin <- ext(ndvi[[1]])[3]
ymax <- ext(ndvi[[1]])[4]


original_hours <- c(18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
# Create a sequence from 0 to 11.5 in steps of 0.5
sequence <- seq(1, 24, by=1)
# Create a mapping from original hours to new sequence
hour_mapping <- setNames(original_hours,sequence)


############################################
########### 2 pair validation ##############
############################################
sim_data_full_list <- 
  list.files("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/",  pattern = "*.csv", full.names = T)

sim_data_all <- grep("2hrm", sim_data_full_list, value = T) %>% 
  map_dfr(read_csv)%>%
  mutate(id = as.factor(id))

sim_data_all %>% group_by(id) %>% slice(1)

###Only do this once
##### Only do this if it wasn't done in the simulation save #
#sim_data_all <- sim_data_all %>% mutate(x_ = x_ + xmin, y_ = y_ + ymin) #####Only do this if it wasn't done in the simulation save
#############################################################

sim_data_all <- sim_data_all %>% ##Changes 0's to 1's.
  mutate(hour = ifelse(hour == 24, 1, hour))
sim_data_all <- sim_data_all %>% ##Changes 0's to 1's.
  mutate(hour = ifelse(hour == 0, 24, hour))
sim_data_all$true_hour <- hour_mapping[as.character(sim_data_all$hour)]
#Update time back to 'real' time 6pm to 6am
sim_data_all <- sim_data_all %>%
  mutate(new_hour = sprintf("%02d:%02d:00", true_hour %/% 1, (true_hour * 60) %% 60))
sim_data_all <- sim_data_all %>%
  mutate(
    date_component = as.Date(t_), # Extract date component
    time_component = format(strptime(new_hour, format = "%H:%M:%S"), "%T"), # Format new_hour as time
    t_ = paste(date_component, time_component) # Combine date and new time
  ) 
sim_data_all <- sim_data_all %>%
  mutate(t_ = lubridate::with_tz(sim_data_all$t_, tzone = "Australia/Perth"))
# Now remove the temporary columns by directly unsetting them
sim_data_all$date_component <- NULL
sim_data_all$time_component <- NULL
sim_data_all$new_hour <- NULL
sim_data_all$true_hour <- NULL
sim_data_all$hour <- NULL


sim_data_all %>% group_by(id) %>% slice(1)

sim_data_all$traj<-sim_data_all$id
# create vector of unique ids for subsetting
sim_data_traj_ids <- unique(sim_data_all$traj)

# convert to track object for step lengths and turning angles etc
sim_data_all <- sim_data_all %>% mk_track(id = id, x_, y_, t_, all_cols = T, crs = 32751)
sim_data_all_nested <- sim_data_all %>% arrange(traj) %>% nest(data = -"traj")

plot(sim_data_all)

# extract covariate values at the end of the all the simulated steps
sim_data_all_nested_steps <- sim_data_all_nested %>%
  mutate(steps = map(data, function(x)
    x %>% steps(keep_cols = "end") %>%
      extract_covariates(ndvi,
                         where = "end") %>%
      mutate(ndvi = factor(Current_habitat_classifications, levels = 1:5, labels = c("dense_veg","grassland", "rocky", "other_disturbed","mine_pit_waste_dump"))) %>%
      extract_covariates(habitat_dist_aligned_log_transformed,
                         where = "end") %>%
      mutate(habitat_distance_end = layer) %>%
      extract_covariates(disturb_dist_aligned_log_transformed,
                         where = "end") %>%
      mutate(disturb_distance_end = layer)))

# sim_data_all_nested_steps$steps[[1]]
sim_data_all_steps <- sim_data_all_nested_steps %>%
  amt::select(traj, steps) %>%
  amt::unnest(cols = steps)

sim_data_all_track <- sim_data_all_steps %>%
  mk_track(id = id, x1_, y1_, t1_, order_by_ts = T, all_cols = T, crs = 32751) %>%
  arrange(traj)
sim_ids <- unique(sim_data_all_steps$id)

#Plot
ggplot() +
  geom_path(data = sim_data_all_steps,
            aes(x = x2_, y = y2_, colour = traj),
            alpha = 0.1) +
  scale_color_viridis_d("id") +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "none")


for(i in 1:4) { # length(sim_data_ids)
  print(ggplot() +
          geom_point(data = sim_data_all_steps %>%
                       filter(traj == sim_data_traj_ids[i]),
                     aes(x = x2_, y = y2_, colour = date(t2_)),
                     size = 0.5, alpha = 0.5) +
          geom_path(data = sim_data_all_steps %>%
                      filter(traj == sim_data_traj_ids[i]),
                    aes(x = x2_, y = y2_, colour = date(t2_)),
                    alpha = 0.5) +
          scale_color_viridis_c("Time", trans = "date") +
          coord_equal() +
          theme_bw())
}


#Now for observed data
all_ssf1 <-read.csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")
all_ssf_breeding <- all_ssf1 %>% 
  filter(season == "breeding") %>% 
  arrange(t1_)
quoll_data_all <- all_ssf_breeding %>% 
  filter(case_ == "1") %>% 
  arrange(t1_)

quoll_data_all$t1_ <- as.POSIXct(quoll_data_all$t1_, format="%Y-%m-%dT%H:%M:%SZ", tz="Australia/Perth")
quoll_data_all$t2_ <- as.POSIXct(quoll_data_all$t2_, format="%Y-%m-%dT%H:%M:%SZ", tz="Australia/Perth")

quoll_data_all <- quoll_data_all %>%
  mutate(t1_ = lubridate::with_tz(quoll_data_all$t1_, tzone = "Australia/Perth"),
         t2_ = lubridate::with_tz(quoll_data_all$t2_, tzone = "Australia/Perth"))

quoll_data_all <- quoll_data_all %>%
  mutate(id_num = as.numeric(factor(id)),
         step_id = step_id_,
         x1 = x1_, x2 = x2_,
         y1 = y1_, y2 = y2_,
         t1 = t1_,
         t1_rounded = round_date(t1_, "hour"),
         hour_t1 = hour(t1_rounded),
         t2 = t2_,
         t2_rounded = round_date(t2_, "hour"),
         hour_t2 = hour(t2_rounded),
         hour = ifelse(hour_t2 == 0, 24, hour_t2),
         yday = yday(t1_),
         year = year(t1_),
         month = month(t1_),
         sl = sl_,
         log_sl = log(sl_),
         ta = ta_,
         cos_ta = cos(ta_),
         spatiotemporal_memory_density_log = kde_ref_spatial_memory_density_log,
         spatiotemporal_memory_density = exp(kde_ref_spatial_memory_density_log))


# ensure that we have just the individuals that were used for model fitting
quoll_year_ids <- unique(quoll_data_all$id)
quoll_data_all <- quoll_data_all %>% filter(id %in% quoll_year_ids)
quoll_ids <- unique(quoll_data_all$id)
# check the initial location and time for each individual
quoll_data_all %>% group_by(id) %>% slice(1)

# make a dataframe of only the presence locations
quoll_data_all_pres <- quoll_data_all %>% filter(y == 1)
# convert to track object for step lengths and turning angles etc
quoll_data_all_pres <- quoll_data_all_pres %>%
  mk_track(id = id, x1_, y1_, t1_, order_by_ts = T, all_cols = T, crs = 32751) %>%
  arrange(id)
#Subset the quoll data to the period that simulations were generated for
quoll_data_nested <- quoll_data_all_pres %>% arrange(id) %>% nest(data = -"id")
#Plot timeline of GPS data
quoll_data_all_pres %>% ggplot(aes(x = t1, y = factor(id), colour = factor(id))) +
  geom_point(alpha = 0.1) +
  scale_y_discrete("quoll ID") +
  scale_x_datetime("Date") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = "none")

for(i in 1:length(unique(quoll_data_all_pres$id))) {
  print(ggplot() +
          geom_point(data = quoll_data_all_pres %>% filter(id == unique(quoll_data_all_pres$id)[i]),
                     aes(x = x2_, y = y2_, colour = date(t2_)),
                     size = 0.5, alpha = 0.5) +
          geom_path(data = quoll_data_all_pres %>% filter(id == unique(quoll_data_all_pres$id)[i]),
                    aes(x = x2_, y = y2_, colour = date(t2_)),
                    alpha = 0.5) +
          scale_color_viridis_c("Time", trans = "date") +
          coord_equal() +
          theme_bw())
}

quoll_data <- quoll_data_all_pres


# points for each individual quoll
# combine list elements into data frame
sim_data_btime <- sim_data_all_track
sim_data_traj_ids <- unique(sim_data_btime$traj)
sim_data_btime$hour<-hour(sim_data_btime$t_)
str(sim_data_btime)
str(quoll_ids)

ndvi_xy <- as.data.frame(ndvi, xy = TRUE)
ndvi_xy$ndvi_discrete <- as.factor(ndvi_xy$Current_habitat_classifications)
n_sims <- 4
buffer <- 2500

for(i in 1:length(quoll_ids)) {
  # Find initial location for quoll data
  quoll_id_initial_x <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    slice(1) %>%
    pull(x_)
  quoll_id_initial_y <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    slice(1) %>%
    pull(y_)
  quoll_id_initial_df <- data.frame("X" = quoll_id_initial_x, "y" = quoll_id_initial_y)
  
  # Set the extent of the plot for quoll data
  extent_quoll <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    summarise(min_x = min(x_), min_y = min(y_), max_x = max(x_), max_y = max(y_))
  
  # Since we want to plot all sim_data_btime ids, no need to filter by quoll_ids for sim_data_btime
  extent_sim <- sim_data_btime %>%
    summarise(min_x = min(x_), min_y = min(y_), max_x = max(x_), max_y = max(y_))
  
  # Combine extents and add a buffer
  combined_min_x <- min(extent_sim$min_x, extent_quoll$min_x) - buffer
  combined_max_x <- max(extent_sim$max_x, extent_quoll$max_x) + buffer
  combined_min_y <- min(extent_sim$min_y, extent_quoll$min_y) - buffer
  combined_max_y <- max(extent_sim$max_y, extent_quoll$max_y) + buffer
  
  # # Plot
  # plot <-  ggplot() +
  #   geom_raster(data = ndvi_xy,
  #               aes(x = x, y = y, fill = ndvi_discrete),
  #               alpha = 0.5) +
  #   scale_fill_brewer("NDVI", palette = "Greens", guide = guide_legend(reverse = TRUE)) +
  #   geom_path(data = sim_data_btime,
  #             aes(x = x2_, y = y2_, colour = traj),
  #             alpha = 0.75,
  #             linewidth = 0.25) +
  #   geom_point(data = sim_data_btime,
  #              aes(x = x2_, y = y2_, colour = traj),
  #              alpha = 0.75,
  #              size = 0.01) +
  #   geom_path(data = quoll_data %>% filter(id == quoll_ids[i]),
  #             aes(x = x2_, y = y2_),
  #             colour = "red",
  #             alpha = 0.75,
  #             linewidth = 0.25) +
  #   geom_point(data = quoll_data %>% filter(id == quoll_ids[i]),
  #              aes(x = x2_, y = y2_),
  #              colour = "red",
  #              alpha = 0.75,
  #              size = 0.01) +
  #   geom_point(data = quoll_id_initial_df, aes(x = X, y = y),
  #              colour = "blue",
  #              alpha = 1,
  #              shape = 4) +
  #   scale_color_viridis_d(guide = "none") +
  #   scale_x_continuous("Easting (m)", limits = c(combined_min_x, combined_max_x)) +
  #   scale_y_continuous("Northing (m)", limits = c(combined_min_y, combined_max_y)) +
  #   ggtitle(paste("Quoll ID:", quoll_ids[i])) +
  #   coord_equal() +
  #   theme_classic() +
  #   theme(legend.position = "right")
  # 
  # # It's a good idea to explicitly print the plot if you're running this in a script
  # print(plot)
}

head(quoll_data)

library(dplyr)
library(tidyr)

quoll_hourly_habitat <- quoll_data %>% 
  group_by(hour, id, NDVI_class_end) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(hour, id) %>%
  mutate(proportion = count / sum(count)) %>%
  pivot_wider(names_from = NDVI_class_end, 
              values_from = proportion, 
              values_fill = list(proportion = 0)) %>%
  left_join(
    quoll_data %>% 
      group_by(hour, id) %>%
      summarise(
        step_length_mean = mean(sl_, na.rm = TRUE),
        step_length_median = median(sl_, na.rm = TRUE),
        step_length_sd = sd(sl_, na.rm = TRUE),
        habitat_dist_mean = mean(habitat_distance_end, na.rm = TRUE),
        habitat_dist_median = median(habitat_distance_end, na.rm = TRUE),
        habitat_dist_sd = sd(habitat_distance_end, na.rm = TRUE),
        disturb_dist_mean = mean(disturb_distance_end, na.rm = TRUE),
        disturb_dist_median = median(disturb_distance_end, na.rm = TRUE),
        disturb_dist_sd = sd(disturb_distance_end, na.rm = TRUE),
        .groups = "drop")
  )



quoll_hourly_habitat <- data.frame("data" = "quoll", quoll_hourly_habitat) %>%
  mutate(id = as.factor(id))
# write.csv(quoll_hourly_habitat,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/quoll_summaries_hourly_habitat_", Sys.Date(), ".csv"))

quoll_hourly_habitat_long <- quoll_hourly_habitat %>%
  pivot_longer(cols = !c(data, hour, id), names_to = "variable", values_to = "value")

##Simulated summary
sim_hourly_habitat <- sim_data_btime %>%
  group_by(hour, id, ndvi) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(hour, id) %>%
  mutate(proportion = count / sum(count)) %>%
  pivot_wider(names_from = ndvi, 
              values_from = proportion, 
              values_fill = list(proportion = 0)) %>%
  left_join(
    sim_data_btime %>% 
      group_by(hour, id) %>%
      summarise(
        step_length_mean = mean(sl_, na.rm = TRUE),
        step_length_median = median(sl_, na.rm = TRUE),
        step_length_sd = sd(sl_, na.rm = TRUE),
        habitat_dist_mean = mean(habitat_distance_end, na.rm = TRUE),
        habitat_dist_median = median(habitat_distance_end, na.rm = TRUE),
        habitat_dist_sd = sd(habitat_distance_end, na.rm = TRUE),
        disturb_dist_mean = mean(disturb_distance_end, na.rm = TRUE),
        disturb_dist_median = median(disturb_distance_end, na.rm = TRUE),
        disturb_dist_sd = sd(disturb_distance_end, na.rm = TRUE),
        .groups = "drop")
  )


sim_hourly_habitat <- data.frame("data" = "2p", sim_hourly_habitat) %>%
  rename(id = id) %>% mutate(id = as.factor(id))
# write.csv(sim_hourly_habitat,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_2p_memALL_summaries_hourly_habitat_", Sys.Date(), ".csv"))

sim_hourly_habitat_long <- sim_hourly_habitat %>%
  pivot_longer(cols = !c(data, hour, id), names_to = "variable", values_to = "value")

# combine the dataframe
hourly_habitat_long <- bind_rows(quoll_hourly_habitat_long, sim_hourly_habitat_long)

#Plot the quoll vs. sim
for(i in 1:length(unique(hourly_habitat_long$variable))) {
  print(ggplot(data = hourly_habitat_long %>%
                 filter(variable == unique(variable)[i]),
               aes(x = factor(hour), y = value, colour = data)) +
          geom_boxplot() +
          ggtitle(unique(hourly_habitat_long$variable)[i]) +
          theme_classic())
}



quoll_mean_sl <- mean(quoll_data$sl_)
quoll_data_ltraj <- amt::as_ltraj(quoll_data, id = quoll_data$id)
quoll_RT <- residenceTime(quoll_data_ltraj,
                          radius = quoll_mean_sl,
                          maxt = 12, units = "hours", addinfo = FALSE)


###Get proportions of ndvi data
# Rocky habitat proportion by ID
rocky_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    rocky_count = sum(NDVI_class_end == "rocky", na.rm = TRUE),
    total_count = n(),
    rocky_prop = rocky_count / total_count,
    .groups = 'drop'
  )
# Grassland habitat proportion by ID
grassland_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    grassland_count = sum(NDVI_class_end == "grassland", na.rm = TRUE),
    total_count = n(),
    grassland_prop = grassland_count / total_count,
    .groups = 'drop'
  )
# Dense vegetation habitat proportion by ID
dense_veg_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    dense_veg_count = sum(NDVI_class_end == "dense_veg", na.rm = TRUE),
    total_count = n(),
    dense_veg_prop = dense_veg_count / total_count,
    .groups = 'drop'
  )
# Mine pit, waste dump habitat proportion by ID
mine_pit_waste_dump_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    mine_pit_waste_dump_count = sum(NDVI_class_end == "mine_pit_waste_dump", na.rm = TRUE),
    total_count = n(),
    mine_pit_waste_dump_prop = mine_pit_waste_dump_count / total_count,
    .groups = 'drop'
  )
# Other disturbed habitat proportion by ID
other_disturbed_prop_by_id <- quoll_data %>%  group_by(id) %>%  summarise(other_disturbed_count = sum(NDVI_class_end == "other_disturbed", na.rm = TRUE), total_count = n(), other_disturbed_prop = other_disturbed_count / total_count, .groups = 'drop')

####Calculate summary statistics for quoll
#Here we’ve created a loop that contains the summary statistics. We subset each individual animal’s trajectory
#and then below each simulated individual’s trajectory and calculate values for all of the summary statistics.
buffer <- 10000
res <- 100
# setup empty objects to store the results
id <- c()
step_length_median <- c()
step_length_mean <- c()
step_length_mean <- c()
step_length_median <- c()
step_length_sd <- c()
rocky_prop<-c()
grassland_prop<-c()
dense_veg_prop<-c()
mine_pit_waste_dump_prop<-c()
other_disturbed_prop<-c()
habitat_dist_mean <- c()
habitat_dist_median <- c()
habitat_dist_sd <- c()
disturb_dist_mean <- c()
disturb_dist_median <- c()
disturb_dist_sd <- c()
gamma_shape <- c()
gamma_scale <- c()
vm_kappa <- c()
straightness <- c()
tot_dist <- c()
intensity_use <- c()
sinuosity <- c()
tac <- c()
#residence_time <- c() # residence time in hours
hr_area_50 <- c()
hr_area_75 <- c()
hr_area_95 <- c()

for(k in 1:length(quoll_ids)) {
  quoll_data_id <- quoll_data %>% filter(id == quoll_ids[k])
  xmin <- min(quoll_data_id$x2_) - buffer
  xmax <- max(quoll_data_id$x2_) + buffer
  ymin <- min(quoll_data_id$y2_) - buffer
  ymax <- max(quoll_data_id$y2_) + buffer
  template_raster <- rast(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax,
                          res = res, crs = crs("epsg:32751"))
  id[k] <- quoll_data_id$id[1]
  step_length_median[k] <- median(quoll_data_id$sl_)
  step_length_mean[k] <- mean(quoll_data_id$sl_)
  rocky_prop<-rocky_prop_by_id$rocky_prop
  grassland_prop<-grassland_prop_by_id$grassland_prop
  dense_veg_prop<-dense_veg_prop_by_id$dense_veg_prop
  mine_pit_waste_dump_prop<-mine_pit_waste_dump_prop_by_id$mine_pit_waste_dump_prop
  other_disturbed_prop<-other_disturbed_prop_by_id$other_disturbed_prop
  habitat_dist_mean<-mean(quoll_data_id$habitat_distance_end)
  habitat_dist_median<-median(quoll_data_id$habitat_distance_end)
  habitat_dist_sd<-sd(quoll_data_id$habitat_distance_end)
  disturb_dist_mean<-mean(quoll_data_id$disturb_distance_end)
  disturb_dist_median<-median(quoll_data_id$disturb_distance_end)
  disturb_dist_sd<-sd(quoll_data_id$disturb_distance_end)
  gamma_fit <- fit_distr(quoll_data_id$sl_, "gamma")
  gamma_shape[k] <- gamma_fit$params$shape
  gamma_scale[k] <- gamma_fit$params$scale
  vM_fit <- fit_distr(quoll_data_id$sl_, "vonmises")
  vm_kappa[k] <- vM_fit$params$kappa
  straightness[k] <- amt::straightness(quoll_data_id)
  tot_dist[k] <- amt::tot_dist(quoll_data_id)
  intensity_use[k] <- amt::intensity_use(quoll_data_id)
  sinuosity[k] <- amt::sinuosity(quoll_data_id)
  tac[k] <- amt::tac(quoll_data_id)
  #residence_time[k] <- mean(quoll_RT[[k]][,2], na.rm = TRUE)/60/60
  quoll_hr_kde <- hr_kde(quoll_data_id, trast = template_raster,
                         levels = c(0.5, 0.75, 0.95))
  quoll_hr_kde_area <- hr_area(quoll_hr_kde)
  hr_area_50[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.5),]$area/1e6
  hr_area_75[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.75),]$area/1e6
  hr_area_95[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.95),]$area/1e6
}

# create a data frame that has traj, id, and all the summaries
quoll_summary_df <- data.frame(
  traj = "obs",
  id = as.numeric(id),
  data = "obs",
  sim = "obs",
  step_length_median = step_length_median,
  step_length_mean = step_length_mean,
  rocky_prop=rocky_prop,
  grassland_prop=grassland_prop,
  dense_veg_prop=dense_veg_prop,
  mine_pit_waste_dump_prop=mine_pit_waste_dump_prop,
  other_disturbed_prop=other_disturbed_prop,
  habitat_dist_mean=habitat_dist_mean,
  habitat_dist_median=habitat_dist_median,
  habitat_dist_sd=habitat_dist_sd,
  disturb_dist_mean=disturb_dist_mean,
  disturb_dist_median=disturb_dist_median,
  disturb_dist_sd = disturb_dist_sd,
  gamma_shape = gamma_shape,
  gamma_scale = gamma_scale,
  vm_kappa = vm_kappa,
  straightness = straightness,
  tot_dist = tot_dist,
  intensity_use = intensity_use,
  sinuosity = sinuosity,
  tac = tac,
  #residence_time = residence_time,
  hr_area_50 = hr_area_50,
  hr_area_75 = hr_area_75,
  hr_area_95 = hr_area_95)

# write_csv(quoll_summary_df,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/quoll_summary_statistics_df_", Sys.Date(), ".csv"))


###Now for simulated data
sim_mean_sl <- mean(sim_data_btime$sl_)
sim_data_ltraj <- amt::as_ltraj(sim_data_btime, id = sim_data_btime$traj)
sim_RT <- residenceTime(sim_data_ltraj, radius = sim_mean_sl, maxt = 12,
                        units = "hours", addinfo = FALSE)


###Get proportions of ndvi data
# Rocky habitat proportion by ID
rocky_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    rocky_count = sum(ndvi == "rocky", na.rm = TRUE),
    total_count = n(),
    rocky_prop = rocky_count / total_count,
    .groups = 'drop'
  )
# Grassland habitat proportion by ID
grassland_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    grassland_count = sum(ndvi == "grassland", na.rm = TRUE),
    total_count = n(),
    grassland_prop = grassland_count / total_count,
    .groups = 'drop'
  )
# Dense vegetation habitat proportion by ID
dense_veg_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    dense_veg_count = sum(ndvi == "dense_veg", na.rm = TRUE),
    total_count = n(),
    dense_veg_prop = dense_veg_count / total_count,
    .groups = 'drop'
  )
# Mine pit, waste dump habitat proportion by ID
mine_pit_waste_dump_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    mine_pit_waste_dump_count = sum(ndvi == "mine_pit_waste_dump", na.rm = TRUE),
    total_count = n(),
    mine_pit_waste_dump_prop = mine_pit_waste_dump_count / total_count,
    .groups = 'drop'
  )
# Other disturbed habitat proportion by ID
other_disturbed_prop_by_id <- sim_data_btime %>%  group_by(id) %>%  summarise(other_disturbed_count = sum(ndvi == "other_disturbed", na.rm = TRUE), total_count = n(), other_disturbed_prop = other_disturbed_count / total_count, .groups = 'drop')

####Calculate summary statistics for quoll
#Here we’ve created a loop that contains the summary statistics. We subset each individual animal’s trajectory
#and then below each simulated individual’s trajectory and calculate values for all of the summary statistics.
buffer <- 10000
res <- 100
# setup empty objects to store the results
id <- c()
step_length_median <- c()
step_length_mean <- c()
step_length_mean <- c()
step_length_median <- c()
step_length_sd <- c()
rocky_prop<-c()
grassland_prop<-c()
dense_veg_prop<-c()
mine_pit_waste_dump_prop<-c()
other_disturbed_prop<-c()
habitat_dist_mean <- c()
habitat_dist_median <- c()
habitat_dist_sd <- c()
disturb_dist_mean <- c()
disturb_dist_median <- c()
disturb_dist_sd <- c()
gamma_shape <- c()
gamma_scale <- c()
vm_kappa <- c()
straightness <- c()
tot_dist <- c()
intensity_use <- c()
sinuosity <- c()
tac <- c()
#residence_time <- c() # residence time in hours
hr_area_50 <- c()
hr_area_75 <- c()
hr_area_95 <- c()

for(k in 1:length(sim_data_traj_ids)) {
  sim_data_traj <- sim_data_btime %>% filter(id == sim_data_traj_ids[k])
  xmin <- min(sim_data_traj$x2_) - buffer
  xmax <- max(sim_data_traj$x2_) + buffer
  ymin <- min(sim_data_traj$y2_) - buffer
  ymax <- max(sim_data_traj$y2_) + buffer
  template_raster <- rast(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax,
                          res = res, crs = crs("epsg:32751"))
  id[k] <- sim_data_traj$id[1]
  step_length_median[k] <- median(sim_data_traj$sl_)
  step_length_mean[k] <- mean(sim_data_traj$sl_)
  rocky_prop<-rocky_prop_by_id$rocky_prop
  grassland_prop<-grassland_prop_by_id$grassland_prop
  dense_veg_prop<-dense_veg_prop_by_id$dense_veg_prop
  mine_pit_waste_dump_prop<-mine_pit_waste_dump_prop_by_id$mine_pit_waste_dump_prop
  other_disturbed_prop<-other_disturbed_prop_by_id$other_disturbed_prop
  habitat_dist_mean<-mean(sim_data_traj$habitat_distance_end)
  habitat_dist_median<-median(sim_data_traj$habitat_distance_end)
  habitat_dist_sd<-sd(sim_data_traj$habitat_distance_end)
  disturb_dist_mean<-mean(sim_data_traj$disturb_distance_end)
  disturb_dist_median<-median(sim_data_traj$disturb_distance_end)
  disturb_dist_sd<-sd(sim_data_traj$disturb_distance_end)
  gamma_fit <- fit_distr(sim_data_traj$sl_, "gamma")
  gamma_shape[k] <- gamma_fit$params$shape
  gamma_scale[k] <- gamma_fit$params$scale
  vM_fit <- fit_distr(sim_data_traj$sl_, "vonmises")
  vm_kappa[k] <- vM_fit$params$kappa
  straightness[k] <- amt::straightness(sim_data_traj)
  tot_dist[k] <- amt::tot_dist(sim_data_traj)
  intensity_use[k] <- amt::intensity_use(sim_data_traj)
  sinuosity[k] <- amt::sinuosity(sim_data_traj)
  tac[k] <- amt::tac(sim_data_traj)
  #residence_time[k] <- mean(quoll_RT[[k]][,2], na.rm = TRUE)/60/60
  sim_hr_kde <- hr_kde(sim_data_traj, trast = template_raster,
                       levels = c(0.5, 0.75, 0.95))
  sim_hr_kde_area <- hr_area(sim_hr_kde)
  hr_area_50[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.5),]$area/1e6
  hr_area_75[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.75),]$area/1e6
  hr_area_95[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.95),]$area/1e6
}

# create a data frame that has traj, id, and all the summaries
sim_summary_df <- data.frame(
  traj = "obs",
  id = as.numeric(id),
  data = "obs",
  sim = "obs",
  step_length_median = step_length_median,
  step_length_mean = step_length_mean,
  rocky_prop=rocky_prop,
  grassland_prop=grassland_prop,
  dense_veg_prop=dense_veg_prop,
  mine_pit_waste_dump_prop=mine_pit_waste_dump_prop,
  other_disturbed_prop=other_disturbed_prop,
  habitat_dist_mean=habitat_dist_mean,
  habitat_dist_median=habitat_dist_median,
  habitat_dist_sd=habitat_dist_sd,
  disturb_dist_mean=disturb_dist_mean,
  disturb_dist_median=disturb_dist_median,
  disturb_dist_sd = disturb_dist_sd,
  gamma_shape = gamma_shape,
  gamma_scale = gamma_scale,
  vm_kappa = vm_kappa,
  straightness = straightness,
  tot_dist = tot_dist,
  intensity_use = intensity_use,
  sinuosity = sinuosity,
  tac = tac,
  #residence_time = residence_time,
  hr_area_50 = hr_area_50,
  hr_area_75 = hr_area_75,
  hr_area_95 = hr_area_95)


# write_csv(sim_summary_df,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_2p_memALL_daily_summary_statistics_df_", Sys.Date(), ".csv"))



#













`

##########3 pair of harmonics
setwd("C:/Users/Mitch/OneDrive/Desktop/Simulation quoll paper/Step selection")
rm(list = ls())

library(tidyverse)
packages <- c("amt", "lubridate", "terra", "tictoc", "beepr", "ks",
              "adehabitatLT", "adehabitatHR", "ggpubr", "patchwork","dplyr")
walk(packages, require, character.only = T)

setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")

#quoll_ids <- unique(all_ssf_breeding_true$id)
quoll_ids <- c(1:4)

# unscaled rasters
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/current mining")
ndvi <- rast("Current_habitat_classifications.tif")
habitat_dist <- rast("Current_habitat_distances.tif")
disturb_dist <- rast("Current_disturb_distances.tif")
# Resample habitat_dist and disturb_dist to match ndvi's extent and resolution
habitat_dist_aligned <- terra::resample(habitat_dist, ndvi, method='bilinear') # or 'near' depending on your data
disturb_dist_aligned <- terra::resample(disturb_dist, ndvi, method='bilinear') # or 'near' depending on your data
##Log transform the distance rasters
constant = 0.1
habitat_dist_aligned_log_transformed <- log(habitat_dist_aligned+constant)
plot(habitat_dist_aligned_log_transformed)
constant = 0.1
disturb_dist_aligned_log_transformed <- log(disturb_dist_aligned+constant)
plot(disturb_dist_aligned_log_transformed)

xmin <- ext(ndvi[[1]])[1]
xmax <- ext(ndvi[[1]])[2]
ymin <- ext(ndvi[[1]])[3]
ymax <- ext(ndvi[[1]])[4]


original_hours <- c(18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
# Create a sequence from 0 to 11.5 in steps of 0.5
sequence <- seq(1, 24, by=1)
# Create a mapping from original hours to new sequence
hour_mapping <- setNames(original_hours,sequence)


############################################
########### 3 pair validation ##############
############################################
sim_data_full_list <- 
  list.files("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/",  pattern = "*.csv", full.names = T)

sim_data_all <- grep("3hrm", sim_data_full_list, value = T) %>% 
  map_dfr(read_csv)%>%
  mutate(id = as.factor(id))

sim_data_all %>% group_by(id) %>% slice(1)

###Only do this once
##### Only do this if it wasn't done in the simulation save #
#sim_data_all <- sim_data_all %>% mutate(x_ = x_ + xmin, y_ = y_ + ymin) #####Only do this if it wasn't done in the simulation save
#############################################################

sim_data_all <- sim_data_all %>% ##Changes 0's to 1's.
  mutate(hour = ifelse(hour == 24, 1, hour))
sim_data_all <- sim_data_all %>% ##Changes 0's to 1's.
  mutate(hour = ifelse(hour == 0, 24, hour))
sim_data_all$true_hour <- hour_mapping[as.character(sim_data_all$hour)]
#Update time back to 'real' time 6pm to 6am
sim_data_all <- sim_data_all %>%
  mutate(new_hour = sprintf("%02d:%02d:00", true_hour %/% 1, (true_hour * 60) %% 60))
sim_data_all <- sim_data_all %>%
  mutate(
    date_component = as.Date(t_), # Extract date component
    time_component = format(strptime(new_hour, format = "%H:%M:%S"), "%T"), # Format new_hour as time
    t_ = paste(date_component, time_component) # Combine date and new time
  ) 
sim_data_all <- sim_data_all %>%
  mutate(t_ = lubridate::with_tz(sim_data_all$t_, tzone = "Australia/Perth"))
# Now remove the temporary columns by directly unsetting them
sim_data_all$date_component <- NULL
sim_data_all$time_component <- NULL
sim_data_all$new_hour <- NULL
sim_data_all$true_hour <- NULL
sim_data_all$hour <- NULL


sim_data_all %>% group_by(id) %>% slice(1)

sim_data_all$traj<-sim_data_all$id
# create vector of unique ids for subsetting
sim_data_traj_ids <- unique(sim_data_all$traj)

# convert to track object for step lengths and turning angles etc
sim_data_all <- sim_data_all %>% mk_track(id = id, x_, y_, t_, all_cols = T, crs = 32751)
sim_data_all_nested <- sim_data_all %>% arrange(traj) %>% nest(data = -"traj")

plot(sim_data_all)

# extract covariate values at the end of the all the simulated steps
sim_data_all_nested_steps <- sim_data_all_nested %>%
  mutate(steps = map(data, function(x)
    x %>% steps(keep_cols = "end") %>%
      extract_covariates(ndvi,
                         where = "end") %>%
      mutate(ndvi = factor(Current_habitat_classifications, levels = 1:5, labels = c("dense_veg","grassland", "rocky", "other_disturbed","mine_pit_waste_dump"))) %>%
      extract_covariates(habitat_dist_aligned_log_transformed,
                         where = "end") %>%
      mutate(habitat_distance_end = layer) %>%
      extract_covariates(disturb_dist_aligned_log_transformed,
                         where = "end") %>%
      mutate(disturb_distance_end = layer)))

# sim_data_all_nested_steps$steps[[1]]
sim_data_all_steps <- sim_data_all_nested_steps %>%
  amt::select(traj, steps) %>%
  amt::unnest(cols = steps)

sim_data_all_track <- sim_data_all_steps %>%
  mk_track(id = id, x1_, y1_, t1_, order_by_ts = T, all_cols = T, crs = 32751) %>%
  arrange(traj)
sim_ids <- unique(sim_data_all_steps$id)

#Plot
ggplot() +
  geom_path(data = sim_data_all_steps,
            aes(x = x2_, y = y2_, colour = traj),
            alpha = 0.1) +
  scale_color_viridis_d("id") +
  coord_equal() +
  theme_bw() +
  theme(legend.position = "none")


for(i in 1:4) { # length(sim_data_ids)
  print(ggplot() +
          geom_point(data = sim_data_all_steps %>%
                       filter(traj == sim_data_traj_ids[i]),
                     aes(x = x2_, y = y2_, colour = date(t2_)),
                     size = 0.5, alpha = 0.5) +
          geom_path(data = sim_data_all_steps %>%
                      filter(traj == sim_data_traj_ids[i]),
                    aes(x = x2_, y = y2_, colour = date(t2_)),
                    alpha = 0.5) +
          scale_color_viridis_c("Time", trans = "date") +
          coord_equal() +
          theme_bw())
}


#Now for observed data
all_ssf1 <-read.csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")
all_ssf_breeding <- all_ssf1 %>% 
  filter(season == "breeding") %>% 
  arrange(t1_)
quoll_data_all <- all_ssf_breeding %>% 
  filter(case_ == "1") %>% 
  arrange(t1_)

quoll_data_all$disturb_distance_end <- log(quoll_data_all$disturb_distance_end+0.1)
quoll_data_all$habitat_distance_end <- log(quoll_data_all$habitat_distance_end+0.1)

quoll_data_all$t1_ <- as.POSIXct(quoll_data_all$t1_, format="%Y-%m-%dT%H:%M:%SZ", tz="Australia/Perth")
quoll_data_all$t2_ <- as.POSIXct(quoll_data_all$t2_, format="%Y-%m-%dT%H:%M:%SZ", tz="Australia/Perth")

quoll_data_all <- quoll_data_all %>%
  mutate(t1_ = lubridate::with_tz(quoll_data_all$t1_, tzone = "Australia/Perth"),
         t2_ = lubridate::with_tz(quoll_data_all$t2_, tzone = "Australia/Perth"))

quoll_data_all <- quoll_data_all %>%
  mutate(id_num = as.numeric(factor(id)),
         step_id = step_id_,
         x1 = x1_, x2 = x2_,
         y1 = y1_, y2 = y2_,
         t1 = t1_,
         t1_rounded = round_date(t1_, "hour"),
         hour_t1 = hour(t1_rounded),
         t2 = t2_,
         t2_rounded = round_date(t2_, "hour"),
         hour_t2 = hour(t2_rounded),
         hour = ifelse(hour_t2 == 0, 24, hour_t2),
         yday = yday(t1_),
         year = year(t1_),
         month = month(t1_),
         sl = sl_,
         log_sl = log(sl_),
         ta = ta_,
         cos_ta = cos(ta_),
         spatiotemporal_memory_density_log = kde_ref_spatial_memory_density_log,
         spatiotemporal_memory_density = exp(kde_ref_spatial_memory_density_log))


# ensure that we have just the individuals that were used for model fitting
quoll_year_ids <- unique(quoll_data_all$id)
quoll_data_all <- quoll_data_all %>% filter(id %in% quoll_year_ids)
quoll_ids <- unique(quoll_data_all$id)
# check the initial location and time for each individual
quoll_data_all %>% group_by(id) %>% slice(1)

# make a dataframe of only the presence locations
quoll_data_all_pres <- quoll_data_all %>% filter(y == 1)
# convert to track object for step lengths and turning angles etc
quoll_data_all_pres <- quoll_data_all_pres %>%
  mk_track(id = id, x1_, y1_, t1_, order_by_ts = T, all_cols = T, crs = 32751) %>%
  arrange(id)
#Subset the quoll data to the period that simulations were generated for
quoll_data_nested <- quoll_data_all_pres %>% arrange(id) %>% nest(data = -"id")
#Plot timeline of GPS data
quoll_data_all_pres %>% ggplot(aes(x = t1, y = factor(id), colour = factor(id))) +
  geom_point(alpha = 0.1) +
  scale_y_discrete("quoll ID") +
  scale_x_datetime("Date") +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = "none")

for(i in 1:length(unique(quoll_data_all_pres$id))) {
  print(ggplot() +
          geom_point(data = quoll_data_all_pres %>% filter(id == unique(quoll_data_all_pres$id)[i]),
                     aes(x = x2_, y = y2_, colour = date(t2_)),
                     size = 0.5, alpha = 0.5) +
          geom_path(data = quoll_data_all_pres %>% filter(id == unique(quoll_data_all_pres$id)[i]),
                    aes(x = x2_, y = y2_, colour = date(t2_)),
                    alpha = 0.5) +
          scale_color_viridis_c("Time", trans = "date") +
          coord_equal() +
          theme_bw())
}

quoll_data <- quoll_data_all_pres


# points for each individual quoll
# combine list elements into data frame
sim_data_btime <- sim_data_all_track
sim_data_traj_ids <- unique(sim_data_btime$traj)
sim_data_btime$hour<-hour(sim_data_btime$t_)
str(sim_data_btime)
str(quoll_ids)

ndvi_xy <- as.data.frame(ndvi, xy = TRUE)
ndvi_xy$ndvi_discrete <- as.factor(ndvi_xy$Current_habitat_classifications)
n_sims <- 4
buffer <- 2500

for(i in 1:length(quoll_ids)) {
  # Find initial location for quoll data
  quoll_id_initial_x <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    slice(1) %>%
    pull(x_)
  quoll_id_initial_y <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    slice(1) %>%
    pull(y_)
  quoll_id_initial_df <- data.frame("X" = quoll_id_initial_x, "y" = quoll_id_initial_y)
  
  # Set the extent of the plot for quoll data
  extent_quoll <- quoll_data %>%
    filter(id == quoll_ids[i]) %>%
    summarise(min_x = min(x_), min_y = min(y_), max_x = max(x_), max_y = max(y_))
  
  # Since we want to plot all sim_data_btime ids, no need to filter by quoll_ids for sim_data_btime
  extent_sim <- sim_data_btime %>%
    summarise(min_x = min(x_), min_y = min(y_), max_x = max(x_), max_y = max(y_))
  
  # Combine extents and add a buffer
  combined_min_x <- min(extent_sim$min_x, extent_quoll$min_x) - buffer
  combined_max_x <- max(extent_sim$max_x, extent_quoll$max_x) + buffer
  combined_min_y <- min(extent_sim$min_y, extent_quoll$min_y) - buffer
  combined_max_y <- max(extent_sim$max_y, extent_quoll$max_y) + buffer
  
  # # Plot
  # plot <-  ggplot() +
  #   geom_raster(data = ndvi_xy,
  #               aes(x = x, y = y, fill = ndvi_discrete),
  #               alpha = 0.5) +
  #   scale_fill_brewer("NDVI", palette = "Greens", guide = guide_legend(reverse = TRUE)) +
  #   geom_path(data = sim_data_btime,
  #             aes(x = x2_, y = y2_, colour = traj),
  #             alpha = 0.75,
  #             linewidth = 0.25) +
  #   geom_point(data = sim_data_btime,
  #              aes(x = x2_, y = y2_, colour = traj),
  #              alpha = 0.75,
  #              size = 0.01) +
  #   geom_path(data = quoll_data %>% filter(id == quoll_ids[i]),
  #             aes(x = x2_, y = y2_),
  #             colour = "red",
  #             alpha = 0.75,
  #             linewidth = 0.25) +
  #   geom_point(data = quoll_data %>% filter(id == quoll_ids[i]),
  #              aes(x = x2_, y = y2_),
  #              colour = "red",
  #              alpha = 0.75,
  #              size = 0.01) +
  #   geom_point(data = quoll_id_initial_df, aes(x = X, y = y),
  #              colour = "blue",
  #              alpha = 1,
  #              shape = 4) +
  #   scale_color_viridis_d(guide = "none") +
  #   scale_x_continuous("Easting (m)", limits = c(combined_min_x, combined_max_x)) +
  #   scale_y_continuous("Northing (m)", limits = c(combined_min_y, combined_max_y)) +
  #   ggtitle(paste("Quoll ID:", quoll_ids[i])) +
  #   coord_equal() +
  #   theme_classic() +
  #   theme(legend.position = "right")
  # 
  # # It's a good idea to explicitly print the plot if you're running this in a script
  # print(plot)
}

head(quoll_data)
library(dplyr)
library(tidyr)

quoll_hourly_habitat <- quoll_data %>% 
  group_by(hour, id, NDVI_class_end) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(hour, id) %>%
  mutate(proportion = count / sum(count)) %>%
  pivot_wider(names_from = NDVI_class_end, 
              values_from = proportion, 
              values_fill = list(proportion = 0)) %>%
  left_join(
    quoll_data %>% 
      group_by(hour, id) %>%
      summarise(
        step_length_mean = mean(sl_, na.rm = TRUE),
        step_length_median = median(sl_, na.rm = TRUE),
        step_length_sd = sd(sl_, na.rm = TRUE),
        habitat_dist_mean = mean(habitat_distance_end, na.rm = TRUE),
        habitat_dist_median = median(habitat_distance_end, na.rm = TRUE),
        habitat_dist_sd = sd(habitat_distance_end, na.rm = TRUE),
        disturb_dist_mean = mean(disturb_distance_end, na.rm = TRUE),
        disturb_dist_median = median(disturb_distance_end, na.rm = TRUE),
        disturb_dist_sd = sd(disturb_distance_end, na.rm = TRUE),
        .groups = "drop")
  )



quoll_hourly_habitat <- data.frame("data" = "quoll", quoll_hourly_habitat) %>%
  mutate(id = as.factor(id))
# write.csv(quoll_hourly_habitat,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/quoll_summaries_hourly_habitat_", Sys.Date(), ".csv"))

quoll_hourly_habitat_long <- quoll_hourly_habitat %>%
  pivot_longer(cols = !c(data, hour, id), names_to = "variable", values_to = "value")

##Simulated summary
sim_hourly_habitat <- sim_data_btime %>%
  group_by(hour, id, ndvi) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(hour, id) %>%
  mutate(proportion = count / sum(count)) %>%
  pivot_wider(names_from = ndvi, 
              values_from = proportion, 
              values_fill = list(proportion = 0)) %>%
  left_join(
    sim_data_btime %>% 
      group_by(hour, id) %>%
      summarise(
        step_length_mean = mean(sl_, na.rm = TRUE),
        step_length_median = median(sl_, na.rm = TRUE),
        step_length_sd = sd(sl_, na.rm = TRUE),
        habitat_dist_mean = mean(habitat_distance_end, na.rm = TRUE),
        habitat_dist_median = median(habitat_distance_end, na.rm = TRUE),
        habitat_dist_sd = sd(habitat_distance_end, na.rm = TRUE),
        disturb_dist_mean = mean(disturb_distance_end, na.rm = TRUE),
        disturb_dist_median = median(disturb_distance_end, na.rm = TRUE),
        disturb_dist_sd = sd(disturb_distance_end, na.rm = TRUE),
        .groups = "drop")
  )


sim_hourly_habitat <- data.frame("data" = "3p", sim_hourly_habitat) %>%
  rename(id = id) %>% mutate(id = as.factor(id))
# write.csv(sim_hourly_habitat,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_3p_memALL_summaries_hourly_habitat_", Sys.Date(), ".csv"))

sim_hourly_habitat_long <- sim_hourly_habitat %>%
  pivot_longer(cols = !c(data, hour, id), names_to = "variable", values_to = "value")

# combine the dataframe
hourly_habitat_long <- bind_rows(quoll_hourly_habitat_long, sim_hourly_habitat_long)

#Plot the quoll vs. sim
for(i in 1:length(unique(hourly_habitat_long$variable))) {
  print(ggplot(data = hourly_habitat_long %>%
                 filter(variable == unique(variable)[i]),
               aes(x = factor(hour), y = value, colour = data)) +
          geom_boxplot() +
          ggtitle(unique(hourly_habitat_long$variable)[i]) +
          theme_classic())
}



quoll_mean_sl <- mean(quoll_data$sl_)
quoll_data_ltraj <- amt::as_ltraj(quoll_data, id = quoll_data$id)
quoll_RT <- residenceTime(quoll_data_ltraj,
                          radius = quoll_mean_sl,
                          maxt = 12, units = "hours", addinfo = FALSE)


###Get proportions of ndvi data
# Rocky habitat proportion by ID
rocky_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    rocky_count = sum(NDVI_class_end == "rocky", na.rm = TRUE),
    total_count = n(),
    rocky_prop = rocky_count / total_count,
    .groups = 'drop'
  )
# Grassland habitat proportion by ID
grassland_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    grassland_count = sum(NDVI_class_end == "grassland", na.rm = TRUE),
    total_count = n(),
    grassland_prop = grassland_count / total_count,
    .groups = 'drop'
  )
# Dense vegetation habitat proportion by ID
dense_veg_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    dense_veg_count = sum(NDVI_class_end == "dense_veg", na.rm = TRUE),
    total_count = n(),
    dense_veg_prop = dense_veg_count / total_count,
    .groups = 'drop'
  )
# Mine pit, waste dump habitat proportion by ID
mine_pit_waste_dump_prop_by_id <- quoll_data %>%
  group_by(id) %>%
  summarise(
    mine_pit_waste_dump_count = sum(NDVI_class_end == "mine_pit_waste_dump", na.rm = TRUE),
    total_count = n(),
    mine_pit_waste_dump_prop = mine_pit_waste_dump_count / total_count,
    .groups = 'drop'
  )
# Other disturbed habitat proportion by ID
other_disturbed_prop_by_id <- quoll_data %>%  group_by(id) %>%  summarise(other_disturbed_count = sum(NDVI_class_end == "other_disturbed", na.rm = TRUE), total_count = n(), other_disturbed_prop = other_disturbed_count / total_count, .groups = 'drop')

####Calculate summary statistics for quoll
#Here we’ve created a loop that contains the summary statistics. We subset each individual animal’s trajectory
#and then below each simulated individual’s trajectory and calculate values for all of the summary statistics.
buffer <- 10000
res <- 100
# setup empty objects to store the results
id <- c()
step_length_median <- c()
step_length_mean <- c()
step_length_mean <- c()
step_length_median <- c()
step_length_sd <- c()
rocky_prop<-c()
grassland_prop<-c()
dense_veg_prop<-c()
mine_pit_waste_dump_prop<-c()
other_disturbed_prop<-c()
habitat_dist_mean <- c()
habitat_dist_median <- c()
habitat_dist_sd <- c()
disturb_dist_mean <- c()
disturb_dist_median <- c()
disturb_dist_sd <- c()
gamma_shape <- c()
gamma_scale <- c()
vm_kappa <- c()
straightness <- c()
tot_dist <- c()
intensity_use <- c()
sinuosity <- c()
tac <- c()
#residence_time <- c() # residence time in hours
hr_area_50 <- c()
hr_area_75 <- c()
hr_area_95 <- c()

for(k in 1:length(quoll_ids)) {
  quoll_data_id <- quoll_data %>% filter(id == quoll_ids[k])
  xmin <- min(quoll_data_id$x2_) - buffer
  xmax <- max(quoll_data_id$x2_) + buffer
  ymin <- min(quoll_data_id$y2_) - buffer
  ymax <- max(quoll_data_id$y2_) + buffer
  template_raster <- rast(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax,
                          res = res, crs = crs("epsg:32751"))
  id[k] <- quoll_data_id$id[1]
  step_length_median[k] <- median(quoll_data_id$sl_)
  step_length_mean[k] <- mean(quoll_data_id$sl_)
  rocky_prop<-rocky_prop_by_id$rocky_prop
  grassland_prop<-grassland_prop_by_id$grassland_prop
  dense_veg_prop<-dense_veg_prop_by_id$dense_veg_prop
  mine_pit_waste_dump_prop<-mine_pit_waste_dump_prop_by_id$mine_pit_waste_dump_prop
  other_disturbed_prop<-other_disturbed_prop_by_id$other_disturbed_prop
  habitat_dist_mean<-mean(quoll_data_id$habitat_distance_end)
  habitat_dist_median<-median(quoll_data_id$habitat_distance_end)
  habitat_dist_sd<-sd(quoll_data_id$habitat_distance_end)
  disturb_dist_mean<-mean(quoll_data_id$disturb_distance_end)
  disturb_dist_median<-median(quoll_data_id$disturb_distance_end)
  disturb_dist_sd<-sd(quoll_data_id$disturb_distance_end)
  gamma_fit <- fit_distr(quoll_data_id$sl_, "gamma")
  gamma_shape[k] <- gamma_fit$params$shape
  gamma_scale[k] <- gamma_fit$params$scale
  vM_fit <- fit_distr(quoll_data_id$sl_, "vonmises")
  vm_kappa[k] <- vM_fit$params$kappa
  straightness[k] <- amt::straightness(quoll_data_id)
  tot_dist[k] <- amt::tot_dist(quoll_data_id)
  intensity_use[k] <- amt::intensity_use(quoll_data_id)
  sinuosity[k] <- amt::sinuosity(quoll_data_id)
  tac[k] <- amt::tac(quoll_data_id)
  #residence_time[k] <- mean(quoll_RT[[k]][,2], na.rm = TRUE)/60/60
  quoll_hr_kde <- hr_kde(quoll_data_id, trast = template_raster,
                         levels = c(0.5, 0.75, 0.95))
  quoll_hr_kde_area <- hr_area(quoll_hr_kde)
  hr_area_50[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.5),]$area/1e6
  hr_area_75[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.75),]$area/1e6
  hr_area_95[k] = quoll_hr_kde_area[which(quoll_hr_kde_area$level == 0.95),]$area/1e6
}

# create a data frame that has traj, id, and all the summaries
quoll_summary_df <- data.frame(
  traj = "obs",
  id = as.numeric(id),
  data = "obs",
  sim = "obs",
  step_length_median = step_length_median,
  step_length_mean = step_length_mean,
  rocky_prop=rocky_prop,
  grassland_prop=grassland_prop,
  dense_veg_prop=dense_veg_prop,
  mine_pit_waste_dump_prop=mine_pit_waste_dump_prop,
  other_disturbed_prop=other_disturbed_prop,
  habitat_dist_mean=habitat_dist_mean,
  habitat_dist_median=habitat_dist_median,
  habitat_dist_sd=habitat_dist_sd,
  disturb_dist_mean=disturb_dist_mean,
  disturb_dist_median=disturb_dist_median,
  disturb_dist_sd = disturb_dist_sd,
  gamma_shape = gamma_shape,
  gamma_scale = gamma_scale,
  vm_kappa = vm_kappa,
  straightness = straightness,
  tot_dist = tot_dist,
  intensity_use = intensity_use,
  sinuosity = sinuosity,
  tac = tac,
  #residence_time = residence_time,
  hr_area_50 = hr_area_50,
  hr_area_75 = hr_area_75,
  hr_area_95 = hr_area_95)

# write_csv(quoll_summary_df,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/quoll_summary_statistics_df_", Sys.Date(), ".csv"))


###Now for simulated data
sim_mean_sl <- mean(sim_data_btime$sl_)
sim_data_ltraj <- amt::as_ltraj(sim_data_btime, id = sim_data_btime$traj)
sim_RT <- residenceTime(sim_data_ltraj, radius = sim_mean_sl, maxt = 12,
                        units = "hours", addinfo = FALSE)


###Get proportions of ndvi data
# Rocky habitat proportion by ID
rocky_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    rocky_count = sum(ndvi == "rocky", na.rm = TRUE),
    total_count = n(),
    rocky_prop = rocky_count / total_count,
    .groups = 'drop'
  )
# Grassland habitat proportion by ID
grassland_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    grassland_count = sum(ndvi == "grassland", na.rm = TRUE),
    total_count = n(),
    grassland_prop = grassland_count / total_count,
    .groups = 'drop'
  )
# Dense vegetation habitat proportion by ID
dense_veg_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    dense_veg_count = sum(ndvi == "dense_veg", na.rm = TRUE),
    total_count = n(),
    dense_veg_prop = dense_veg_count / total_count,
    .groups = 'drop'
  )
# Mine pit, waste dump habitat proportion by ID
mine_pit_waste_dump_prop_by_id <- sim_data_btime %>%
  group_by(id) %>%
  summarise(
    mine_pit_waste_dump_count = sum(ndvi == "mine_pit_waste_dump", na.rm = TRUE),
    total_count = n(),
    mine_pit_waste_dump_prop = mine_pit_waste_dump_count / total_count,
    .groups = 'drop'
  )
# Other disturbed habitat proportion by ID
other_disturbed_prop_by_id <- sim_data_btime %>%  group_by(id) %>%  summarise(other_disturbed_count = sum(ndvi == "other_disturbed", na.rm = TRUE), total_count = n(), other_disturbed_prop = other_disturbed_count / total_count, .groups = 'drop')

####Calculate summary statistics for quoll
#Here we’ve created a loop that contains the summary statistics. We subset each individual animal’s trajectory
#and then below each simulated individual’s trajectory and calculate values for all of the summary statistics.
buffer <- 10000
res <- 100
# setup empty objects to store the results
id <- c()
step_length_median <- c()
step_length_mean <- c()
step_length_mean <- c()
step_length_median <- c()
step_length_sd <- c()
rocky_prop<-c()
grassland_prop<-c()
dense_veg_prop<-c()
mine_pit_waste_dump_prop<-c()
other_disturbed_prop<-c()
habitat_dist_mean <- c()
habitat_dist_median <- c()
habitat_dist_sd <- c()
disturb_dist_mean <- c()
disturb_dist_median <- c()
disturb_dist_sd <- c()
gamma_shape <- c()
gamma_scale <- c()
vm_kappa <- c()
straightness <- c()
tot_dist <- c()
intensity_use <- c()
sinuosity <- c()
tac <- c()
#residence_time <- c() # residence time in hours
hr_area_50 <- c()
hr_area_75 <- c()
hr_area_95 <- c()

for(k in 1:length(sim_data_traj_ids)) {
  sim_data_traj <- sim_data_btime %>% filter(id == sim_data_traj_ids[k])
  xmin <- min(sim_data_traj$x2_) - buffer
  xmax <- max(sim_data_traj$x2_) + buffer
  ymin <- min(sim_data_traj$y2_) - buffer
  ymax <- max(sim_data_traj$y2_) + buffer
  template_raster <- rast(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax,
                          res = res, crs = crs("epsg:32751"))
  id[k] <- sim_data_traj$id[1]
  step_length_median[k] <- median(sim_data_traj$sl_)
  step_length_mean[k] <- mean(sim_data_traj$sl_)
  rocky_prop<-rocky_prop_by_id$rocky_prop
  grassland_prop<-grassland_prop_by_id$grassland_prop
  dense_veg_prop<-dense_veg_prop_by_id$dense_veg_prop
  mine_pit_waste_dump_prop<-mine_pit_waste_dump_prop_by_id$mine_pit_waste_dump_prop
  other_disturbed_prop<-other_disturbed_prop_by_id$other_disturbed_prop
  habitat_dist_mean<-mean(sim_data_traj$habitat_distance_end)
  habitat_dist_median<-median(sim_data_traj$habitat_distance_end)
  habitat_dist_sd<-sd(sim_data_traj$habitat_distance_end)
  disturb_dist_mean<-mean(sim_data_traj$disturb_distance_end)
  disturb_dist_median<-median(sim_data_traj$disturb_distance_end)
  disturb_dist_sd<-sd(sim_data_traj$disturb_distance_end)
  gamma_fit <- fit_distr(sim_data_traj$sl_, "gamma")
  gamma_shape[k] <- gamma_fit$params$shape
  gamma_scale[k] <- gamma_fit$params$scale
  vM_fit <- fit_distr(sim_data_traj$sl_, "vonmises")
  vm_kappa[k] <- vM_fit$params$kappa
  straightness[k] <- amt::straightness(sim_data_traj)
  tot_dist[k] <- amt::tot_dist(sim_data_traj)
  intensity_use[k] <- amt::intensity_use(sim_data_traj)
  sinuosity[k] <- amt::sinuosity(sim_data_traj)
  tac[k] <- amt::tac(sim_data_traj)
  #residence_time[k] <- mean(quoll_RT[[k]][,2], na.rm = TRUE)/60/60
  sim_hr_kde <- hr_kde(sim_data_traj, trast = template_raster,
                       levels = c(0.5, 0.75, 0.95))
  sim_hr_kde_area <- hr_area(sim_hr_kde)
  hr_area_50[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.5),]$area/1e6
  hr_area_75[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.75),]$area/1e6
  hr_area_95[k] = sim_hr_kde_area[which(sim_hr_kde_area$level == 0.95),]$area/1e6
}

# create a data frame that has traj, id, and all the summaries
sim_summary_df <- data.frame(
  traj = "obs",
  id = as.numeric(id),
  data = "obs",
  sim = "obs",
  step_length_median = step_length_median,
  step_length_mean = step_length_mean,
  rocky_prop=rocky_prop,
  grassland_prop=grassland_prop,
  dense_veg_prop=dense_veg_prop,
  mine_pit_waste_dump_prop=mine_pit_waste_dump_prop,
  other_disturbed_prop=other_disturbed_prop,
  habitat_dist_mean=habitat_dist_mean,
  habitat_dist_median=habitat_dist_median,
  habitat_dist_sd=habitat_dist_sd,
  disturb_dist_mean=disturb_dist_mean,
  disturb_dist_median=disturb_dist_median,
  disturb_dist_sd = disturb_dist_sd,
  gamma_shape = gamma_shape,
  gamma_scale = gamma_scale,
  vm_kappa = vm_kappa,
  straightness = straightness,
  tot_dist = tot_dist,
  intensity_use = intensity_use,
  sinuosity = sinuosity,
  tac = tac,
  #residence_time = residence_time,
  hr_area_50 = hr_area_50,
  hr_area_75 = hr_area_75,
  hr_area_95 = hr_area_95)


# write_csv(sim_summary_df,
#           paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_3p_memALL_daily_summary_statistics_df_", Sys.Date(), ".csv"))




#####################################################
#### Summarise all statistics and fit PCA ###########
#####################################################

library(tidyverse)
packages <- c("ggh4x", "patchwork", "terra", "ggExtra",
              "cowplot", "ggpubr", "viridis", "scales","readr")
walk(packages, require, character.only = T)
rm(list = ls())


setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/current mining")
ndvi <- rast("Current_habitat_classifications.tif")
habitat_dist <- rast("Current_habitat_distances.tif")
disturb_dist <- rast("Current_disturb_distances.tif")
# Resample habitat_dist and disturb_dist to match ndvi's extent and resolution
habitat_dist_aligned <- terra::resample(habitat_dist, ndvi, method='bilinear') # or 'near' depending on your data
disturb_dist_aligned <- terra::resample(disturb_dist, ndvi, method='bilinear') # or 'near' depending on your data
##Log transform the distance rasters
constant = 0.1
habitat_dist_aligned_log_transformed <- log(habitat_dist_aligned+constant)
plot(habitat_dist_aligned_log_transformed)
constant = 0.1
disturb_dist_aligned_log_transformed <- log(disturb_dist_aligned+constant)
plot(disturb_dist_aligned_log_transformed)

xmin <- ext(ndvi[[1]])[1]
xmax <- ext(ndvi[[1]])[2]
ymin <- ext(ndvi[[1]])[3]
ymax <- ext(ndvi[[1]])[4]

summaries_quoll <-
  read_csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/quoll_summary_statistics_df_2024-03-11.csv")


# read in the simulated data
summaries_0p <-
  read_csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_0p_memALL_daily_summary_statistics_df_2024-03-11.csv")
summaries_0p <- summaries_0p %>% mutate(sim = "0p")
summaries_1p <-
  read_csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_1p_memALL_daily_summary_statistics_df_2024-03-11.csv")
summaries_1p <- summaries_1p %>% mutate(sim = "1p")
summaries_2p <-
  read_csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_2p_memALL_daily_summary_statistics_df_2024-03-11.csv")
summaries_2p <- summaries_2p %>% mutate(sim = "2p")
summaries_3p <-
  read_csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/sim_3p_memALL_daily_summary_statistics_df_2024-03-11.csv")
summaries_3p <- summaries_3p %>% mutate(sim = "3p")

summaries_all <- bind_rows(summaries_quoll,
                                  summaries_0p,
                                  summaries_1p,
                                  summaries_2p,
                                  summaries_3p)
#summaries_all <- summaries_all %>% mutate(residence_time_log = log(residence_time))

summaries_all_long <- summaries_all %>%
  pivot_longer(cols = !c(traj, id, data, sim),
               names_to = "summary", values_to = "value") %>%
  mutate(sim = factor(sim, levels = c("obs", "0p", "1p", "2p", "3p"),
                      labels = c("quoll", "0p", "1p", "2p", "3p")),
         # cateogrise the summary statistics for plotting
         category = case_when(summary == "step_length_median" ~ "step",
                              summary == "step_length_mean" ~ "step",
                              summary == "gamma_shape" ~ "step",
                              summary == "gamma_scale" ~ "step",
                              summary == "vm_kappa" ~ "step",
                              summary == "straightness" ~ "path",
                              summary == "sinuosity" ~ "path",
                              summary == "tac" ~ "life_history",
                              summary == "tot_dist" ~ "life_history",
                              summary == "intensity_use" ~ "life_history",
                              #summary == "residence_time" ~ "life_history",
                              #summary == "residence_time_log" ~ "life_history",
                              summary == "hr_area_50" ~ "life_history",
                              summary == "hr_area_75" ~ "life_history",
                              summary == "hr_area_95" ~ "life_history",
                              summary == "rocky_prop" ~ "habitat",
                              summary == "grassland_prop" ~ "habitat",
                              summary == "dense_veg_prop" ~ "habitat",
                              summary == "mine_pit_waste_dump_prop" ~ "habitat",
                              summary == "other_disturbed_prop" ~ "habitat",
                              summary == "disturb_dist_mean" ~ "habitat",
                              summary == "disturb_dist_median" ~ "habitat",
                              summary == "disturb_dist_sd" ~ "habitat",
                              summary == "habitat_dist_mean" ~ "habitat",
                              summary == "habitat_dist_median" ~ "habitat",
                              summary == "habitat_dist_sd" ~ "habitat"
                              
         ))
head(summaries_all_long)
tail(summaries_all_long)



#####################
######### PCA #######
#####################
# remove some of the columns that should not be included in the PCA, such as id etc,
# and some that are covered by another summary, such as the monthly overlap with
# volume of intersection, as this information is largely contained in the Bhattacharyya's
# affinity overlap
summaries_all_pca <- summaries_all %>%
  dplyr::select(-traj, -id, -data, -sim,-intensity_use,-tac#,-step_length_median,-habitat_dist_median,-disturb_dist_median,-habitat_dist_sd,-disturb_dist_sd,-intensity_use,-tac,-gamma_shape #-residence_time_log, -residence_time
                ) %>%
  na.omit() %>%
  prcomp(center = TRUE, scale. = TRUE)

# summaries_all_pca
#PCA diagnostics
summary(summaries_all_pca)

# biplot(summaries_all_pca)
# screeplot(summaries_all_pca, type = "l")
# proportion of variance explained
prop_var <- summaries_all_pca$sdevˆ2 / sum(summaries_all_pca$sdevˆ2)
# to create of vector of the simulation labels
# this should remove the same columns as above, besides sim
sims <- summaries_all %>%
  dplyr::select(-traj, -id, -data, #-residence_time_log, -residence_time
                ) %>%
  na.omit() %>%
  pull(sim)
# Create a data frame for ggplot
pca_data <- data.frame(sim = sims, summaries_all_pca$x)
# reorder levels for the legend
pca_data$sim <- factor(pca_data$sim, levels = c("obs", "0p", "1p", "2p", "3p"),
                       labels = c("Quoll", "0p", "1p", "2p", "3p"))
# Biplot
ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = pca_data %>% filter(sims != "obs"),
             aes(PC1, PC2, colour = sim),
             size = 1, alpha = 0.5) +
  geom_point(data = pca_data %>% filter(sims == "obs"),
             aes(PC1, PC2),
             colour = "red") +
  scale_colour_viridis_d("Data") +
  ggtitle("PCA Biplot") +
  xlab("Principal Component 1") +
  ylab("Principal Component 2") +
  theme_bw()



# Function to compute convex hull and return a data frame
get_hull_coordinates <- function(data, group) {
  points <- data[data$sim == group,]
  hull <- chull(points$PC1, points$PC2)
  return(data.frame(points[hull, ]))
}

# Compute hulls for each group
hulls <- do.call(rbind, lapply(unique(pca_data$sim), get_hull_coordinates, data = pca_data))
#Add density to the x and y axes.
# Generate Viridis colors for three groups
viridis_colors <- viridis::viridis(4)
# Create a named vector of colors for all groups
colors <- setNames(c('red', viridis_colors), c("Quoll", "0p", "1p", "2p", "3p"))

# Adjusted y-axis limit
pca_plot <- ggplot() +
  geom_point(data = pca_data,
             aes(PC1, PC2, color = sim),
             alpha = 0.5) +
  geom_polygon(data = hulls %>% filter(sim == "Quoll"),
               aes(x = PC1, y = PC2, group = sim),
               alpha = 0.5, fill = "red") +
  theme_bw() +
  ggtitle("PCA of all summary statistics") +
  xlab(paste0("PC1")) +
  ylab(paste0("PC2")) +
  scale_color_manual("Data", values = colors) +
  scale_y_continuous(limits = c(-3, 4)) + # Updated y-axis limit here
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "bottom")

# Marginal density plots remain the same
xdens <-
  axis_canvas(pca_plot, axis = "x") +
  geom_density(data = pca_data, aes(x = PC1, fill = sim), size = 0.25, alpha = 0.25) +
  scale_fill_manual(values = colors)
ydens <-
  axis_canvas(pca_plot, axis = "y", coord_flip = TRUE) +
  geom_density(data = pca_data, aes(x = PC2, fill = sim), size = 0.25, alpha = 0.25) +
  scale_fill_manual(values = colors) +
  coord_flip()

# Combine plots with new y-axis limit
pca_plot %>%
  insert_xaxis_grob(xdens, grid::unit(1.5, "cm"), position = "top") %>%
  insert_yaxis_grob(ydens, grid::unit(1.5, "cm"), position = "right") %>%
  ggdraw()



# #str(hulls)
# 
# # save the plot using ggplot
# ggsave(paste0("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/First simulations for choosing harmonic pairs/Validation summaries/plots/memALL_PCA_all_summaries_",
#               Sys.Date(), ".pdf"), width=150, height=120, units="mm", dpi = 300)
# 
# 


library(ggplot2)
library(grid)

# Extract loadings
loadings <- summaries_all_pca$rotation

# Create a data frame for loadings
loadings_df <- data.frame(Variable = rownames(loadings), PC1 = loadings[, 1], PC2 = loadings[, 2])

# Create the base plot with labels
loadings_plot <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = Variable)) #+
 # geom_text(color = "red", size = 5)  # Adjust the size for better visibility
  
  # Add arrows from the origin to the variable label

loadings_plot <- loadings_plot + geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
                                                arrow = arrow(length = unit(0.2, "inches")),
                                                color = "blue", size = 1)   # Adjust size for better visibility
  

# # Customize the plot
# loadings_plot <- loadings_plot +
#   theme_minimal() +
#   labs(title = "PCA loadings", x = "PC1", y = "PC2") +
#   theme(text = element_text(size = 20),  # Adjust the text size
#         axis.title = element_text(size = 20))  # Adjust the axis title size
# 
# # Print the plot
# print(loadings_plot)
library(ggrepel)
# Assuming loadings_plot is a ggplot object with PCA loadings data created before adding any labels
# Customize the plot with non-overlapping labels
loadings_plot <- loadings_plot +
  theme_minimal() +
  labs(title = "PCA loadings", x = "PC1", y = "PC2") +
  theme(text = element_text(size = 20),  # Adjust the overall text size
        axis.title = element_text(size = 20)) +  # Adjust the axis title size
  geom_text_repel(
    aes(label = Variable),  # Make sure 'Variable' is the name of the column with labels
    size = 5,  # Adjust label text size
    color = "black",  # Ensure labels are black
    max.overlaps = Inf  # Allow for more overlaps if necessary
  )

# Print the plot
print(loadings_plot)

#2
library(ggplot2)
library(ggrepel)

# Create a data frame for loadings
loadings_df <- data.frame(Variable = rownames(loadings), PC1 = loadings[, 1], PC2 = loadings[, 2])

# Calculate the magnitude of the loadings for the color mapping
loadings_df$Magnitude <- sqrt(loadings_df$PC1^2 + loadings_df$PC2^2)

# Create the base plot with labels
loadings_plot <- ggplot(loadings_df, aes(x = PC1, y = PC2, label = Variable))

# Add arrows from the origin to the variable label, with color based on the magnitude
loadings_plot <- loadings_plot + geom_segment(
  aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Magnitude),
  arrow = arrow(length = unit(0.2, "inches")),
  size = 1)  # Adjust size for better visibility

# Map the color to the magma scale
loadings_plot <- loadings_plot + scale_color_gradientn(colors = viridis::magma(10))

# Customize the plot with non-overlapping labels
loadings_plot <- loadings_plot +
  theme_minimal() +
  labs(title = "", x = "PC1", y = "PC2") +
  theme(
    text = element_text(size = 20),  # Adjust the overall text size
    axis.title = element_text(size = 20),  # Adjust the axis title size
    legend.title = element_blank()  # Remove the legend title if desired
  ) +
  geom_text_repel(
    aes(label = Variable),  # Make sure 'Variable' is the name of the column with labels
    size = 5,  # Adjust label text size
    color = "black",  # Ensure labels are black
    max.overlaps = Inf  # Allow for more overlaps if necessary
  )

# Print the plot
print(loadings_plot)