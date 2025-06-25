##################
######Distance raster for temporal effects
##################
###Modelling for SSF
library(amt)
library(MuMIn)
library(gdata)
library(survival)
library(ggplot2)
library(tidyverse)
library(AICcmodavg)
library(lubridate)

rm(list = ls())
#dev.off()

setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")

all_ssf1 <-read.csv("Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")
names(all_ssf1)


########################
#One pair of harmonics
########################
library(amt)
library(MuMIn)
library(gdata)
library(survival)
library(ggplot2)
library(tidyverse)
library(AICcmodavg)
library(lubridate)

#dev.off()

setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")
#all_ssf1 <-read.csv("Woodie_SSF_w_perc.csv")
all_ssf1 <-read.csv("Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")
names(all_ssf1)

# Assuming all_ssf1 is your dataframe and t1_ is the column with timestamps
all_ssf1$t1_ <- as.POSIXct(all_ssf1$t1_, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")
all_ssf1$t2_ <- as.POSIXct(all_ssf1$t2_, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")

# Load the lubridate package
library(lubridate)
all_ssf1$realtime <- all_ssf1$t2_
# Extract the hour
hour_component <- hour(all_ssf1$t2_)
# Extract the minute and determine if .5 should be added
minute_component <- ifelse(minute(all_ssf1$t2_) == 30, 0.5, 0)
# Combine hour and minute components
all_ssf1$hour_from_t2_ <- hour_component + minute_component
# View the modified dataframe
head(all_ssf1)
# 
# # Define the original hours in order
# original_hours <- c(18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
# # Create a sequence from 0 to 11.5 in steps of 0.5
# new_sequence <- seq(1, 24, by=1)
# # Create a mapping from original hours to new sequence
# hour_mapping <- setNames(new_sequence, original_hours)
# # Assuming hourly_coefs is your dataframe
# # Map the original hour values to the new sequence
# all_ssf1$newhour <- hour_mapping[as.character(all_ssf1$hour_from_t2_)]
constant=0.1
#all_ssf1$t1_<-dmy_hm(all_ssf1$t1_)
#all_ssf1$t2_<-dmy_hm(all_ssf1$t2_)
all_ssf2 <- all_ssf1 %>%
  mutate(id_num = as.numeric(factor(id)),
         step_id = step_id_,
         burst_ = burst_,
         x1 = x1_, x2 = x2_,
         y1 = y1_, y2 = y2_,
         t1 = t1_,
         t1_rounded = round_date(all_ssf1$t1_, "hour"),
         hour_t1 = hour(t1_rounded),
         t2 = t2_,
         t2_rounded = round_date(all_ssf1$t2_, "hour"),
         hour_t2 = hour(t2_rounded),
         hour = hour_t2,
         yday = yday(t1_),
         year = year(t1_),
         month = month(t1_),
         sl = sl_,
         log_sl = log(sl_),
         ta = ta_,
         cos_ta = cos(ta_),
         NDVI_class_end = NDVI_class_end,
         disturb_distance_end = disturb_distance_end,
         log_disturb = log(disturb_distance_end+constant),
         log_habitat = log(habitat_distance_end+constant),
         habitat_distance_end = habitat_distance_end,
         name=name,orig_name=orig_name,season=season,age=age,
         spatial_memory_ref = kde_ref_spatial_memory_density_log,
         spatial_memory_pi = kde_pi_spatial_memory_density_log,
         memory_log=kde_ref_spatial_memory_density_log,
         memory=exp(kde_ref_spatial_memory_density_log),
         hour_s1 = sin(2*pi*hour/24),
         hour_s2 = sin(4*pi*hour/24),
         hour_s3 = sin(6*pi*hour/24),
         hour_s4 = sin(8*pi*hour/24),
         hour_c1 = cos(2*pi*hour/24),
         hour_c2 = cos(4*pi*hour/24),
         hour_c3 = cos(6*pi*hour/24),
         hour_c4 = cos(8*pi*hour/24))


##Look at habitat distance through time (days)
ggplot(data = all_ssf2 %>% filter(case_ == 1,season=="breeding"),
       aes(x = t1_, y = habitat_distance_end, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic()

##Look at habitat distance through time (days)
ggplot(data = all_ssf2 %>% filter(case_ == 1,season=="breeding"),
       aes(x = t1_, y = memory, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic()

##Look at habitat distance through time (hours)
all_ssf2 <- all_ssf2 %>%
  mutate(hour_transformed = ifelse(hour < 6, hour + 24, hour))

ggplot(data = all_ssf2 %>% filter(case_ == 1, season == "breeding"),
       aes(x = hour_transformed, y = habitat_distance_end, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic() +
  coord_cartesian(xlim = c(18, 30)) # Adjusted based on the transformed hours


ggplot(data = all_ssf2 %>% filter(case_ == 1, season == "breeding"),
       aes(x = hour_transformed, y = memory, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic() +
  coord_cartesian(xlim = c(18, 30)) # Adjusted based on the transformed hours


##Separate seasons for separate modelling
################For breeding season
all_ssf_breeding <- all_ssf2 %>% 
  filter(season == "breeding") %>% 
  arrange(t1_)

all_ssf_breeding$NDVI_class_end <-factor(all_ssf_breeding$NDVI_class_end, levels = c("rocky", "grassland", "dense_veg","mine_pit_waste_dump","other_disturbed"))

unique(all_ssf_breeding$NDVI_class_end)


quoll_ids <- unique(all_ssf_breeding$id)


quoll_data_matrix_unscaled <- all_ssf_breeding %>% transmute(
  
  ndvi = NDVI_class_end,
  
  habitat_distance = log_habitat,
  habitat_distance_s1 = log_habitat * hour_s1,
  habitat_distance_c1 = log_habitat * hour_c1,
  # habitat_distance_s2 = log_habitat * hour_s2, # New
  # habitat_distance_c2 = log_habitat * hour_c2, # New
  
  disturb_distance = log_disturb,
  disturb_distance_s1 = log_disturb * hour_s1,
  disturb_distance_c1 = log_disturb * hour_c1,
  # disturb_distance_s2 = log_disturb * hour_s2, # New
  # disturb_distance_c2 = log_disturb * hour_c2, # New
  
  memory = memory_log,
  memory_s1 = memory_log * hour_s1,
  memory_c1 = memory_log * hour_c1,
  # memory_s2 = memory * hour_s2, # New
  # memory_c2 = memory * hour_c2, # New
  
  #spatial_memory_pi = spatial_memory_pi ^ 2,
  #spatial_memory_pi_s1 = (spatial_memory_pi ^ 2) * hour_s1,
  #spatial_memory_pi_c1 = (spatial_memory_pi ^ 2) * hour_c1,
  #spatial_memory_pi_s2 = (spatial_memory_pi ^ 2) * hour_s2, # New
  #spatial_memory_pi_c2 = (spatial_memory_pi ^ 2) * hour_c2, # New
  
  step_l = sl,
  # step_l_s1 = sl * hour_s1,
  # step_l_c1 = sl * hour_c1,
  # step_l_s2 = sl * hour_s2, # New
  # step_l_c2 = sl * hour_c2, # New
  
  log_step_l = log_sl,
  # log_step_l_s1 = log_sl * hour_s1,
  # log_step_l_c1 = log_sl * hour_c1,
  # log_step_l_s2 = log_sl * hour_s2, # New
  # log_step_l_c2 = log_sl * hour_c2, # New
  # 
  cos_turn_a = cos_ta,
  # cos_turn_a_s1 = cos_ta * hour_s1,
  # cos_turn_a_c1 = cos_ta * hour_c1,
  # cos_turn_a_s2 = cos_ta * hour_s2, # New
  # cos_turn_a_c2 = cos_ta * hour_c2 # New
)


head(quoll_data_matrix_unscaled)
quoll_data_matrix_unscaled<-as.data.frame(quoll_data_matrix_unscaled)
## Identify numeric columns
#numeric_columns <- sapply(quoll_data_matrix_unscaled, is.numeric)
## Scale only numeric columns
#scaled_data <- scale(quoll_data_matrix_unscaled[, numeric_columns])
## Combine scaled numeric data with non-numeric columns
#quoll_data_matrix_scaled <- cbind(quoll_data_matrix_unscaled[!numeric_columns], scaled_data)
## Convert row names to the first column if necessary
#quoll_data_matrix_scaled <- data.frame(quoll_data_matrix_scaled, row.names = NULL)
## View the head of the scaled dataset
#head(quoll_data_matrix_scaled)

# Using base R
quoll_data_matrix_unscaled_no_ndvi <- quoll_data_matrix_unscaled[ , !(names(quoll_data_matrix_unscaled) %in% c("ndvi","step_l","log_step_l","cos_turn_a"))]
quoll_data_matrix_scaled_no_ndvi <- scale(quoll_data_matrix_unscaled_no_ndvi)

quoll_data_matrix_scaled <- cbind(ndvi = quoll_data_matrix_unscaled$ndvi, quoll_data_matrix_scaled_no_ndvi)


mean_vals <- attr(quoll_data_matrix_scaled_no_ndvi, "scaled:center")
sd_vals <- attr(quoll_data_matrix_scaled_no_ndvi, "scaled:scale")
scaling_attributes <- data.frame(variable = names(quoll_data_matrix_unscaled_no_ndvi), mean = mean_vals, sd = sd_vals)

quoll_data_scaled_1p <- data.frame(id = all_ssf_breeding$id,  step_id = all_ssf_breeding$step_id, case_ = all_ssf_breeding$case_, quoll_data_matrix_scaled)
head(quoll_data_scaled_1p)

# Ensure the row order and number of rows in both datasets match
# Replace the ndvi column in quoll_data_scaled_1p with that from quoll_data_matrix_unscaled
quoll_data_scaled_1p$ndvi <- quoll_data_matrix_unscaled$ndvi
quoll_data_scaled_1p$step_l <- quoll_data_matrix_unscaled$step_l
quoll_data_scaled_1p$log_step_l <- quoll_data_matrix_unscaled$log_step_l
quoll_data_scaled_1p$cos_turn_a <- quoll_data_matrix_unscaled$cos_turn_a


# View the first few rows to confirm the change
head(quoll_data_scaled_1p)



formula_twostep <- case_ ~ 
  
  ndvi +
  
  habitat_distance +
  habitat_distance_s1 +
  habitat_distance_c1 +
  # habitat_distance_s2 + # Added
  # habitat_distance_c2 + # Added
  
  disturb_distance +
  disturb_distance_s1 +
  disturb_distance_c1 +
  # disturb_distance_s2 + # Added
  # disturb_distance_c2 + # Added
  
  memory +
  memory_s1 +
  memory_c1 +
  # memory_s2 + # Added
  # memory_c2 + # Added
  
  #spatial_memory_pi +
  #spatial_memory_pi_s1 +
  #spatial_memory_pi_c1 +
  #spatial_memory_pi_s2 + # Added
  #spatial_memory_pi_c2 + # Added
  
  step_l +
  # step_l_s1 +
  # step_l_c1 +
  # step_l_s2 + # Added
  # step_l_c2 + # Added
  
  log_step_l +
  # log_step_l_s1 +
  # log_step_l_c1 +
  # log_step_l_s2 + # Added
  # log_step_l_c2 + # Added
  
  cos_turn_a +
  # cos_turn_a_s1 +
  # cos_turn_a_c1 +
  # cos_turn_a_s2 + # Added
  # cos_turn_a_c2 + # Added
  
  strata(step_id) +
  cluster(id)


library(tictoc)
library(TwoStepCLogit)

  tic()
  model_twostep_1p_harms <- Ts.estim(formula = formula_twostep,
                                     data = quoll_data_scaled_1p,
                                     all.m.1 = TRUE,
                                     D = "UN(1)",
                                     itermax = 10000)
  toc()
# 

#summary(model_twostep_1p_harms)
model_twostep_1p_harms
model_twostep_1p_harms$beta
model_twostep_1p_harms$se
model_twostep_1p_harms$vcov
diag(model_twostep_1p_harms$D) # between cluster variance
model_twostep_1p_harms$r.effect # individual estimates

hist(model_twostep_1p_harms$r.effect[,6])


coefs_clr <- data.frame(coefs = names(model_twostep_1p_harms$beta), value = model_twostep_1p_harms$beta)
coefs_clr_no_ndvi <- coefs_clr[!coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), ]

coefs_clr_no_ndvi$scale_sd <- scaling_attributes$sd
coefs_clr_no_ndvi <- coefs_clr_no_ndvi %>% mutate(value_nat = value / scale_sd)

head(coefs_clr)
head(coefs_clr_no_ndvi)

# Rows to add back
rows_to_add <- coefs_clr[coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), ]
rows_to_add$scale_sd <- NA
rows_to_add$value_nat <- coefs_clr[coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), "value"]
coefs_clr_upd <- rbind(coefs_clr_no_ndvi, rows_to_add)
print(coefs_clr_upd)

####Reconstructing coefficients with two pairs of harmonics, with quadratic terms

# we want to reconstruct the temporally dynamic coefficients for the full cycle, and then we can subset to the active hours after
hour <- seq(0.0,23.5,0.5)

all_ssf_breeding$hour30 <- format(as.POSIXct(all_ssf_breeding$t2), format = "%H:%M")
all_ssf_breeding$hour30 <- gsub(":", ".", all_ssf_breeding$hour30)
all_ssf_breeding$hour30 <- sub("30", "50", all_ssf_breeding$hour30)
all_ssf_breeding$hour30 <- as.numeric(all_ssf_breeding$hour30)

unique_hours<-unique(all_ssf_breeding$hour30)


# this creates a matrix of the hour harmonics, which is 24 x 3 (hours x harmonics, including the linear term)
hour_harmonics_df <- data.frame("linear_term" = rep(1, length(hour)),
                                "hour_s1" = sin(2*pi*hour/24),
                                "hour_s2" = sin(4*pi*hour/24))

hour_harmonics_linear_df <- data.frame("linear_term" = rep(1, length(hour)))


# we'll do the filtering after reconstructing the harmonics - i think it's best to do the whole cycle initially
# hour_harmonics_df <- hour_harmonics_df[hour_harmonics_df$hour %in% unique_hours, ]

hour_transformed = ifelse(hour < 6, hour + 24, hour)

# i spaced the coefficients out so i could see them a bit easier
harmonics_scaled_df_1p <- data.frame(
  "hour" = hour,
  # can add this in here now
  "hour_transformed" = ifelse(hour < 6, hour + 24, hour),
  
  # note that for the linear only terms i've added in the LINEAR ONLY data frame here, as there is only the linear term, and the matrix multiplication will work (now it's just an 1 x 1 %*% 1 x n, which will result in a 1 x n matrix (i.e. a vector))
  "ndvigrassland" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvigrassland", coefs)) %>%
                                 ###### under here ######
                               pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvidense_veg" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvidense_veg", coefs)) %>%
                                 pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvimine_pit_waste_dump" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvimine_pit_waste_dump", coefs)) %>%
                                           pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndviother_disturbed" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndviother_disturbed", coefs)) %>%
                                       pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  # For these ones, we will pull out the three coefficients for each variable (linear, sin1,cos1), which R considers to be COLUMN, and therefore a 3 x 1 matrix (or more simply a (column) vector). We transpose (the t() function) it to make it a 1 x 3. We then have the 24 x 3 matrix of the harmonics, which we TRANSPOSE to turn it into a 3 x 24 matrix. Then when we multiply them together, we have 1 x 3 %*% 3 x 24, which results in a 1 x 24 matrix (i.e. a VECTOR of the temporally dynamic coefficient throughout the day), which we add into the dataframe
  
  "step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("step_l", coefs) & !grepl("log", coefs)) %>% 
                          pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  # add in the log_sl term
  "log_step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("log_step_l", coefs)) %>% 
                              pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "cos_turn_a" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("cos_turn_a", coefs)) %>% 
                              pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "disturb_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("disturb_distance", coefs) & !grepl("sq", coefs)) %>%
                                    pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "habitat_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("habitat_distance", coefs)) %>% 
                                    pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "memory" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("memory", coefs)) %>% 
                          pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
  #,
  
  #"spatial_memory_pi" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("spatial_memory_pi", coefs)) %>% 
  #                                   pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
)



harmonics_nat_df_1p <- data.frame(
  "hour" = hour,
  # can add this in here now
  "hour_transformed" = ifelse(hour < 6, hour + 24, hour),
  
  # note that for the linear only terms i've added in the LINEAR ONLY data frame here, as there is only the linear term, and the matrix multiplication will work (now it's just an 1 x 1 %*% 1 x n, which will result in a 1 x n matrix (i.e. a vector))
  "ndvigrassland" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvigrassland", coefs)) %>%
                                 ###### under here ######
                               pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvidense_veg" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvidense_veg", coefs)) %>%
                                 pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvimine_pit_waste_dump" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvimine_pit_waste_dump", coefs)) %>%
                                           pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndviother_disturbed" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndviother_disturbed", coefs)) %>%
                                       pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  # For these ones, we will pull out the three coefficients for each variable (linear, sin1,cos1), which R considers to be COLUMN, and therefore a 3 x 1 matrix (or more simply a (column) vector). We transpose (the t() function) it to make it a 1 x 3. We then have the 24 x 3 matrix of the harmonics, which we TRANSPOSE to turn it into a 3 x 24 matrix. Then when we multiply them together, we have 1 x 3 %*% 3 x 24, which results in a 1 x 24 matrix (i.e. a VECTOR of the temporally dynamic coefficient throughout the day), which we add into the dataframe
  
  "step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("step_l", coefs) & !grepl("log", coefs)) %>% 
                          pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  # add in the log_sl term
  "log_step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("log_step_l", coefs)) %>%
                              pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "cos_turn_a" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("cos_turn_a", coefs)) %>%
                              pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "disturb_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("disturb_distance", coefs) & !grepl("sq", coefs)) %>%
                                    pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "habitat_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("habitat_distance", coefs)) %>% 
                                    pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "memory" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("memory", coefs)) %>% 
                          pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df)))
  #,
  
  #"spatial_memory_pi" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("spatial_memory_pi", coefs)) %>% 
  #                                   pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
)




##Get the tentative distributions
all_ssf_breed_true <- all_ssf_breeding %>% 
  filter(case_ == "1") %>% 
  arrange(t1_)
# Fit a gamma distribution to the step length of true steps
fit_sl <- fit_distr(all_ssf_breed_true$sl_, "gamma")
# Extract the shape and rate parameters of the gamma distribution
shape_sl <- fit_sl$params[1]
scale_sl <- fit_sl$params[2]
# Load the circular package
library(circular)
# Fit a von Mises distribution to the ta_ values
fit_ta <- fit_distr(all_ssf_breed_true$ta_, "vonmises")
# Extract the kappa parameter of the von Mises distribution
kappa_ta <- fit_ta$params["kappa"]


tentative_shape <- shape_sl$shape
tentative_scale <- scale_sl$scale
tentative_kappa <- kappa_ta$kappa

hist(rgamma(1e4,
            shape = tentative_shape,
            scale = tentative_scale),
     breaks = 100)

harmonics_nat_scaled_df_1p <- harmonics_nat_df_1p %>% mutate(shape = tentative_shape + log_step_l,
                                                             scale = 1/((1/tentative_scale) - step_l),
                                                             kappa = tentative_kappa + cos_turn_a)

hist(rgamma(1e4,
            shape = harmonics_nat_scaled_df_1p$shape[1],
            scale = harmonics_nat_scaled_df_1p$scale[1]),
     breaks = 100)




harmonics_nat_scaled_df_1p
# Add a new column 'rocky' with 0s to the dataframe 'harmonics_nat_scaled_df_1p'
harmonics_nat_scaled_df_1p$ndvirocky <- 0

harmonics_scaled_long_1p <- pivot_longer(harmonics_nat_scaled_df_1p, cols = !1:2, names_to = "coef")

ggplot() +
  geom_path(data = harmonics_scaled_long_1p,
            aes(x = hour, y = value, colour = coef)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour") +
  scale_color_discrete("Estimate") +
  theme_classic() +
  theme(legend.position = "bottom")


# when removing the daytime hours
harmonics_scaled_df_1p_active <- harmonics_nat_scaled_df_1p %>% dplyr::filter(hour %in% unique_hours)
harmonics_scaled_long_1p_active <- pivot_longer(harmonics_scaled_df_1p_active, cols = !1:2, names_to = "coef")
##write.csv(harmonics_scaled_df_1p_active,"C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_3hrm.csv")
movement_summary <- all_ssf_breeding %>% filter(case_ == 1) %>%  group_by(hour30) %>% summarise(mean_sl = mean(sl), median_sl = median(sl))
movement_summary$hour <- movement_summary$hour30

final1<-cbind(harmonics_scaled_df_1p_active,movement_summary)
final1<-final1[1:20]
final1<-final1 %>% dplyr::select(-hour30)



#write.csv(final1,"C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_1hrm.csv",row.names = FALSE)

harmonics_scaled_df_1p_active<-read.csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_1hrm.csv")

str(harmonics_scaled_long_1p_active)
head(harmonics_scaled_long_1p_active)


library(ggplot2)
# Generate a set of random colors
# First, find the number of unique levels in the 'coef' column
num_levels <- length(unique(harmonics_scaled_long_1p_active$coef))

# Then, generate random colors for each level
# You can adjust the seed for reproducibility
set.seed(123) # Remove or change the seed value for different colors
random_colors <- grDevices::rainbow(num_levels)

# Now, create the plot
ggplot(harmonics_scaled_long_1p_active, aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(size = 2.5) +  # Increased line thickness
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5","6"),
                     limits = c(19, 30)) +
  scale_color_manual("Estimate", values = setNames(random_colors, unique(harmonics_scaled_long_1p_active$coef))) +
  theme_classic() +
  theme(legend.position = "bottom")

# Assuming 'filtered_data' is your dataset and it's already loaded into your R session
# Now, 'hour_transformed' will have 30 instead of 6 wherever applicable


# Load necessary libraries
library(dplyr)
library(ggplot2)

# Step 1: Filter out specific levels from the dataset
filtered_data <- harmonics_scaled_long_1p_active %>%
  filter(!coef %in% c("shape", "scale", "kappa", "memory", "step_l", "cos_turn_a", "log_step_l","disturb_distance","habitat_distance"))

# Assuming 'filtered_data' is your pre-processed dataset

# Load necessary libraries
library(ggplot2)
library(viridis) # for magma palette

# Modify 'hour_transformed' so that every instance of 6 becomes 30
filtered_data <- filtered_data %>%
  mutate(hour_transformed = ifelse(hour_transformed == 6, 30, hour_transformed))

# Assuming 'filtered_data' and 'magma_colors' are already defined and you've loaded the dplyr library

# Update 'coef' values with new category names
filtered_data <- filtered_data %>%
  mutate(coef = case_when(
    coef == "ndvidense_veg"           ~ "Riparian habitat",
    coef == "ndvigrassland"           ~ "Spinifex grassland",
    coef == "ndvirocky"               ~ "Rocky habitat",
    coef == "ndvimine_pit_waste_dump" ~ "Mine pits and waste dumps",
    coef == "ndviother_disturbed"     ~ "Other disturbed land",
    TRUE                              ~ coef
  )) %>%
  # Convert 'coef' to a factor with levels in the specified order
  mutate(coef = factor(coef, levels = c(
    "Riparian habitat", 
    "Spinifex grassland", 
    "Rocky habitat", 
    "Mine pits and waste dumps", 
    "Other disturbed land"
  )))

# Now that 'coef' values are updated, regenerate the magma_colors mapping for the updated levels
num_levels_updated <- length(unique(filtered_data$coef))
magma_colors_updated <- setNames(viridis(num_levels_updated), unique(filtered_data$coef))

#Create the plot with updated 'coef' values and corresponding colors, doubling text and line sizes
ggplot(filtered_data, aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(linewidth = 5) +  # Double the line thickness
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 2) +  # Double the line thickness for hline
  scale_y_continuous(expression(Coefficients~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("18", "19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5", "6"),
                     limits = c(18, 30)) +
  scale_color_manual(name = "Covariate", values = magma_colors_updated) +
  theme_classic(base_size = 44) +  # Adjusted the base text size
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 88),  # Double the plot title size based on the new base size
    axis.title = element_text(size = 44),  # Adjusted axis title size
    axis.text = element_text(size = 44),  # Adjusted axis text size
    legend.title = element_text(size = 40),  # Adjusted legend title size
    legend.text = element_text(size = 36),  # Adjusted legend text size
    legend.key.size = unit(2, "lines")  # Adjust legend key size if necessary
  ) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE)) # Adjust the legend to spread over two rows





# Step 1: Filter out specific levels from the dataset
filtered_data <- harmonics_scaled_long_1p_active %>%
  filter(!coef %in% c("ndvigrassland", "ndvirocky", "ndvidense_veg", "memory", "ndvimine_pit_waste_dump", "ndviother_disturbed", "log_step_l","cos_turn_a","step_l","kappa","scale","shape"))

# Assuming 'filtered_data' is your pre-processed dataset

# Load necessary libraries
library(ggplot2)
library(viridis) # for magma palette

# Modify 'hour_transformed' so that every instance of 6 becomes 30
filtered_data <- filtered_data %>%
  mutate(hour_transformed = ifelse(hour_transformed == 6, 30, hour_transformed))

# Assuming 'filtered_data' and 'magma_colors' are already defined and you've loaded the dplyr library

#Update 'coef' values with new category names
filtered_data <- filtered_data %>%
  mutate(coef = case_when(
    coef == "habitat_distance"           ~ "Distance from rocky habitat",
    coef == "disturb_distance"           ~ "Distance from mining disturbance",
    TRUE                              ~ coef
  )) %>%
  # Convert 'coef' to a factor with levels in the specified order
  mutate(coef = factor(coef, levels = c(
    "Distance from rocky habitat",
    "Distance from mining disturbance"
  )))

# Now that 'coef' values are updated, regenerate the magma_colors mapping for the updated levels
num_levels_updated <- length(unique(filtered_data$coef))
magma_colors_updated <- setNames(viridis(num_levels_updated), unique(filtered_data$coef))

# Create the plot with updated 'coef' values and corresponding colors
ggplot(filtered_data, aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(linewidth = 2.5) +  # Adjusted line thickness
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  scale_y_continuous(expression(Coefficients~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("18", "19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5", "6"),
                     limits = c(18, 30)) +
  scale_color_manual(name = "Estimate", values = magma_colors_updated) +
  theme_classic(base_size = 22) +  # Adjusted the base text size
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 44),  # Adjusted plot title size
    axis.title = element_text(size = 22),  # Adjusted axis title size
    axis.text = element_text(size = 22),  # Adjusted axis text size
    legend.title = element_text(size = 20),  # Adjusted legend title size
    legend.text = element_text(size = 18)  # Adjusted legend text size
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) # Adjust the legend to spread over two rows









# Step 1: Filter out specific levels from the dataset
filtered_data <- harmonics_scaled_long_1p_active %>%
  filter(!coef %in% c("ndvigrassland", "ndvirocky", "ndvidense_veg", "habitat_distance", "disturb_distance", "ndvimine_pit_waste_dump", "ndviother_disturbed", "log_step_l","cos_turn_a","step_l","kappa","scale","shape"))

# Assuming 'filtered_data' is your pre-processed dataset

# Load necessary libraries
library(ggplot2)
library(viridis) # for magma palette

# Modify 'hour_transformed' so that every instance of 6 becomes 30
filtered_data <- filtered_data %>%
  mutate(hour_transformed = ifelse(hour_transformed == 6, 30, hour_transformed))

# Assuming 'filtered_data' and 'magma_colors' are already defined and you've loaded the dplyr library

#Update 'coef' values with new category names
filtered_data <- filtered_data %>%
  mutate(coef = case_when(
    coef == "memory"           ~ "Memory process",
    TRUE                              ~ coef
  )) %>%
  # Convert 'coef' to a factor with levels in the specified order
  mutate(coef = factor(coef, levels = c(
    "Memory process"
      )))

# Now that 'coef' values are updated, regenerate the magma_colors mapping for the updated levels
num_levels_updated <- length(unique(filtered_data$coef))
magma_colors_updated <- setNames(viridis(num_levels_updated), unique(filtered_data$coef))

# Create the plot with updated 'coef' values and corresponding colors
ggplot(filtered_data, aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(linewidth = 2.5) +  # Adjusted line thickness
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 1) +
  scale_y_continuous(expression(Coefficients~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("18", "19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5", "6"),
                     limits = c(18, 30)) +
  scale_color_manual(name = "Estimate", values = magma_colors_updated) +
  theme_classic(base_size = 22) +  # Adjusted the base text size
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 44),  # Adjusted plot title size
    axis.title = element_text(size = 22),  # Adjusted axis title size
    axis.text = element_text(size = 22),  # Adjusted axis text size
    legend.title = element_text(size = 20),  # Adjusted legend title size
    legend.text = element_text(size = 18)  # Adjusted legend text size
  ) +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) # Adjust the legend to spread over two rows







# Step 1: Filter out specific levels from the dataset
filtered_data <- harmonics_scaled_long_1p_active %>%
  filter(!coef %in% c("ndvigrassland", "ndvirocky", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","scale","shape","kappa","cos_turn_a","step_l","log_step_l"))

# Assuming 'filtered_data' is your pre-processed dataset

# Load necessary libraries
library(ggplot2)
library(viridis) # for magma palette

# Modify 'hour_transformed' so that every instance of 6 becomes 30
filtered_data <- filtered_data %>%
  mutate(hour_transformed = ifelse(hour_transformed == 6, 30, hour_transformed))

# Assuming 'filtered_data' and 'magma_colors' are already defined and you've loaded the dplyr library

#Update 'coef' values with new category names
filtered_data <- filtered_data %>%
  mutate(coef = case_when(
    coef == "memory"   ~ "Memory process",
    coef == "habitat_distance"           ~ "Distance from rocky habitat",
    coef == "disturb_distance"           ~ "Distance from mining",
    TRUE                              ~ coef
  )) %>%
  # Convert 'coef' to a factor with levels in the specified order
  mutate(coef = factor(coef, levels = c(
    "Distance from rocky habitat",
    "Distance from mining",
    "Memory process"
  )))

filtered_data$scaled <- scale(filtered_data$value)

# Now that 'coef' values are updated, regenerate the magma_colors mapping for the updated levels
num_levels_updated <- length(unique(filtered_data$coef))
magma_colors_updated <- setNames(viridis(num_levels_updated), unique(filtered_data$coef))

#Create the plot with updated 'coef' values and corresponding colors, doubling text and line sizes
ggplot(filtered_data, aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(linewidth = 5) +  # Double the line thickness
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 2) +  # Double the line thickness for hline
  scale_y_continuous(expression(Coefficients~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("18", "19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5", "6"),
                     limits = c(18, 30)) +
  scale_color_manual(name = "Covariate", values = magma_colors_updated) +
  theme_classic(base_size = 44) +  # Adjusted the base text size
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 88),  # Double the plot title size based on the new base size
    axis.title = element_text(size = 44),  # Adjusted axis title size
    axis.text = element_text(size = 44),  # Adjusted axis text size
    legend.title = element_text(size = 40),  # Adjusted legend title size
    legend.text = element_text(size = 36),  # Adjusted legend text size
    legend.key.size = unit(2, "lines")  # Adjust legend key size if necessary
  ) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE)) # Adjust the legend to spread over two rows

#



















########################
#Two pairs of harmonics
########################
library(amt)
library(MuMIn)
library(gdata)
library(survival)
library(ggplot2)
library(tidyverse)
library(AICcmodavg)
library(lubridate)

#dev.off()

setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")
#all_ssf1 <-read.csv("Woodie_SSF_w_perc.csv")
all_ssf1 <-read.csv("Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")
names(all_ssf1)

# Assuming all_ssf1 is your dataframe and t1_ is the column with timestamps
all_ssf1$t1_ <- as.POSIXct(all_ssf1$t1_, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")
all_ssf1$t2_ <- as.POSIXct(all_ssf1$t2_, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")

# Load the lubridate package
library(lubridate)
all_ssf1$realtime <- all_ssf1$t2_
# Extract the hour
hour_component <- hour(all_ssf1$t2_)
# Extract the minute and determine if .5 should be added
minute_component <- ifelse(minute(all_ssf1$t2_) == 30, 0.5, 0)
# Combine hour and minute components
all_ssf1$hour_from_t2_ <- hour_component + minute_component
# View the modified dataframe
head(all_ssf1)
# 
# # Define the original hours in order
# original_hours <- c(18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
# # Create a sequence from 0 to 11.5 in steps of 0.5
# new_sequence <- seq(1, 24, by=1)
# # Create a mapping from original hours to new sequence
# hour_mapping <- setNames(new_sequence, original_hours)
# # Assuming hourly_coefs is your dataframe
# # Map the original hour values to the new sequence
# all_ssf1$newhour <- hour_mapping[as.character(all_ssf1$hour_from_t2_)]
constant=0.1
#all_ssf1$t1_<-dmy_hm(all_ssf1$t1_)
#all_ssf1$t2_<-dmy_hm(all_ssf1$t2_)
all_ssf2 <- all_ssf1 %>%
  mutate(id_num = as.numeric(factor(id)),
         step_id = step_id_,
         burst_ = burst_,
         x1 = x1_, x2 = x2_,
         y1 = y1_, y2 = y2_,
         t1 = t1_,
         t1_rounded = round_date(all_ssf1$t1_, "hour"),
         hour_t1 = hour(t1_rounded),
         t2 = t2_,
         t2_rounded = round_date(all_ssf1$t2_, "hour"),
         hour_t2 = hour(t2_rounded),
         hour = hour_t2,
         yday = yday(t1_),
         year = year(t1_),
         month = month(t1_),
         sl = sl_,
         log_sl = log(sl_),
         ta = ta_,
         cos_ta = cos(ta_),
         NDVI_class_end = NDVI_class_end,
         disturb_distance_end = disturb_distance_end,
         log_disturb = log(disturb_distance_end+constant),
         log_habitat = log(habitat_distance_end+constant),
         habitat_distance_end = habitat_distance_end,
         name=name,orig_name=orig_name,season=season,age=age,
         spatial_memory_ref = kde_ref_spatial_memory_density_log,
         spatial_memory_pi = kde_pi_spatial_memory_density_log,
         memory_log=kde_ref_spatial_memory_density_log,
         memory=exp(kde_ref_spatial_memory_density_log),
         hour_s1 = sin(2*pi*hour/24),
         hour_s2 = sin(4*pi*hour/24),
         hour_s3 = sin(6*pi*hour/24),
         hour_s4 = sin(8*pi*hour/24),
         hour_c1 = cos(2*pi*hour/24),
         hour_c2 = cos(4*pi*hour/24),
         hour_c3 = cos(6*pi*hour/24),
         hour_c4 = cos(8*pi*hour/24))


##Look at habitat distance through time (days)
ggplot(data = all_ssf2 %>% filter(case_ == 1,season=="breeding"),
       aes(x = t1_, y = habitat_distance_end, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic()

##Look at habitat distance through time (days)
ggplot(data = all_ssf2 %>% filter(case_ == 1,season=="breeding"),
       aes(x = t1_, y = memory, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic()

##Look at habitat distance through time (hours)
all_ssf2 <- all_ssf2 %>%
  mutate(hour_transformed = ifelse(hour < 6, hour + 24, hour))

ggplot(data = all_ssf2 %>% filter(case_ == 1, season == "breeding"),
       aes(x = hour_transformed, y = habitat_distance_end, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic() +
  coord_cartesian(xlim = c(18, 30)) # Adjusted based on the transformed hours


ggplot(data = all_ssf2 %>% filter(case_ == 1, season == "breeding"),
       aes(x = hour_transformed, y = memory, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic() +
  coord_cartesian(xlim = c(18, 30)) # Adjusted based on the transformed hours


##Separate seasons for separate modelling
################For breeding season
all_ssf_breeding <- all_ssf2 %>% 
  filter(season == "breeding") %>% 
  arrange(t1_)

all_ssf_breeding$NDVI_class_end <-factor(all_ssf_breeding$NDVI_class_end, levels = c("rocky", "grassland", "dense_veg","mine_pit_waste_dump","other_disturbed"))

unique(all_ssf_breeding$NDVI_class_end)


quoll_ids <- unique(all_ssf_breeding$id)


quoll_data_matrix_unscaled <- all_ssf_breeding %>% transmute(
  
  ndvi = NDVI_class_end,
  
  habitat_distance = log_habitat,
  habitat_distance_s1 = log_habitat * hour_s1,
  habitat_distance_c1 = log_habitat * hour_c1,
  habitat_distance_s2 = log_habitat * hour_s2, # New
  habitat_distance_c2 = log_habitat * hour_c2, # New
  
  disturb_distance = log_disturb,
  disturb_distance_s1 = log_disturb * hour_s1,
  disturb_distance_c1 = log_disturb * hour_c1,
  disturb_distance_s2 = log_disturb * hour_s2, # New
  disturb_distance_c2 = log_disturb * hour_c2, # New
  
  memory = memory_log,
  memory_s1 = memory_log * hour_s1,
  memory_c1 = memory_log * hour_c1,
  memory_s2 = memory_log * hour_s2, # New
  memory_c2 = memory_log * hour_c2, # New
  
  #spatial_memory_pi = spatial_memory_pi ^ 2,
  #spatial_memory_pi_s1 = (spatial_memory_pi ^ 2) * hour_s1,
  #spatial_memory_pi_c1 = (spatial_memory_pi ^ 2) * hour_c1,
  #spatial_memory_pi_s2 = (spatial_memory_pi ^ 2) * hour_s2, # New
  #spatial_memory_pi_c2 = (spatial_memory_pi ^ 2) * hour_c2, # New
  
  step_l = sl,
  # step_l_s1 = sl * hour_s1,
  # step_l_c1 = sl * hour_c1,
  # step_l_s2 = sl * hour_s2, # New
  # step_l_c2 = sl * hour_c2, # New
  
  log_step_l = log_sl,
  # log_step_l_s1 = log_sl * hour_s1,
  # log_step_l_c1 = log_sl * hour_c1,
  # log_step_l_s2 = log_sl * hour_s2, # New
  # log_step_l_c2 = log_sl * hour_c2, # New
  # 
  cos_turn_a = cos_ta,
  # cos_turn_a_s1 = cos_ta * hour_s1,
  # cos_turn_a_c1 = cos_ta * hour_c1,
  # cos_turn_a_s2 = cos_ta * hour_s2, # New
  # cos_turn_a_c2 = cos_ta * hour_c2 # New
)


head(quoll_data_matrix_unscaled)
quoll_data_matrix_unscaled<-as.data.frame(quoll_data_matrix_unscaled)
## Identify numeric columns
#numeric_columns <- sapply(quoll_data_matrix_unscaled, is.numeric)
## Scale only numeric columns
#scaled_data <- scale(quoll_data_matrix_unscaled[, numeric_columns])
## Combine scaled numeric data with non-numeric columns
#quoll_data_matrix_scaled <- cbind(quoll_data_matrix_unscaled[!numeric_columns], scaled_data)
## Convert row names to the first column if necessary
#quoll_data_matrix_scaled <- data.frame(quoll_data_matrix_scaled, row.names = NULL)
## View the head of the scaled dataset
#head(quoll_data_matrix_scaled)

# Using base R
quoll_data_matrix_unscaled_no_ndvi <- quoll_data_matrix_unscaled[ , !(names(quoll_data_matrix_unscaled) %in% c("ndvi","step_l","log_step_l","cos_turn_a"))]
quoll_data_matrix_scaled_no_ndvi <- scale(quoll_data_matrix_unscaled_no_ndvi)

quoll_data_matrix_scaled <- cbind(ndvi = quoll_data_matrix_unscaled$ndvi, quoll_data_matrix_scaled_no_ndvi)


mean_vals <- attr(quoll_data_matrix_scaled_no_ndvi, "scaled:center")
sd_vals <- attr(quoll_data_matrix_scaled_no_ndvi, "scaled:scale")
scaling_attributes <- data.frame(variable = names(quoll_data_matrix_unscaled_no_ndvi), mean = mean_vals, sd = sd_vals)

quoll_data_scaled_2p <- data.frame(id = all_ssf_breeding$id,  step_id = all_ssf_breeding$step_id, case_ = all_ssf_breeding$case_, quoll_data_matrix_scaled)
head(quoll_data_scaled_2p)

# Ensure the row order and number of rows in both datasets match
# Replace the ndvi column in quoll_data_scaled_1p with that from quoll_data_matrix_unscaled
quoll_data_scaled_2p$ndvi <- quoll_data_matrix_unscaled$ndvi
quoll_data_scaled_2p$step_l <- quoll_data_matrix_unscaled$step_l
quoll_data_scaled_2p$log_step_l <- quoll_data_matrix_unscaled$log_step_l
quoll_data_scaled_2p$cos_turn_a <- quoll_data_matrix_unscaled$cos_turn_a


# View the first few rows to confirm the change
head(quoll_data_scaled_2p)



formula_twostep <- case_ ~ 
  
  ndvi +
  
  habitat_distance +
  habitat_distance_s1 +
  habitat_distance_c1 +
  habitat_distance_s2 + # Added
  habitat_distance_c2 + # Added
  
  disturb_distance +
  disturb_distance_s1 +
  disturb_distance_c1 +
  disturb_distance_s2 + # Added
  disturb_distance_c2 + # Added
  
  memory +
  memory_s1 +
  memory_c1 +
  memory_s2 + # Added
  memory_c2 + # Added
  
  #spatial_memory_pi +
  #spatial_memory_pi_s1 +
  #spatial_memory_pi_c1 +
  #spatial_memory_pi_s2 + # Added
  #spatial_memory_pi_c2 + # Added
  
  step_l +
  # step_l_s1 +
  # step_l_c1 +
  # step_l_s2 + # Added
  # step_l_c2 + # Added
  
  log_step_l +
  # log_step_l_s1 +
  # log_step_l_c1 +
  # log_step_l_s2 + # Added
  # log_step_l_c2 + # Added
  
  cos_turn_a +
  # cos_turn_a_s1 +
  # cos_turn_a_c1 +
  # cos_turn_a_s2 + # Added
  # cos_turn_a_c2 + # Added
  
  strata(step_id) +
  cluster(id)


library(tictoc)
library(TwoStepCLogit)

  tic()
  model_twostep_2p_harms <- Ts.estim(formula = formula_twostep,
                                     data = quoll_data_scaled_2p,
                                     all.m.1 = TRUE,
                                     D = "UN(1)",
                                     itermax = 10000)
  toc()
# 

#summary(model_twostep_2p_harms)
model_twostep_2p_harms
model_twostep_2p_harms$beta
model_twostep_2p_harms$se
model_twostep_2p_harms$vcov
diag(model_twostep_2p_harms$D) # between cluster variance
model_twostep_2p_harms$r.effect # individual estimates

hist(model_twostep_2p_harms$r.effect[,6])


coefs_clr <- data.frame(coefs = names(model_twostep_2p_harms$beta), value = model_twostep_2p_harms$beta)
coefs_clr_no_ndvi <- coefs_clr[!coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), ]

coefs_clr_no_ndvi$scale_sd <- scaling_attributes$sd
coefs_clr_no_ndvi <- coefs_clr_no_ndvi %>% mutate(value_nat = value / scale_sd)

head(coefs_clr)
head(coefs_clr_no_ndvi)

# Rows to add back
rows_to_add <- coefs_clr[coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), ]
rows_to_add$scale_sd <- NA
rows_to_add$value_nat <- coefs_clr[coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), "value"]
coefs_clr_upd <- rbind(coefs_clr_no_ndvi, rows_to_add)
print(coefs_clr_upd)

####Reconstructing coefficients with two pairs of harmonics, with quadratic terms

# we want to reconstruct the temporally dynamic coefficients for the full cycle, and then we can subset to the active hours after
hour <- seq(0.0,23.5,0.5)

all_ssf_breeding$hour30 <- format(as.POSIXct(all_ssf_breeding$t2), format = "%H:%M")
all_ssf_breeding$hour30 <- gsub(":", ".", all_ssf_breeding$hour30)
all_ssf_breeding$hour30 <- sub("30", "50", all_ssf_breeding$hour30)
all_ssf_breeding$hour30 <- as.numeric(all_ssf_breeding$hour30)

unique_hours<-unique(all_ssf_breeding$hour30)


# this creates a matrix of the hour harmonics, which is 24 x 3 (hours x harmonics, including the linear term)
hour_harmonics_df <- data.frame("linear_term" = rep(1, length(hour)),
                                "hour_s1" = sin(2*pi*hour/24),
                                "hour_s2" = sin(4*pi*hour/24),
                                "hour_c1" = cos(2*pi*hour/24),
                                "hour_c2" = cos(4*pi*hour/24))

hour_harmonics_linear_df <- data.frame("linear_term" = rep(1, length(hour)))


# we'll do the filtering after reconstructing the harmonics - i think it's best to do the whole cycle initially
# hour_harmonics_df <- hour_harmonics_df[hour_harmonics_df$hour %in% unique_hours, ]

hour_transformed = ifelse(hour < 6, hour + 24, hour)

# i spaced the coefficients out so i could see them a bit easier
harmonics_scaled_df_2p <- data.frame(
  "hour" = hour,
  # can add this in here now
  "hour_transformed" = ifelse(hour < 6, hour + 24, hour),
  
  # note that for the linear only terms i've added in the LINEAR ONLY data frame here, as there is only the linear term, and the matrix multiplication will work (now it's just an 1 x 1 %*% 1 x n, which will result in a 1 x n matrix (i.e. a vector))
  "ndvigrassland" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvigrassland", coefs)) %>%
                                 ###### under here ######
                               pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvidense_veg" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvidense_veg", coefs)) %>%
                                 pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvimine_pit_waste_dump" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvimine_pit_waste_dump", coefs)) %>%
                                           pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndviother_disturbed" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndviother_disturbed", coefs)) %>%
                                       pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  # For these ones, we will pull out the three coefficients for each variable (linear, sin1,cos1), which R considers to be COLUMN, and therefore a 3 x 1 matrix (or more simply a (column) vector). We transpose (the t() function) it to make it a 1 x 3. We then have the 24 x 3 matrix of the harmonics, which we TRANSPOSE to turn it into a 3 x 24 matrix. Then when we multiply them together, we have 1 x 3 %*% 3 x 24, which results in a 1 x 24 matrix (i.e. a VECTOR of the temporally dynamic coefficient throughout the day), which we add into the dataframe
  
  "step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("step_l", coefs) & !grepl("log", coefs)) %>% 
                          pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  # add in the log_sl term
  "log_step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("log_step_l", coefs)) %>% 
                              pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "cos_turn_a" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("cos_turn_a", coefs)) %>% 
                              pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "disturb_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("disturb_distance", coefs) & !grepl("sq", coefs)) %>%
                                    pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "habitat_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("habitat_distance", coefs)) %>% 
                                    pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "memory" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("memory", coefs)) %>% 
                                      pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
  #,
  
  #"spatial_memory_pi" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("spatial_memory_pi", coefs)) %>% 
  #                                   pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
  )



harmonics_nat_df_2p <- data.frame(
  "hour" = hour,
  # can add this in here now
  "hour_transformed" = ifelse(hour < 6, hour + 24, hour),
  
  # note that for the linear only terms i've added in the LINEAR ONLY data frame here, as there is only the linear term, and the matrix multiplication will work (now it's just an 1 x 1 %*% 1 x n, which will result in a 1 x n matrix (i.e. a vector))
  "ndvigrassland" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvigrassland", coefs)) %>%
                                 ###### under here ######
                               pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvidense_veg" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvidense_veg", coefs)) %>%
                                 pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvimine_pit_waste_dump" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvimine_pit_waste_dump", coefs)) %>%
                                           pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndviother_disturbed" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndviother_disturbed", coefs)) %>%
                                       pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  # For these ones, we will pull out the three coefficients for each variable (linear, sin1,cos1), which R considers to be COLUMN, and therefore a 3 x 1 matrix (or more simply a (column) vector). We transpose (the t() function) it to make it a 1 x 3. We then have the 24 x 3 matrix of the harmonics, which we TRANSPOSE to turn it into a 3 x 24 matrix. Then when we multiply them together, we have 1 x 3 %*% 3 x 24, which results in a 1 x 24 matrix (i.e. a VECTOR of the temporally dynamic coefficient throughout the day), which we add into the dataframe
  
  "step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("step_l", coefs) & !grepl("log", coefs)) %>% 
                          pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  # add in the log_sl term
  "log_step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("log_step_l", coefs)) %>%
                                pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "cos_turn_a" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("cos_turn_a", coefs)) %>%
                              pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "disturb_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("disturb_distance", coefs) & !grepl("sq", coefs)) %>%
                                    pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "habitat_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("habitat_distance", coefs)) %>% 
                                    pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "memory" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("memory", coefs)) %>% 
                                      pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df)))
  #,
  
  #"spatial_memory_pi" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("spatial_memory_pi", coefs)) %>% 
  #                                   pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
)




##Get the tentative distributions
all_ssf_breed_true <- all_ssf_breeding %>% 
  filter(case_ == "1") %>% 
  arrange(t1_)
# Fit a gamma distribution to the step length of true steps
fit_sl <- fit_distr(all_ssf_breed_true$sl_, "gamma")
# Extract the shape and rate parameters of the gamma distribution
shape_sl <- fit_sl$params[1]
scale_sl <- fit_sl$params[2]
# Load the circular package
library(circular)
# Fit a von Mises distribution to the ta_ values
fit_ta <- fit_distr(all_ssf_breed_true$ta_, "vonmises")
# Extract the kappa parameter of the von Mises distribution
kappa_ta <- fit_ta$params["kappa"]


tentative_shape <- shape_sl$shape
tentative_scale <- scale_sl$scale
tentative_kappa <- kappa_ta$kappa

hist(rgamma(1e4,
       shape = tentative_shape,
       scale = tentative_scale),
       breaks = 100)

harmonics_nat_scaled_df_2p <- harmonics_nat_df_2p %>% mutate(shape = tentative_shape + log_step_l,
                                                                 scale = 1/((1/tentative_scale) - step_l),
                                                                 kappa = tentative_kappa + cos_turn_a)

hist(rgamma(1e4,
            shape = harmonics_nat_scaled_df_2p$shape[1],
            scale = harmonics_nat_scaled_df_2p$scale[1]),
     breaks = 100)




harmonics_nat_scaled_df_2p
# Add a new column 'rocky' with 0s to the dataframe 'harmonics_nat_scaled_df_2p'
harmonics_nat_scaled_df_2p$ndvirocky <- 0

harmonics_scaled_long_2p <- pivot_longer(harmonics_nat_scaled_df_2p, cols = !1:2, names_to = "coef")

ggplot() +
  geom_path(data = harmonics_scaled_long_2p,
            aes(x = hour, y = value, colour = coef)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour") +
  scale_color_discrete("Estimate") +
  theme_classic() +
  theme(legend.position = "bottom")


# when removing the daytime hours
harmonics_scaled_df_2p_active <- harmonics_nat_scaled_df_2p %>% dplyr::filter(hour %in% unique_hours)
harmonics_scaled_long_2p_active <- pivot_longer(harmonics_scaled_df_2p_active, cols = !1:2, names_to = "coef")
##write.csv(harmonics_scaled_df_2p_active,"C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_3hrm.csv")
movement_summary <- all_ssf_breeding %>% filter(case_ == 1) %>%  group_by(hour30) %>% summarise(mean_sl = mean(sl), median_sl = median(sl))
movement_summary$hour <- movement_summary$hour30

final2<-cbind(harmonics_scaled_df_2p_active,movement_summary)
final2<-final2[1:20]
final2<-final2 %>% dplyr::select(-hour30)

#write.csv(final2,"C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_2hrm.csv",row.names = FALSE)


harmonics_scaled_df_2p_active<-read.csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_2hrm.csv")


library(ggplot2)
# Generate a set of random colors
# First, find the number of unique levels in the 'coef' column
num_levels <- length(unique(harmonics_scaled_long_2p_active$coef))

# Then, generate random colors for each level
# You can adjust the seed for reproducibility
set.seed(123) # Remove or change the seed value for different colors
random_colors <- grDevices::rainbow(num_levels)

# Now, create the plot
ggplot(harmonics_scaled_long_2p_active, aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(size = 2.5) +  # Increased line thickness
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5","6"),
                     limits = c(19, 30)) +
  scale_color_manual("Estimate", values = setNames(random_colors, unique(harmonics_scaled_long_2p_active$coef))) +
  theme_classic() +
  theme(legend.position = "bottom")



ggplot() +
  geom_path(data = harmonics_scaled_long_2p %>%
              filter(!coef %in% c("step_l", "log_step_l", "cos_turn_a", "memory")),
            aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5","6"),
                     limits = c(19, 30)) +
  scale_color_manual("Estimate", values = setNames(random_colors, unique(harmonics_scaled_long_2p_active$coef))) +
  theme_classic() +
  theme(legend.position = "bottom")




ggplot() +
  geom_path(data = harmonics_scaled_long_2p %>%
              filter(coef %in% c("step_l", "log_step_l", "cos_turn_a")),
            aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5","6"),
                     limits = c(19, 30)) +
  scale_color_manual("Estimate", values = setNames(random_colors, unique(harmonics_scaled_long_2p_active$coef))) +
  theme_classic() +
  theme(legend.position = "bottom")









########################
#Three pairs of harmonics
########################
library(amt)
library(MuMIn)
library(gdata)
library(survival)
library(ggplot2)
library(tidyverse)
library(AICcmodavg)
library(lubridate)

#dev.off()

setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")
#all_ssf1 <-read.csv("Woodie_SSF_w_perc.csv")
all_ssf1 <-read.csv("Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")
names(all_ssf1)

# Assuming all_ssf1 is your dataframe and t1_ is the column with timestamps
all_ssf1$t1_ <- as.POSIXct(all_ssf1$t1_, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")
all_ssf1$t2_ <- as.POSIXct(all_ssf1$t2_, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")

# Load the lubridate package
library(lubridate)
all_ssf1$realtime <- all_ssf1$t2_
# Extract the hour
hour_component <- hour(all_ssf1$t2_)
# Extract the minute and determine if .5 should be added
minute_component <- ifelse(minute(all_ssf1$t2_) == 30, 0.5, 0)
# Combine hour and minute components
all_ssf1$hour_from_t2_ <- hour_component + minute_component
# View the modified dataframe
head(all_ssf1)
# 
# # Define the original hours in order
# original_hours <- c(18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
# # Create a sequence from 0 to 11.5 in steps of 0.5
# new_sequence <- seq(1, 24, by=1)
# # Create a mapping from original hours to new sequence
# hour_mapping <- setNames(new_sequence, original_hours)
# # Assuming hourly_coefs is your dataframe
# # Map the original hour values to the new sequence
# all_ssf1$newhour <- hour_mapping[as.character(all_ssf1$hour_from_t2_)]
constant=0.1
#all_ssf1$t1_<-dmy_hm(all_ssf1$t1_)
#all_ssf1$t2_<-dmy_hm(all_ssf1$t2_)
all_ssf2 <- all_ssf1 %>%
  mutate(id_num = as.numeric(factor(id)),
         step_id = step_id_,
         burst_ = burst_,
         x1 = x1_, x2 = x2_,
         y1 = y1_, y2 = y2_,
         t1 = t1_,
         t1_rounded = round_date(all_ssf1$t1_, "hour"),
         hour_t1 = hour(t1_rounded),
         t2 = t2_,
         t2_rounded = round_date(all_ssf1$t2_, "hour"),
         hour_t2 = hour(t2_rounded),
         hour = hour_t2,
         yday = yday(t1_),
         year = year(t1_),
         month = month(t1_),
         sl = sl_,
         log_sl = log(sl_),
         ta = ta_,
         cos_ta = cos(ta_),
         NDVI_class_end = NDVI_class_end,
         disturb_distance_end = disturb_distance_end,
         log_disturb = log(disturb_distance_end+constant),
         log_habitat = log(habitat_distance_end+constant),
         habitat_distance_end = habitat_distance_end,
         name=name,orig_name=orig_name,season=season,age=age,
         spatial_memory_ref = kde_ref_spatial_memory_density_log,
         spatial_memory_pi = kde_pi_spatial_memory_density_log,
         memory_log=kde_ref_spatial_memory_density_log,
         memory=exp(kde_ref_spatial_memory_density_log),
         hour_s1 = sin(2*pi*hour/24),
         hour_s2 = sin(4*pi*hour/24),
         hour_s3 = sin(6*pi*hour/24),
         hour_s4 = sin(8*pi*hour/24),
         hour_c1 = cos(2*pi*hour/24),
         hour_c2 = cos(4*pi*hour/24),
         hour_c3 = cos(6*pi*hour/24),
         hour_c4 = cos(8*pi*hour/24))


##Look at habitat distance through time (days)
ggplot(data = all_ssf2 %>% filter(case_ == 1,season=="breeding"),
       aes(x = t1_, y = habitat_distance_end, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic()

##Look at habitat distance through time (days)
ggplot(data = all_ssf2 %>% filter(case_ == 1,season=="breeding"),
       aes(x = t1_, y = memory, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic()

##Look at habitat distance through time (hours)
all_ssf2 <- all_ssf2 %>%
  mutate(hour_transformed = ifelse(hour < 6, hour + 24, hour))

ggplot(data = all_ssf2 %>% filter(case_ == 1, season == "breeding"),
       aes(x = hour_transformed, y = habitat_distance_end, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic() +
  coord_cartesian(xlim = c(18, 30)) # Adjusted based on the transformed hours


ggplot(data = all_ssf2 %>% filter(case_ == 1, season == "breeding"),
       aes(x = hour_transformed, y = memory, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic() +
  coord_cartesian(xlim = c(18, 30)) # Adjusted based on the transformed hours


##Separate seasons for separate modelling
################For breeding season
all_ssf_breeding <- all_ssf2 %>% 
  filter(season == "breeding") %>% 
  arrange(t1_)

all_ssf_breeding$NDVI_class_end <-factor(all_ssf_breeding$NDVI_class_end, levels = c("rocky", "grassland", "dense_veg","mine_pit_waste_dump","other_disturbed"))

unique(all_ssf_breeding$NDVI_class_end)


quoll_ids <- unique(all_ssf_breeding$id)


quoll_data_matrix_unscaled <- all_ssf_breeding %>% transmute(
  
  ndvi = NDVI_class_end,
  
  habitat_distance = log_habitat,
  habitat_distance_s1 = log_habitat * hour_s1,
  habitat_distance_c1 = log_habitat * hour_c1,
  habitat_distance_s2 = log_habitat * hour_s2, # New
  habitat_distance_c2 = log_habitat * hour_c2,
  habitat_distance_s3 = log_habitat * hour_s3, # New
  habitat_distance_c3 = log_habitat * hour_c3,# New
  
  disturb_distance = log_disturb,
  disturb_distance_s1 = log_disturb * hour_s1,
  disturb_distance_c1 = log_disturb * hour_c1,
  disturb_distance_s2 = log_disturb * hour_s2, # New
  disturb_distance_c2 = log_disturb * hour_c2,
  disturb_distance_s3 = log_disturb * hour_s3, # New
  disturb_distance_c3 = log_disturb * hour_c3,# New
  
  memory = memory_log,
  memory_s1 = memory_log * hour_s1,
  memory_c1 = memory_log * hour_c1,
  memory_s2 = memory_log * hour_s2, # New
  memory_c2 = memory_log * hour_c2,
  memory_s3 = memory_log * hour_s3, # New
  memory_c3 = memory_log * hour_c3,# New
  
  #spatial_memory_pi = spatial_memory_pi ^ 2,
  #spatial_memory_pi_s1 = (spatial_memory_pi ^ 2) * hour_s1,
  #spatial_memory_pi_c1 = (spatial_memory_pi ^ 2) * hour_c1,
  #spatial_memory_pi_s2 = (spatial_memory_pi ^ 2) * hour_s2, # New
  #spatial_memory_pi_c2 = (spatial_memory_pi ^ 2) * hour_c2, # New
  
  step_l = sl,
  # step_l_s1 = sl * hour_s1,
  # step_l_c1 = sl * hour_c1,
  # step_l_s2 = sl * hour_s2, # New
  # step_l_c2 = sl * hour_c2, # New
  
  log_step_l = log_sl,
  # log_step_l_s1 = log_sl * hour_s1,
  # log_step_l_c1 = log_sl * hour_c1,
  # log_step_l_s2 = log_sl * hour_s2, # New
  # log_step_l_c2 = log_sl * hour_c2, # New
  # 
  cos_turn_a = cos_ta,
  # cos_turn_a_s1 = cos_ta * hour_s1,
  # cos_turn_a_c1 = cos_ta * hour_c1,
  # cos_turn_a_s2 = cos_ta * hour_s2, # New
  # cos_turn_a_c2 = cos_ta * hour_c2 # New
)


head(quoll_data_matrix_unscaled)
quoll_data_matrix_unscaled<-as.data.frame(quoll_data_matrix_unscaled)
## Identify numeric columns
#numeric_columns <- sapply(quoll_data_matrix_unscaled, is.numeric)
## Scale only numeric columns
#scaled_data <- scale(quoll_data_matrix_unscaled[, numeric_columns])
## Combine scaled numeric data with non-numeric columns
#quoll_data_matrix_scaled <- cbind(quoll_data_matrix_unscaled[!numeric_columns], scaled_data)
## Convert row names to the first column if necessary
#quoll_data_matrix_scaled <- data.frame(quoll_data_matrix_scaled, row.names = NULL)
## View the head of the scaled dataset
#head(quoll_data_matrix_scaled)

# Using base R
quoll_data_matrix_unscaled_no_ndvi <- quoll_data_matrix_unscaled[ , !(names(quoll_data_matrix_unscaled) %in% c("ndvi","step_l","log_step_l","cos_turn_a"))]
quoll_data_matrix_scaled_no_ndvi <- scale(quoll_data_matrix_unscaled_no_ndvi)

quoll_data_matrix_scaled <- cbind(ndvi = quoll_data_matrix_unscaled$ndvi, quoll_data_matrix_scaled_no_ndvi)


mean_vals <- attr(quoll_data_matrix_scaled_no_ndvi, "scaled:center")
sd_vals <- attr(quoll_data_matrix_scaled_no_ndvi, "scaled:scale")
scaling_attributes <- data.frame(variable = names(quoll_data_matrix_unscaled_no_ndvi), mean = mean_vals, sd = sd_vals)

quoll_data_scaled_3p <- data.frame(id = all_ssf_breeding$id,  step_id = all_ssf_breeding$step_id, case_ = all_ssf_breeding$case_, quoll_data_matrix_scaled)
head(quoll_data_scaled_3p)

# Ensure the row order and number of rows in both datasets match
# Replace the ndvi column in quoll_data_scaled_1p with that from quoll_data_matrix_unscaled
quoll_data_scaled_3p$ndvi <- quoll_data_matrix_unscaled$ndvi
quoll_data_scaled_3p$step_l <- quoll_data_matrix_unscaled$step_l
quoll_data_scaled_3p$log_step_l <- quoll_data_matrix_unscaled$log_step_l
quoll_data_scaled_3p$cos_turn_a <- quoll_data_matrix_unscaled$cos_turn_a


# View the first few rows to confirm the change
head(quoll_data_scaled_3p)



formula_twostep <- case_ ~ 
  
  ndvi +
  
  habitat_distance +
  habitat_distance_s1 +
  habitat_distance_c1 +
  habitat_distance_s2 + # Added
  habitat_distance_c2 +
  habitat_distance_s3 + # Added
  habitat_distance_c3 +# Added
  
  disturb_distance +
  disturb_distance_s1 +
  disturb_distance_c1 +
  disturb_distance_s2 + # Added
  disturb_distance_c2 +
  disturb_distance_s3 + # Added
  disturb_distance_c3 +# Added
  
  memory +
  memory_s1 +
  memory_c1 +
  memory_s2 + # Added
  memory_c2 +
  memory_s3 + # Added
  memory_c3 +# Added
  
  #spatial_memory_pi +
  #spatial_memory_pi_s1 +
  #spatial_memory_pi_c1 +
  #spatial_memory_pi_s2 + # Added
  #spatial_memory_pi_c2 + # Added
  
  step_l +
  # step_l_s1 +
  # step_l_c1 +
  # step_l_s2 + # Added
  # step_l_c2 + # Added
  
  log_step_l +
  # log_step_l_s1 +
  # log_step_l_c1 +
  # log_step_l_s2 + # Added
  # log_step_l_c2 + # Added
  
  cos_turn_a +
  # cos_turn_a_s1 +
  # cos_turn_a_c1 +
  # cos_turn_a_s2 + # Added
  # cos_turn_a_c2 + # Added
  
  strata(step_id) +
  cluster(id)


library(tictoc)
library(TwoStepCLogit)

  tic()
  model_twostep_3p_harms <- Ts.estim(formula = formula_twostep,
                                     data = quoll_data_scaled_3p,
                                     all.m.1 = TRUE,
                                     D = "UN(1)",
                                     itermax = 10000)
  toc()


#summary(model_twostep_3p_harms)
model_twostep_3p_harms
model_twostep_3p_harms$beta
model_twostep_3p_harms$se
model_twostep_3p_harms$vcov
diag(model_twostep_3p_harms$D) # between cluster variance
model_twostep_3p_harms$r.effect # individual estimates

hist(model_twostep_3p_harms$r.effect[,6])


coefs_clr <- data.frame(coefs = names(model_twostep_3p_harms$beta), value = model_twostep_3p_harms$beta)
coefs_clr_no_ndvi <- coefs_clr[!coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), ]

coefs_clr_no_ndvi$scale_sd <- scaling_attributes$sd
coefs_clr_no_ndvi <- coefs_clr_no_ndvi %>% mutate(value_nat = value / scale_sd)

head(coefs_clr)
head(coefs_clr_no_ndvi)

# Rows to add back
rows_to_add <- coefs_clr[coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), ]
rows_to_add$scale_sd <- NA
rows_to_add$value_nat <- coefs_clr[coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), "value"]
coefs_clr_upd <- rbind(coefs_clr_no_ndvi, rows_to_add)
print(coefs_clr_upd)

####Reconstructing coefficients with two pairs of harmonics, with quadratic terms

# we want to reconstruct the temporally dynamic coefficients for the full cycle, and then we can subset to the active hours after
hour <- seq(0.0,23.5,0.5)

all_ssf_breeding$hour30 <- format(as.POSIXct(all_ssf_breeding$t2), format = "%H:%M")
all_ssf_breeding$hour30 <- gsub(":", ".", all_ssf_breeding$hour30)
all_ssf_breeding$hour30 <- sub("30", "50", all_ssf_breeding$hour30)
all_ssf_breeding$hour30 <- as.numeric(all_ssf_breeding$hour30)

unique_hours<-unique(all_ssf_breeding$hour30)


# this creates a matrix of the hour harmonics, which is 24 x 3 (hours x harmonics, including the linear term)
hour_harmonics_df <- data.frame("linear_term" = rep(1, length(hour)),
                                "hour_s1" = sin(2*pi*hour/24),
                                "hour_s2" = sin(4*pi*hour/24),
                                "hour_s3" = sin(2*pi*hour/24),
                                "hour_c1" = cos(2*pi*hour/24),
                                "hour_c2" = cos(4*pi*hour/24),
                                "hour_c3" = cos(4*pi*hour/24))

hour_harmonics_linear_df <- data.frame("linear_term" = rep(1, length(hour)))


# we'll do the filtering after reconstructing the harmonics - i think it's best to do the whole cycle initially
# hour_harmonics_df <- hour_harmonics_df[hour_harmonics_df$hour %in% unique_hours, ]

hour_transformed = ifelse(hour < 6, hour + 24, hour)

# i spaced the coefficients out so i could see them a bit easier
harmonics_scaled_df_3p <- data.frame(
  "hour" = hour,
  # can add this in here now
  "hour_transformed" = ifelse(hour < 6, hour + 24, hour),
  
  # note that for the linear only terms i've added in the LINEAR ONLY data frame here, as there is only the linear term, and the matrix multiplication will work (now it's just an 1 x 1 %*% 1 x n, which will result in a 1 x n matrix (i.e. a vector))
  "ndvigrassland" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvigrassland", coefs)) %>%
                                 ###### under here ######
                               pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvidense_veg" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvidense_veg", coefs)) %>%
                                 pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvimine_pit_waste_dump" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvimine_pit_waste_dump", coefs)) %>%
                                           pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndviother_disturbed" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndviother_disturbed", coefs)) %>%
                                       pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  # For these ones, we will pull out the three coefficients for each variable (linear, sin1,cos1), which R considers to be COLUMN, and therefore a 3 x 1 matrix (or more simply a (column) vector). We transpose (the t() function) it to make it a 1 x 3. We then have the 24 x 3 matrix of the harmonics, which we TRANSPOSE to turn it into a 3 x 24 matrix. Then when we multiply them together, we have 1 x 3 %*% 3 x 24, which results in a 1 x 24 matrix (i.e. a VECTOR of the temporally dynamic coefficient throughout the day), which we add into the dataframe
  
  "step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("step_l", coefs) & !grepl("log", coefs)) %>% 
                          pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  # add in the log_sl term
  "log_step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("log_step_l", coefs)) %>% 
                              pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "cos_turn_a" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("cos_turn_a", coefs)) %>% 
                              pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "disturb_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("disturb_distance", coefs) & !grepl("sq", coefs)) %>%
                                    pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "habitat_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("habitat_distance", coefs)) %>% 
                                    pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "memory" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("memory", coefs)) %>% 
                          pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
  #,
  
  #"spatial_memory_pi" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("spatial_memory_pi", coefs)) %>% 
  #                                   pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
)



harmonics_nat_df_3p <- data.frame(
  "hour" = hour,
  # can add this in here now
  "hour_transformed" = ifelse(hour < 6, hour + 24, hour),
  
  # note that for the linear only terms i've added in the LINEAR ONLY data frame here, as there is only the linear term, and the matrix multiplication will work (now it's just an 1 x 1 %*% 1 x n, which will result in a 1 x n matrix (i.e. a vector))
  "ndvigrassland" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvigrassland", coefs)) %>%
                                 ###### under here ######
                               pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvidense_veg" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvidense_veg", coefs)) %>%
                                 pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvimine_pit_waste_dump" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvimine_pit_waste_dump", coefs)) %>%
                                           pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndviother_disturbed" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndviother_disturbed", coefs)) %>%
                                       pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  # For these ones, we will pull out the three coefficients for each variable (linear, sin1,cos1), which R considers to be COLUMN, and therefore a 3 x 1 matrix (or more simply a (column) vector). We transpose (the t() function) it to make it a 1 x 3. We then have the 24 x 3 matrix of the harmonics, which we TRANSPOSE to turn it into a 3 x 24 matrix. Then when we multiply them together, we have 1 x 3 %*% 3 x 24, which results in a 1 x 24 matrix (i.e. a VECTOR of the temporally dynamic coefficient throughout the day), which we add into the dataframe
  
  "step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("step_l", coefs) & !grepl("log", coefs)) %>% 
                          pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  # add in the log_sl term
  "log_step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("log_step_l", coefs)) %>%
                              pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "cos_turn_a" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("cos_turn_a", coefs)) %>%
                              pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "disturb_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("disturb_distance", coefs) & !grepl("sq", coefs)) %>%
                                    pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "habitat_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("habitat_distance", coefs)) %>% 
                                    pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "memory" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("memory", coefs)) %>% 
                          pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df)))
  #,
  
  #"spatial_memory_pi" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("spatial_memory_pi", coefs)) %>% 
  #                                   pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
)




##Get the tentative distributions
all_ssf_breed_true <- all_ssf_breeding %>% 
  filter(case_ == "1") %>% 
  arrange(t1_)
# Fit a gamma distribution to the step length of true steps
fit_sl <- fit_distr(all_ssf_breed_true$sl_, "gamma")
# Extract the shape and rate parameters of the gamma distribution
shape_sl <- fit_sl$params[1]
scale_sl <- fit_sl$params[2]
# Load the circular package
library(circular)
# Fit a von Mises distribution to the ta_ values
fit_ta <- fit_distr(all_ssf_breed_true$ta_, "vonmises")
# Extract the kappa parameter of the von Mises distribution
kappa_ta <- fit_ta$params["kappa"]


tentative_shape <- shape_sl$shape
tentative_scale <- scale_sl$scale
tentative_kappa <- kappa_ta$kappa

hist(rgamma(1e4,
            shape = tentative_shape,
            scale = tentative_scale),
     breaks = 100)

harmonics_nat_scaled_df_3p <- harmonics_nat_df_3p %>% mutate(shape = tentative_shape + log_step_l,
                                                             scale = 1/((1/tentative_scale) - step_l),
                                                             kappa = tentative_kappa + cos_turn_a)

hist(rgamma(1e4,
            shape = harmonics_nat_scaled_df_3p$shape[1],
            scale = harmonics_nat_scaled_df_3p$scale[1]),
     breaks = 100)




harmonics_nat_scaled_df_3p
# Add a new column 'rocky' with 0s to the dataframe 'harmonics_nat_scaled_df_3p'
harmonics_nat_scaled_df_3p$ndvirocky <- 0

harmonics_scaled_long_3p <- pivot_longer(harmonics_nat_scaled_df_3p, cols = !1:2, names_to = "coef")

ggplot() +
  geom_path(data = harmonics_scaled_long_3p,
            aes(x = hour, y = value, colour = coef)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour") +
  scale_color_discrete("Estimate") +
  theme_classic() +
  theme(legend.position = "bottom")


# when removing the daytime hours
harmonics_scaled_df_3p_active <- harmonics_nat_scaled_df_3p %>% dplyr::filter(hour %in% unique_hours)
harmonics_scaled_long_3p_active <- pivot_longer(harmonics_scaled_df_3p_active, cols = !1:2, names_to = "coef")
##write.csv(harmonics_scaled_df_3p_active,"C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_3hrm.csv")
movement_summary <- all_ssf_breeding %>% filter(case_ == 1) %>%  group_by(hour30) %>% summarise(mean_sl = mean(sl), median_sl = median(sl))
movement_summary$hour <- movement_summary$hour30

final3<-cbind(harmonics_scaled_df_3p_active,movement_summary)
final3<-final3[1:20]
final3<-final3 %>% dplyr::select(-hour30)

#write.csv(final3,"C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_3hrm.csv",row.names = FALSE)


harmonics_scaled_df_3p_active<-read.csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_3hrm.csv")


library(ggplot2)
# Generate a set of random colors
# First, find the number of unique levels in the 'coef' column
num_levels <- length(unique(harmonics_scaled_long_3p_active$coef))

# Then, generate random colors for each level
# You can adjust the seed for reproducibility
set.seed(123) # Remove or change the seed value for different colors
random_colors <- grDevices::rainbow(num_levels)

# Now, create the plot
ggplot(harmonics_scaled_long_3p_active, aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(size = 2.5) +  # Increased line thickness
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5","6"),
                     limits = c(19, 30)) +
  scale_color_manual("Estimate", values = setNames(random_colors, unique(harmonics_scaled_long_3p_active$coef))) +
  theme_classic() +
  theme(legend.position = "bottom")



ggplot() +
  geom_path(data = harmonics_scaled_long_3p %>%
              filter(!coef %in% c("step_l", "log_step_l", "cos_turn_a", "memory")),
            aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5","6"),
                     limits = c(19, 30)) +
  scale_color_manual("Estimate", values = setNames(random_colors, unique(harmonics_scaled_long_3p_active$coef))) +
  theme_classic() +
  theme(legend.position = "bottom")




ggplot() +
  geom_path(data = harmonics_scaled_long_3p %>%
              filter(coef %in% c("step_l", "log_step_l", "cos_turn_a")),
            aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5","6"),
                     limits = c(19, 30)) +
  scale_color_manual("Estimate", values = setNames(random_colors, unique(harmonics_scaled_long_3p_active$coef))) +
  theme_classic() +
  theme(legend.position = "bottom")





##############################
########## 0 harmonics #######
##############################
library(amt)
library(MuMIn)
library(gdata)
library(survival)
library(ggplot2)
library(tidyverse)
library(AICcmodavg)
library(lubridate)

#dev.off()

setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data")
#all_ssf1 <-read.csv("Woodie_SSF_w_perc.csv")
all_ssf1 <-read.csv("Cowan_2024_observed_quoll_locations_and_habitat_sampling_for_SSF.csv")
names(all_ssf1)

# Assuming all_ssf1 is your dataframe and t1_ is the column with timestamps
all_ssf1$t1_ <- as.POSIXct(all_ssf1$t1_, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")
all_ssf1$t2_ <- as.POSIXct(all_ssf1$t2_, format="%Y-%m-%dT%H:%M:%SZ", tz="UTC")

# Load the lubridate package
library(lubridate)
all_ssf1$realtime <- all_ssf1$t2_
# Extract the hour
hour_component <- hour(all_ssf1$t2_)
# Extract the minute and determine if .5 should be added
minute_component <- ifelse(minute(all_ssf1$t2_) == 30, 0.5, 0)
# Combine hour and minute components
all_ssf1$hour_from_t2_ <- hour_component + minute_component
# View the modified dataframe
head(all_ssf1)
# 
# # Define the original hours in order
# original_hours <- c(18.5, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
# # Create a sequence from 0 to 11.5 in steps of 0.5
# new_sequence <- seq(1, 24, by=1)
# # Create a mapping from original hours to new sequence
# hour_mapping <- setNames(new_sequence, original_hours)
# # Assuming hourly_coefs is your dataframe
# # Map the original hour values to the new sequence
# all_ssf1$newhour <- hour_mapping[as.character(all_ssf1$hour_from_t2_)]
constant=0.1
#all_ssf1$t1_<-dmy_hm(all_ssf1$t1_)
#all_ssf1$t2_<-dmy_hm(all_ssf1$t2_)
all_ssf2 <- all_ssf1 %>%
  mutate(id_num = as.numeric(factor(id)),
         step_id = step_id_,
         burst_ = burst_,
         x1 = x1_, x2 = x2_,
         y1 = y1_, y2 = y2_,
         t1 = t1_,
         t1_rounded = round_date(all_ssf1$t1_, "hour"),
         hour_t1 = hour(t1_rounded),
         t2 = t2_,
         t2_rounded = round_date(all_ssf1$t2_, "hour"),
         hour_t2 = hour(t2_rounded),
         hour = hour_t2,
         yday = yday(t1_),
         year = year(t1_),
         month = month(t1_),
         sl = sl_,
         log_sl = log(sl_),
         ta = ta_,
         cos_ta = cos(ta_),
         NDVI_class_end = NDVI_class_end,
         disturb_distance_end = disturb_distance_end,
         log_disturb = log(disturb_distance_end+constant),
         log_habitat = log(habitat_distance_end+constant),
         habitat_distance_end = habitat_distance_end,
         name=name,orig_name=orig_name,season=season,age=age,
         spatial_memory_ref = kde_ref_spatial_memory_density_log,
         spatial_memory_pi = kde_pi_spatial_memory_density_log,
         memory_log=kde_ref_spatial_memory_density_log,
         memory=exp(kde_ref_spatial_memory_density_log),
         hour_s1 = sin(2*pi*hour/24),
         hour_s2 = sin(4*pi*hour/24),
         hour_s3 = sin(6*pi*hour/24),
         hour_s4 = sin(8*pi*hour/24),
         hour_c1 = cos(2*pi*hour/24),
         hour_c2 = cos(4*pi*hour/24),
         hour_c3 = cos(6*pi*hour/24),
         hour_c4 = cos(8*pi*hour/24))


##Look at habitat distance through time (days)
ggplot(data = all_ssf2 %>% filter(case_ == 1,season=="breeding"),
       aes(x = t1_, y = habitat_distance_end, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic()

##Look at habitat distance through time (days)
ggplot(data = all_ssf2 %>% filter(case_ == 1,season=="breeding"),
       aes(x = t1_, y = memory, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic()

##Look at habitat distance through time (hours)
all_ssf2 <- all_ssf2 %>%
  mutate(hour_transformed = ifelse(hour < 6, hour + 24, hour))

ggplot(data = all_ssf2 %>% filter(case_ == 1, season == "breeding"),
       aes(x = hour_transformed, y = habitat_distance_end, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic() +
  coord_cartesian(xlim = c(18, 30)) # Adjusted based on the transformed hours


ggplot(data = all_ssf2 %>% filter(case_ == 1, season == "breeding"),
       aes(x = hour_transformed, y = memory, colour = factor(id))) +
  geom_point(alpha = 0.25) +
  scale_colour_viridis_d() +
  geom_smooth() +
  theme_classic() +
  coord_cartesian(xlim = c(18, 30)) # Adjusted based on the transformed hours


##Separate seasons for separate modelling
################For breeding season
all_ssf_breeding <- all_ssf2 %>% 
  filter(season == "breeding") %>% 
  arrange(t1_)

all_ssf_breeding$NDVI_class_end <-factor(all_ssf_breeding$NDVI_class_end, levels = c("rocky", "grassland", "dense_veg","mine_pit_waste_dump","other_disturbed"))

unique(all_ssf_breeding$NDVI_class_end)


quoll_ids <- unique(all_ssf_breeding$id)


quoll_data_matrix_unscaled <- all_ssf_breeding %>% transmute(
  
  ndvi = NDVI_class_end,
  
  habitat_distance = log_habitat,
  # habitat_distance_s1 = log_habitat * hour_s1,
  # habitat_distance_c1 = log_habitat * hour_c1,
  # habitat_distance_s2 = log_habitat * hour_s2, # New
  # habitat_distance_c2 = log_habitat * hour_c2,
  # habitat_distance_s3 = log_habitat * hour_s3, # New
  # habitat_distance_c3 = log_habitat * hour_c3,# New
  
  disturb_distance = log_disturb,
  # disturb_distance_s1 = log_disturb * hour_s1,
  # disturb_distance_c1 = log_disturb * hour_c1,
  # disturb_distance_s2 = log_disturb * hour_s2, # New
  # disturb_distance_c2 = log_disturb * hour_c2,
  # disturb_distance_s3 = log_disturb * hour_s3, # New
  # disturb_distance_c3 = log_disturb * hour_c3,# New
  
  memory = memory_log,
  # memory_s1 = memory_log * hour_s1,
  # memory_c1 = memory_log * hour_c1,
  # memory_s2 = memory_log * hour_s2, # New
  # memory_c2 = memory_log * hour_c2,
  # memory_s3 = memory_log * hour_s3, # New
  # memory_c3 = memory_log * hour_c3,# New
  
  #spatial_memory_pi = spatial_memory_pi ^ 2,
  #spatial_memory_pi_s1 = (spatial_memory_pi ^ 2) * hour_s1,
  #spatial_memory_pi_c1 = (spatial_memory_pi ^ 2) * hour_c1,
  #spatial_memory_pi_s2 = (spatial_memory_pi ^ 2) * hour_s2, # New
  #spatial_memory_pi_c2 = (spatial_memory_pi ^ 2) * hour_c2, # New
  
  step_l = sl,
  # step_l_s1 = sl * hour_s1,
  # step_l_c1 = sl * hour_c1,
  # step_l_s2 = sl * hour_s2, # New
  # step_l_c2 = sl * hour_c2, # New
  
  log_step_l = log_sl,
  # log_step_l_s1 = log_sl * hour_s1,
  # log_step_l_c1 = log_sl * hour_c1,
  # log_step_l_s2 = log_sl * hour_s2, # New
  # log_step_l_c2 = log_sl * hour_c2, # New
  # 
  cos_turn_a = cos_ta,
  # cos_turn_a_s1 = cos_ta * hour_s1,
  # cos_turn_a_c1 = cos_ta * hour_c1,
  # cos_turn_a_s2 = cos_ta * hour_s2, # New
  # cos_turn_a_c2 = cos_ta * hour_c2 # New
)


head(quoll_data_matrix_unscaled)
quoll_data_matrix_unscaled<-as.data.frame(quoll_data_matrix_unscaled)
## Identify numeric columns
#numeric_columns <- sapply(quoll_data_matrix_unscaled, is.numeric)
## Scale only numeric columns
#scaled_data <- scale(quoll_data_matrix_unscaled[, numeric_columns])
## Combine scaled numeric data with non-numeric columns
#quoll_data_matrix_scaled <- cbind(quoll_data_matrix_unscaled[!numeric_columns], scaled_data)
## Convert row names to the first column if necessary
#quoll_data_matrix_scaled <- data.frame(quoll_data_matrix_scaled, row.names = NULL)
## View the head of the scaled dataset
#head(quoll_data_matrix_scaled)

# Using base R
quoll_data_matrix_unscaled_no_ndvi <- quoll_data_matrix_unscaled[ , !(names(quoll_data_matrix_unscaled) %in% c("ndvi","step_l","log_step_l","cos_turn_a"))]
quoll_data_matrix_scaled_no_ndvi <- scale(quoll_data_matrix_unscaled_no_ndvi)

quoll_data_matrix_scaled <- cbind(ndvi = quoll_data_matrix_unscaled$ndvi, quoll_data_matrix_scaled_no_ndvi)


mean_vals <- attr(quoll_data_matrix_scaled_no_ndvi, "scaled:center")
sd_vals <- attr(quoll_data_matrix_scaled_no_ndvi, "scaled:scale")
scaling_attributes <- data.frame(variable = names(quoll_data_matrix_unscaled_no_ndvi), mean = mean_vals, sd = sd_vals)

quoll_data_scaled_0p <- data.frame(id = all_ssf_breeding$id,  step_id = all_ssf_breeding$step_id, case_ = all_ssf_breeding$case_, quoll_data_matrix_scaled)
head(quoll_data_scaled_0p)

# Ensure the row order and number of rows in both datasets match
# Replace the ndvi column in quoll_data_scaled_1p with that from quoll_data_matrix_unscaled
quoll_data_scaled_0p$ndvi <- quoll_data_matrix_unscaled$ndvi
quoll_data_scaled_0p$step_l <- quoll_data_matrix_unscaled$step_l
quoll_data_scaled_0p$log_step_l <- quoll_data_matrix_unscaled$log_step_l
quoll_data_scaled_0p$cos_turn_a <- quoll_data_matrix_unscaled$cos_turn_a


# View the first few rows to confirm the change
head(quoll_data_scaled_0p)



formula_twostep <- case_ ~ 
  
  ndvi +
  
  habitat_distance +
  # habitat_distance_s1 +
  # habitat_distance_c1 +
  # habitat_distance_s2 + # Added
  # habitat_distance_c2 +
  # habitat_distance_s3 + # Added
  # habitat_distance_c3 +# Added
  
  disturb_distance +
  # disturb_distance_s1 +
  # disturb_distance_c1 +
  # disturb_distance_s2 + # Added
  # disturb_distance_c2 +
  # disturb_distance_s3 + # Added
  # disturb_distance_c3 +# Added
  
  memory +
  # memory_s1 +
  # memory_c1 +
  # memory_s2 + # Added
  # memory_c2 +
  # memory_s3 + # Added
  # memory_c3 +# Added
  
  #spatial_memory_pi +
  #spatial_memory_pi_s1 +
  #spatial_memory_pi_c1 +
  #spatial_memory_pi_s2 + # Added
  #spatial_memory_pi_c2 + # Added
  
  step_l +
  # step_l_s1 +
  # step_l_c1 +
  # step_l_s2 + # Added
  # step_l_c2 + # Added
  
  log_step_l +
  # log_step_l_s1 +
  # log_step_l_c1 +
  # log_step_l_s2 + # Added
  # log_step_l_c2 + # Added
  
  cos_turn_a +
  # cos_turn_a_s1 +
  # cos_turn_a_c1 +
  # cos_turn_a_s2 + # Added
  # cos_turn_a_c2 + # Added
  
  strata(step_id) +
  cluster(id)


library(tictoc)
library(TwoStepCLogit)

  tic()
  model_twostep_0p_harms <- Ts.estim(formula = formula_twostep,
                                     data = quoll_data_scaled_0p,
                                     all.m.1 = TRUE,
                                     D = "UN(1)",
                                     itermax = 10000)
  toc()

#summary(model_twostep_0p_harms)
model_twostep_0p_harms
model_twostep_0p_harms$beta
model_twostep_0p_harms$se
model_twostep_0p_harms$vcov
diag(model_twostep_0p_harms$D) # between cluster variance
model_twostep_0p_harms$r.effect # individual estimates

hist(model_twostep_0p_harms$r.effect[,6])


coefs_clr <- data.frame(coefs = names(model_twostep_0p_harms$beta), value = model_twostep_0p_harms$beta)
coefs_clr_no_ndvi <- coefs_clr[!coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), ]

coefs_clr_no_ndvi$scale_sd <- scaling_attributes$sd
coefs_clr_no_ndvi <- coefs_clr_no_ndvi %>% mutate(value_nat = value / scale_sd)

head(coefs_clr)
head(coefs_clr_no_ndvi)

# Rows to add back
rows_to_add <- coefs_clr[coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), ]
rows_to_add$scale_sd <- NA
rows_to_add$value_nat <- coefs_clr[coefs_clr$coefs %in% c("ndvigrassland", "ndvidense_veg", "ndvimine_pit_waste_dump", "ndviother_disturbed","step_l","log_step_l","cos_turn_a"), "value"]
coefs_clr_upd <- rbind(coefs_clr_no_ndvi, rows_to_add)
print(coefs_clr_upd)

####Reconstructing coefficients with two pairs of harmonics, with quadratic terms

# we want to reconstruct the temporally dynamic coefficients for the full cycle, and then we can subset to the active hours after
hour <- seq(0.0,23.5,0.5)

all_ssf_breeding$hour30 <- format(as.POSIXct(all_ssf_breeding$t2), format = "%H:%M")
all_ssf_breeding$hour30 <- gsub(":", ".", all_ssf_breeding$hour30)
all_ssf_breeding$hour30 <- sub("30", "50", all_ssf_breeding$hour30)
all_ssf_breeding$hour30 <- as.numeric(all_ssf_breeding$hour30)

unique_hours<-unique(all_ssf_breeding$hour30)


# this creates a matrix of the hour harmonics, which is 24 x 3 (hours x harmonics, including the linear term)
hour_harmonics_df <- data.frame("linear_term" = rep(1, length(hour)))

hour_harmonics_linear_df <- data.frame("linear_term" = rep(1, length(hour)))


# we'll do the filtering after reconstructing the harmonics - i think it's best to do the whole cycle initially
# hour_harmonics_df <- hour_harmonics_df[hour_harmonics_df$hour %in% unique_hours, ]

hour_transformed = ifelse(hour < 6, hour + 24, hour)

# i spaced the coefficients out so i could see them a bit easier
harmonics_scaled_df_0p <- data.frame(
  "hour" = hour,
  # can add this in here now
  "hour_transformed" = ifelse(hour < 6, hour + 24, hour),
  
  # note that for the linear only terms i've added in the LINEAR ONLY data frame here, as there is only the linear term, and the matrix multiplication will work (now it's just an 1 x 1 %*% 1 x n, which will result in a 1 x n matrix (i.e. a vector))
  "ndvigrassland" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvigrassland", coefs)) %>%
                                 ###### under here ######
                               pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvidense_veg" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvidense_veg", coefs)) %>%
                                 pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvimine_pit_waste_dump" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvimine_pit_waste_dump", coefs)) %>%
                                           pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndviother_disturbed" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndviother_disturbed", coefs)) %>%
                                       pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  # For these ones, we will pull out the three coefficients for each variable (linear, sin1,cos1), which R considers to be COLUMN, and therefore a 3 x 1 matrix (or more simply a (column) vector). We transpose (the t() function) it to make it a 1 x 3. We then have the 24 x 3 matrix of the harmonics, which we TRANSPOSE to turn it into a 3 x 24 matrix. Then when we multiply them together, we have 1 x 3 %*% 3 x 24, which results in a 1 x 24 matrix (i.e. a VECTOR of the temporally dynamic coefficient throughout the day), which we add into the dataframe
  
  "step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("step_l", coefs) & !grepl("log", coefs)) %>% 
                          pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  # add in the log_sl term
  "log_step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("log_step_l", coefs)) %>% 
                              pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "cos_turn_a" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("cos_turn_a", coefs)) %>% 
                              pull(value) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "disturb_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("disturb_distance", coefs) & !grepl("sq", coefs)) %>%
                                    pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "habitat_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("habitat_distance", coefs)) %>% 
                                    pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "memory" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("memory", coefs)) %>% 
                          pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
  #,
  
  #"spatial_memory_pi" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("spatial_memory_pi", coefs)) %>% 
  #                                   pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
)



harmonics_nat_df_0p <- data.frame(
  "hour" = hour,
  # can add this in here now
  "hour_transformed" = ifelse(hour < 6, hour + 24, hour),
  
  # note that for the linear only terms i've added in the LINEAR ONLY data frame here, as there is only the linear term, and the matrix multiplication will work (now it's just an 1 x 1 %*% 1 x n, which will result in a 1 x n matrix (i.e. a vector))
  "ndvigrassland" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvigrassland", coefs)) %>%
                                 ###### under here ######
                               pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvidense_veg" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvidense_veg", coefs)) %>%
                                 pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndvimine_pit_waste_dump" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndvimine_pit_waste_dump", coefs)) %>%
                                           pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "ndviother_disturbed" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("ndviother_disturbed", coefs)) %>%
                                       pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  # For these ones, we will pull out the three coefficients for each variable (linear, sin1,cos1), which R considers to be COLUMN, and therefore a 3 x 1 matrix (or more simply a (column) vector). We transpose (the t() function) it to make it a 1 x 3. We then have the 24 x 3 matrix of the harmonics, which we TRANSPOSE to turn it into a 3 x 24 matrix. Then when we multiply them together, we have 1 x 3 %*% 3 x 24, which results in a 1 x 24 matrix (i.e. a VECTOR of the temporally dynamic coefficient throughout the day), which we add into the dataframe
  
  "step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("step_l", coefs) & !grepl("log", coefs)) %>% 
                          pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  # add in the log_sl term
  "log_step_l" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("log_step_l", coefs)) %>%
                              pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "cos_turn_a" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("cos_turn_a", coefs)) %>%
                              pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_linear_df))),
  
  "disturb_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("disturb_distance", coefs) & !grepl("sq", coefs)) %>%
                                    pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "habitat_distance" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("habitat_distance", coefs)) %>% 
                                    pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df))),
  
  "memory" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("memory", coefs)) %>% 
                          pull(value_nat) %>% t() %*% t(as.matrix(hour_harmonics_df)))
  #,
  
  #"spatial_memory_pi" = as.numeric(coefs_clr_upd %>% dplyr::filter(grepl("spatial_memory_pi", coefs)) %>% 
  #                                   pull(value) %>% t() %*% t(as.matrix(hour_harmonics_df)))
)




##Get the tentative distributions
all_ssf_breed_true <- all_ssf_breeding %>% 
  filter(case_ == "1") %>% 
  arrange(t1_)
# Fit a gamma distribution to the step length of true steps
fit_sl <- fit_distr(all_ssf_breed_true$sl_, "gamma")
# Extract the shape and rate parameters of the gamma distribution
shape_sl <- fit_sl$params[1]
scale_sl <- fit_sl$params[2]
# Load the circular package
library(circular)
# Fit a von Mises distribution to the ta_ values
fit_ta <- fit_distr(all_ssf_breed_true$ta_, "vonmises")
# Extract the kappa parameter of the von Mises distribution
kappa_ta <- fit_ta$params["kappa"]


tentative_shape <- shape_sl$shape
tentative_scale <- scale_sl$scale
tentative_kappa <- kappa_ta$kappa

hist(rgamma(1e4,
            shape = tentative_shape,
            scale = tentative_scale),
     breaks = 100)

harmonics_nat_scaled_df_0p <- harmonics_nat_df_0p %>% mutate(shape = tentative_shape + log_step_l,
                                                             scale = 1/((1/tentative_scale) - step_l),
                                                             kappa = tentative_kappa + cos_turn_a)

hist(rgamma(1e4,
            shape = harmonics_nat_scaled_df_0p$shape[1],
            scale = harmonics_nat_scaled_df_0p$scale[1]),
     breaks = 100)




harmonics_nat_scaled_df_0p
# Add a new column 'rocky' with 0s to the dataframe 'harmonics_nat_scaled_df_0p'
harmonics_nat_scaled_df_0p$ndvirocky <- 0

harmonics_scaled_long_0p <- pivot_longer(harmonics_nat_scaled_df_0p, cols = !1:2, names_to = "coef")

ggplot() +
  geom_path(data = harmonics_scaled_long_0p,
            aes(x = hour, y = value, colour = coef)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour") +
  scale_color_discrete("Estimate") +
  theme_classic() +
  theme(legend.position = "bottom")


# when removing the daytime hours
harmonics_scaled_df_0p_active <- harmonics_nat_scaled_df_0p %>% dplyr::filter(hour %in% unique_hours)
harmonics_scaled_long_0p_active <- pivot_longer(harmonics_scaled_df_0p_active, cols = !1:2, names_to = "coef")
##write.csv(harmonics_scaled_df_0p_active,"C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_0hrm.csv")
movement_summary <- all_ssf_breeding %>% filter(case_ == 1) %>%  group_by(hour30) %>% summarise(mean_sl = mean(sl), median_sl = median(sl))
movement_summary$hour <- movement_summary$hour30

final0<-cbind(harmonics_scaled_df_0p_active,movement_summary)
final0<-final0[1:20]
final0<-final0 %>% dplyr::select(-hour30)

#write.csv(final0,"C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_0hrm.csv",row.names = FALSE)


harmonics_scaled_df_0p_active<-read.csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/iSSF coefficients for each harmonic pair/hours_temporal_coefficients_0hrm.csv")


library(ggplot2)
# Generate a set of random colors
# First, find the number of unique levels in the 'coef' column
num_levels <- length(unique(harmonics_scaled_long_0p_active$coef))

# Then, generate random colors for each level
# You can adjust the seed for reproducibility
set.seed(123) # Remove or change the seed value for different colors
random_colors <- grDevices::rainbow(num_levels)

# Now, create the plot
ggplot(harmonics_scaled_long_0p_active, aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(size = 2.5) +  # Increased line thickness
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5","6"),
                     limits = c(19, 30)) +
  scale_color_manual("Estimate", values = setNames(random_colors, unique(harmonics_scaled_long_0p_active$coef))) +
  theme_classic() +
  theme(legend.position = "bottom")



ggplot() +
  geom_path(data = harmonics_scaled_long_0p %>%
              filter(!coef %in% c("step_l", "log_step_l", "cos_turn_a", "memory")),
            aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5","6"),
                     limits = c(19, 30)) +
  scale_color_manual("Estimate", values = setNames(random_colors, unique(harmonics_scaled_long_0p_active$coef))) +
  theme_classic() +
  theme(legend.position = "bottom")




ggplot() +
  geom_path(data = harmonics_scaled_long_0p %>%
              filter(coef %in% c("step_l", "log_step_l", "cos_turn_a")),
            aes(x = hour_transformed, y = value, colour = coef)) +
  geom_line(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(expression(Time-varying~parameter~values~beta)) +
  scale_x_continuous("Hour", 
                     breaks = c(19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
                     labels = c("19", "20", "21", "22", "23", "0", "1", "2", "3", "4", "5","6"),
                     limits = c(19, 30)) +
  scale_color_manual("Estimate", values = setNames(random_colors, unique(harmonics_scaled_long_0p_active$coef))) +
  theme_classic() +
  theme(legend.position = "bottom")



