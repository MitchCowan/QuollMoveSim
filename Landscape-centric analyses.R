################################################################################
# This code is adapted from the methods described in:
#
#   Hofmann, D. D., Cozzi, G., McNutt, J. W., Ozgul, A., & Behr, D. M. (2023).
#   A three-step approach for assessing landscape connectivity via simulated 
#   dispersal: African wild dog case study.
#   *Landscape Ecology*, 38, 981–998.
#   https://doi.org/10.1007/s10980-023-01602-4
#
# Core components of this workflow—trajectory interpolation, cell transition 
# mapping, and betweenness centrality—were adapted to suit a northern quoll 
# case study across contrasting mining landscape scenarios.
################################################################################




########################
#### Rocky Habitat Connectivity
########################

# 1. Load required packages
library(tidyverse)
library(adehabitatLT)
library(raster)
library(sp)

# 2. Clear environment
rm(list = ls())

# 3. Read all simulation trajectory data
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Simulation trajectories")
data_all <- read.csv("all_simulation_trajectories_without_memory.csv", stringsAsFactors = FALSE)

# 4. Format timestamps and ID variables
data_all$t_ <- ifelse(
  nchar(as.character(data_all$t_)) <= 10,
  paste0(as.character(data_all$t_), " 00:00:00"),
  as.character(data_all$t_)
)
data_all$t_ <- as.POSIXct(data_all$t_, format = "%Y-%m-%d %H:%M")
data_all$id <- as.factor(data_all$id)


# =====================
# 5. Process NON-MINING landscape
# =====================

# A. Filter and order tracks
data_nm <- data_all %>% 
  filter(landscape == "non_mining") %>% 
  arrange(id, t_)

# B. Nest individual tracks
animals_nm <- data_nm %>% 
  group_by(id) %>% 
  summarise(track = list(tibble(x_ = x_, y_ = y_, t_ = t_))) %>% 
  ungroup()

# C. Load habitat raster
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/non_mining")
h_cat_nm <- raster(paste("Non_mining_habitat_classifications", "tif", sep = "."))

# D. Identify rocky patches (class 3)
rocky_binary_nm  <- reclassify(h_cat_nm, cbind(c(1,2,3,4,5), c(NA,NA,1,NA,NA)))
rocky_patches_nm <- clump(rocky_binary_nm, directions = 8)

# E. Extract patch IDs for steps
tracks_nm <- animals_nm %>% unnest(track)
tracks_nm$patch_id <- raster::extract(rocky_patches_nm, cbind(tracks_nm$x_, tracks_nm$y_))
tracks_nm <- tracks_nm %>% arrange(t_) %>% mutate(row_number = row_number())

# F. Calculate transitions between rocky patches
tracks_next_nm <- tracks_nm %>%
  filter(!is.na(patch_id)) %>%
  mutate(
    next_patch_id = lead(patch_id),
    next_row_num  = lead(row_number)
  ) %>%
  na.omit()

tracks_diff_nm <- tracks_next_nm %>%
  group_by(patch_id, next_patch_id) %>%
  summarise(
    avg_steps = mean(next_row_num - row_number, na.rm = TRUE),
    frequency = n(),
    .groups = "drop"
  ) %>%
  mutate(landscape = "non_mining")


# =====================
# 6. Process ORIGINAL (current mining) landscape
# =====================

data_orig <- data_all %>% 
  filter(landscape == "original") %>% 
  arrange(id, t_)
animals_orig <- data_orig %>% 
  group_by(id) %>% 
  summarise(track = list(tibble(x_ = x_, y_ = y_, t_ = t_))) %>% 
  ungroup()
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/current mining")
h_cat_orig <- raster(paste("Current_habitat_classifications", "tif", sep = "."))
rocky_binary_orig  <- reclassify(h_cat_orig, cbind(c(1,2,3,4,5), c(NA,NA,1,NA,NA)))
rocky_patches_orig <- clump(rocky_binary_orig, directions = 8)
tracks_orig <- animals_orig %>% unnest(track)
tracks_orig$patch_id <- raster::extract(rocky_patches_orig, cbind(tracks_orig$x_, tracks_orig$y_))
tracks_orig <- tracks_orig %>% arrange(t_) %>% mutate(row_number = row_number())
tracks_next_orig <- tracks_orig %>%
  filter(!is.na(patch_id)) %>%
  mutate(
    next_patch_id = lead(patch_id),
    next_row_num  = lead(row_number)
  ) %>%
  na.omit()
tracks_diff_orig <- tracks_next_orig %>%
  group_by(patch_id, next_patch_id) %>%
  summarise(
    avg_steps = mean(next_row_num - row_number, na.rm = TRUE),
    frequency = n(),
    .groups = "drop"
  ) %>%
  mutate(landscape = "original")


# =====================
# 7. Process FRAGMENTED (dispersed mining) landscape
# =====================

data_frag <- data_all %>% 
  filter(landscape == "fragmented") %>% 
  arrange(id, t_)
animals_frag <- data_frag %>% 
  group_by(id) %>% 
  summarise(track = list(tibble(x_ = x_, y_ = y_, t_ = t_))) %>% 
  ungroup()
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/dispersed mining")
h_cat_frag <- raster(paste("Dispersed_habitat_classifications", "tif", sep = "."))
rocky_binary_frag  <- reclassify(h_cat_frag, cbind(c(1,2,3,4,5), c(NA,NA,1,NA,NA)))
rocky_patches_frag <- clump(rocky_binary_frag, directions = 8)
tracks_frag <- animals_frag %>% unnest(track)
tracks_frag$patch_id <- raster::extract(rocky_patches_frag, cbind(tracks_frag$x_, tracks_frag$y_))
tracks_frag <- tracks_frag %>% arrange(t_) %>% mutate(row_number = row_number())
tracks_next_frag <- tracks_frag %>%
  filter(!is.na(patch_id)) %>%
  mutate(
    next_patch_id = lead(patch_id),
    next_row_num  = lead(row_number)
  ) %>%
  na.omit()
tracks_diff_frag <- tracks_next_frag %>%
  group_by(patch_id, next_patch_id) %>%
  summarise(
    avg_steps = mean(next_row_num - row_number, na.rm = TRUE),
    frequency = n(),
    .groups = "drop"
  ) %>%
  mutate(landscape = "dispersed")


# =====================
# 8. Process DENSE (aggregated mining) landscape
# =====================

data_dense <- data_all %>% 
  filter(landscape == "dense") %>% 
  arrange(id, t_)
animals_dense <- data_dense %>% 
  group_by(id) %>% 
  summarise(track = list(tibble(x_ = x_, y_ = y_, t_ = t_))) %>% 
  ungroup()
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/aggregated mining")
h_cat_dense <- raster(paste("Aggregated habitat classification", "tif", sep = "."))
rocky_binary_dense  <- reclassify(h_cat_dense, cbind(c(1,2,3,4,5), c(NA,NA,1,NA,NA)))
rocky_patches_dense <- clump(rocky_binary_dense, directions = 8)
tracks_dense <- animals_dense %>% unnest(track)
tracks_dense$patch_id <- raster::extract(rocky_patches_dense, cbind(tracks_dense$x_, tracks_dense$y_))
tracks_dense <- tracks_dense %>% arrange(t_) %>% mutate(row_number = row_number())
tracks_next_dense <- tracks_dense %>%
  filter(!is.na(patch_id)) %>%
  mutate(
    next_patch_id = lead(patch_id),
    next_row_num  = lead(row_number)
  ) %>%
  na.omit()
tracks_diff_dense <- tracks_next_dense %>%
  group_by(patch_id, next_patch_id) %>%
  summarise(
    avg_steps = mean(next_row_num - row_number, na.rm = TRUE),
    frequency = n(),
    .groups = "drop"
  ) %>%
  mutate(landscape = "aggregated")


# =====================
# 9. Combine all landscapes
# =====================

rocky_habitat_connectivity <- bind_rows(
  tracks_diff_nm,
  tracks_diff_orig,
  tracks_diff_frag,
  tracks_diff_dense
)

rocky_habitat_connectivity$X <-row_number(rocky_habitat_connectivity)

# Optional: export to CSV
# write.csv(
#   rocky_habitat_connectivity,
#   "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Data used for models/rocky habitat connectivity.csv",
#   row.names = FALSE
# )
# 
# message("✅ rocky_habitat_connectivity.csv written with all landscapes.")


#############################
### Connectivity regression #
#############################

# 10. Load data for regression
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Data used for models")
data <- read.csv("rocky habitat connectivity.csv", stringsAsFactors = FALSE)

# 11. Clean & format variables
library(dplyr)
data <- data %>%
  mutate(
    landscape = relevel(factor(landscape), ref = "non_mining"),
    avg_steps = as.numeric(avg_steps)
  )

# 12. Explore distribution of response variable
library(ggplot2)
ggplot(data, aes(x = avg_steps)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  labs(title = "Distribution of avg_steps", x = "avg_steps", y = "Density")

hist(data$avg_steps)
plot(density(data$avg_steps))
lines(density(rnorm(1000, mean = mean(data$avg_steps), sd = sd(data$avg_steps))), col = "blue")

# 13. Fit GLMM (Gamma, log-link)
library(glmmTMB)
avg_steps_model <- glmmTMB(avg_steps ~ landscape, data = data, family = Gamma(link = "log"))
summary(avg_steps_model)

# 14. Extract and format model coefficients
extract_coefficients <- function(model, model_name) {
  s <- summary(model)$coeff$cond
  df <- as.data.frame(s)
  df$term <- rownames(df)
  df <- df %>%
    mutate(
      conf.low  = Estimate - 1.96 * `Std. Error`,
      conf.high = Estimate + 1.96 * `Std. Error`,
      color     = ifelse(`Pr(>|z|)` < 0.05, "#9CC065", "#FF9F5F"),
      model     = model_name
    )
  rownames(df) <- NULL
  df
}

# 15. Visualise coefficients
plot_coefficients <- function(coefs, title = "") {
  coefs <- coefs %>% filter(term != "(Intercept)")
  ggplot(coefs, aes(x = Estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, size = 1) +
    geom_point(aes(color = color), size = 3) +
    scale_color_identity() +
    labs(title = title, x = "Revisitations coefficient (β)", y = "Landscape") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10),
      plot.title = element_text(hjust = 0.5)
    )
}

# 16. Plot all landscapes
avg_steps_coeff <- extract_coefficients(avg_steps_model, "All landscapes") %>%
  mutate(term = recode(term,
                       "landscapeoriginal"   = "Current mining",
                       "landscapeaggregated" = "Aggregated mining",
                       "landscapedispersed"  = "Dispersed mining"
  )) %>%
  mutate(term = factor(term, levels = c("Aggregated mining", "Current mining", "Dispersed mining")))

avg_steps_plot <- plot_coefficients(avg_steps_coeff, "Mining Composition Analysis")
print(avg_steps_plot)

# 17. Secondary model: mining-only subset
library(forcats)
data_mining <- data %>%
  filter(landscape != "non_mining") %>%
  mutate(landscape = fct_relevel(landscape, "original"))

avg_steps_mining_model <- glmmTMB(avg_steps ~ landscape, data = data_mining, family = Gamma(link = "log"))
summary(avg_steps_mining_model)

# 18. Plot mining-only model
mining_coeff <- extract_coefficients(avg_steps_mining_model, "Mining-only") %>%
  mutate(term = recode(term,
                       "(Intercept)"           = "Current mining",
                       "landscapeaggregated"   = "Aggregated mining",
                       "landscapedispersed"    = "Dispersed mining"
  )) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = factor(term, levels = c("Aggregated mining", "Current mining", "Dispersed mining")))

mining_plot <- plot_coefficients(mining_coeff, "Mining-only Composition")
print(mining_plot)





################################################################################
#### Heatmap Generator for Simulated Quoll Trajectories
################################################################################

# -----------------------------
# 1. Setup Environment
# -----------------------------
rm(list = ls())

library(tidyverse)
library(raster)
library(viridis)
library(rgeos)
library(sf)
library(spatstat)
library(maptools)
library(igraph)
library(ggpubr)
library(ggnetwork)
library(adehabitatLT)
library(lubridate)
library(tibble)
library(viridis)

# -----------------------------
# 2. Load Simulation Data (All Landscapes)
# -----------------------------
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Simulation trajectories")
data_all <- read.csv("all_simulation_trajectories_without_memory.csv", stringsAsFactors = FALSE)

# Format datetime and ID
data_all$t_ <- ifelse(nchar(as.character(data_all$t_)) <= 10,
                      paste0(data_all$t_, " 00:00:00"),
                      as.character(data_all$t_))
data_all$t_ <- as.POSIXct(data_all$t_, format = "%Y-%m-%d %H:%M")
data_all$id <- as.factor(data_all$id)

# -----------------------------
# 3. Define Utility Function
# -----------------------------
rasterizeSpatstat <- function(l, r){
  values(r) <- 0
  im <- as.im.RasterLayer(r)
  summed <- im
  for (y in 1:length(l)){
    line    <- as.psp(l[y, ], window = im)
    line    <- as.mask.psp(line)
    line_r  <- as.im.owin(line, na.replace = 0)
    summed  <- Reduce("+", list(summed, line_r))
  }
  return(raster(summed))
}

# -----------------------------
# 4. Define Landscapes and Updated Raster Paths
# -----------------------------
landscapes <- list(
  original = "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/current mining/Current_habitat_classifications.tif",
  fragmented = "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/dispersed mining/Dispersed_habitat_classifications.tif",
  dense = "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/aggregated mining/Aggregated habitat classification.tif",
  non_mining = "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/non_mining/Non_mining_habitat_classifications.tif"
)

output_dir <- "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps"

# -----------------------------
# 5. Loop Over Landscapes
# -----------------------------
for (land in names(landscapes)) {
  
  cat("\nProcessing:", land, "landscape\n")
  
  # Subset data
  data <- data_all %>% filter(landscape == land) %>% arrange(id, t_)
  
  # Convert each individual’s trajectory to a spatial line
  tracks <- lapply(unique(data$id), function(x) {
    sub <- data[data$id == x, ]
    coordinates(sub) <- c("x_", "y_")
    spLines(sub)
  })
  tracks <- do.call(rbind, tracks)
  tracks$ID <- 1:length(tracks)
  
  # Load correct habitat raster and create template raster
  h_cat <- raster(landscapes[[land]])
  heatmap <- raster(h_cat)
  
  # Rasterize the movement tracks
  heatmap <- rasterizeSpatstat(tracks, heatmap)
  heatmap <- crop(heatmap, h_cat)
  
  # Plot the heatmap using magma palette
  plot_heatmap <- ggplot(as.data.frame(heatmap, xy = TRUE)) +
    geom_raster(aes(x = x, y = y, fill = layer)) +
    scale_fill_gradientn(
      colors = viridis::magma(100),
      name = "Traversal Frequency",
      guide = guide_colorbar(
        show.limits = TRUE,
        title.position = "top",
        title.hjust = 0.5,
        ticks = TRUE,
        barheight = unit(0.3, "cm"),
        barwidth = unit(5.0, "cm")
      )
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    ggtitle(paste("Heatmap -", land)) +
    coord_fixed()  # Use coord_fixed() instead of coord_sf() for raster data
  
  print(plot_heatmap)
  
  # Write raster to disk
  # writeRaster(
  #   heatmap,
  #   filename = file.path(output_dir, paste0("heatmap_", land, ".tif")),
  #   format = "GTiff",
  #   overwrite = TRUE
  # )
}





##############################################
#### Betweenness Maps
##############################################

# Load packages
library(tidyverse)
library(raster)
library(sp)
library(igraph)
library(adehabitatLT)
library(lubridate)
library(Matrix)
library(viridis)
library(purrr)

# Clear environment
rm(list = ls())

# Read all trajectory data
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Simulation trajectories")
data_all <- read.csv("all_simulation_trajectories_without_memory.csv", stringsAsFactors = FALSE)

# Format timestamps and IDs
data_all$t_ <- ifelse(nchar(as.character(data_all$t_)) <= 10,
                      paste0(data_all$t_, " 00:00:00"),
                      as.character(data_all$t_))
data_all$t_ <- as.POSIXct(data_all$t_, format = "%Y-%m-%d %H:%M")
data_all$id <- as.factor(data_all$id)

# Function to interpolate steps every ~10m
interpolate_steps_metric_with_time <- function(track) {
  interpolated_points <- list()
  track_df <- as.data.frame(track)
  
  for (i in 1:(nrow(track_df) - 1)) {
    start_point <- as.numeric(track_df[i, c("x_", "y_")])
    end_point <- as.numeric(track_df[i + 1, c("x_", "y_")])
    start_time <- as.POSIXct(track_df$t_[i])
    
    dist <- sqrt(sum((end_point - start_point)^2))
    npts <- max(1, floor(dist / 10)) - 1
    
    if(npts > 0) {
      x_seq <- seq(from = start_point[1], to = end_point[1], length.out = npts + 2)[-c(1, npts + 2)]
      y_seq <- seq(from = start_point[2], to = end_point[2], length.out = npts + 2)[-c(1, npts + 2)]
      time_seq <- seq(from = start_time, by = "min", length.out = npts + 2)[-c(1, npts + 2)]
      interpolated_points[[i]] <- data.frame(x_ = x_seq, y_ = y_seq, t_ = time_seq)
    }
  }
  
  do.call(rbind, interpolated_points)
}

# Function to extract cell transitions from interpolated track
extract_transitions_sparse <- function(track, grid_raster) {
  coords <- cbind(track$x_, track$y_)
  cells <- cellFromXY(grid_raster, coords)
  cells <- cells[!is.na(cells)]
  if (length(cells) < 2) {
    return(sparseMatrix(dims = c(ncell(grid_raster), ncell(grid_raster))))
  }
  sparseMatrix(i = head(cells, -1), j = tail(cells, -1), x = 1, dims = c(ncell(grid_raster), ncell(grid_raster)))
}

############################
##### Non-mining ###########
############################

# Filter and prepare
data <- data_all %>% filter(landscape == "non_mining") %>% arrange(id, t_)
animals <- data %>%
  group_by(id) %>% summarise(track = list(tibble(x_ = x_, y_ = y_, t_ = t_)), .groups = "drop") %>%
  mutate(track = map(track, interpolate_steps_metric_with_time))

# Load habitat raster and create grid
h_cat <- raster("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/non_mining/Non_mining_habitat_classifications.tif")
grid_raster <- aggregate(h_cat, fact = 257 / 10, fun = max)

# Build transition matrix
matrices <- lapply(animals$track, extract_transitions_sparse, grid_raster)
transition_matrix <- Reduce(`+`, matrices)

# Build edge list and graph
nonzero <- which(transition_matrix != 0, arr.ind = TRUE)
edge_list <- data.frame(from = nonzero[, 1], to = nonzero[, 2], weight = transition_matrix[nonzero])
g <- graph_from_data_frame(edge_list, directed = FALSE)

# Compute betweenness
bt <- betweenness(g, weights = E(g)$weight)
grid_raster[] <- 0
grid_raster[as.numeric(names(bt))] <- bt

# Save and plot
plot(grid_raster, col = magma(100))
#writeRaster(grid_raster, "Interpolated_betweenness_non_mining.tif", format = "GTiff", overwrite = TRUE)

############################
##### Current mining #######
############################
data <- data_all %>% filter(landscape == "original") %>% arrange(id, t_)
animals <- data %>%
  group_by(id) %>% summarise(track = list(tibble(x_ = x_, y_ = y_, t_ = t_)), .groups = "drop") %>%
  mutate(track = map(track, interpolate_steps_metric_with_time))

h_cat <- raster("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/current mining/Current_habitat_classifications.tif")
grid_raster <- aggregate(h_cat, fact = 257 / 10, fun = max)

matrices <- lapply(animals$track, extract_transitions_sparse, grid_raster)
transition_matrix <- Reduce(`+`, matrices)

nonzero <- which(transition_matrix != 0, arr.ind = TRUE)
edge_list <- data.frame(from = nonzero[, 1], to = nonzero[, 2], weight = transition_matrix[nonzero])
g <- graph_from_data_frame(edge_list, directed = FALSE)

bt <- betweenness(g, weights = E(g)$weight)
grid_raster[] <- 0
grid_raster[as.numeric(names(bt))] <- bt

plot(grid_raster, col = magma(100))
#writeRaster(grid_raster, "Interpolated_betweenness_original.tif", format = "GTiff", overwrite = TRUE)

############################
##### Dispersed mining #####
############################

data <- data_all %>% filter(landscape == "fragmented") %>% arrange(id, t_)
animals <- data %>%
  group_by(id) %>% summarise(track = list(tibble(x_ = x_, y_ = y_, t_ = t_)), .groups = "drop") %>%
  mutate(track = map(track, interpolate_steps_metric_with_time))

h_cat <- raster("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/dispersed mining/Dispersed_habitat_classifications.tif")
grid_raster <- aggregate(h_cat, fact = 257 / 10, fun = max)

matrices <- lapply(animals$track, extract_transitions_sparse, grid_raster)
transition_matrix <- Reduce(`+`, matrices)

nonzero <- which(transition_matrix != 0, arr.ind = TRUE)
edge_list <- data.frame(from = nonzero[, 1], to = nonzero[, 2], weight = transition_matrix[nonzero])
g <- graph_from_data_frame(edge_list, directed = FALSE)

bt <- betweenness(g, weights = E(g)$weight)
grid_raster[] <- 0
grid_raster[as.numeric(names(bt))] <- bt

plot(grid_raster, col = magma(100))
#writeRaster(grid_raster, "Interpolated_betweenness_fragmented.tif", format = "GTiff", overwrite = TRUE)

############################
##### Aggregated mining ####
############################

data <- data_all %>% filter(landscape == "dense") %>% arrange(id, t_)
animals <- data %>%
  group_by(id) %>% summarise(track = list(tibble(x_ = x_, y_ = y_, t_ = t_)), .groups = "drop") %>%
  mutate(track = map(track, interpolate_steps_metric_with_time))

h_cat <- raster("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/aggregated mining/Aggregated habitat classification.tif")
grid_raster <- aggregate(h_cat, fact = 257 / 10, fun = max)

matrices <- lapply(animals$track, extract_transitions_sparse, grid_raster)
transition_matrix <- Reduce(`+`, matrices)

nonzero <- which(transition_matrix != 0, arr.ind = TRUE)
edge_list <- data.frame(from = nonzero[, 1], to = nonzero[, 2], weight = transition_matrix[nonzero])
g <- graph_from_data_frame(edge_list, directed = FALSE)

bt <- betweenness(g, weights = E(g)$weight)
grid_raster[] <- 0
grid_raster[as.numeric(names(bt))] <- bt

plot(grid_raster, col = magma(100))
#writeRaster(grid_raster, "Interpolated_betweenness_dense.tif", format = "GTiff", overwrite = TRUE)
