# ============================================================
# Home Range Analysis and Coefficient Plotting - Simulated Data
# ============================================================

# -----------------------------
# 0. Load Required Libraries
# -----------------------------
library(tidyverse)
packages <- c("amt", "lubridate", "terra", "tictoc", "beepr", "ks",
              "adehabitatLT", "adehabitatHR", "ggpubr", "patchwork", "dplyr")
walk(packages, require, character.only = TRUE)

rm(list = ls())  # Clear the workspace

# -----------------------------
# 1. Load Simulation Trajectory Data
# -----------------------------
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Simulation trajectories")

# Read in the combined simulation trajectories
data <- read.csv("all_simulation_trajectories_with_memory.csv")
head(data)
str(data)
unique(data$landscape)
unique(data$id)

# -----------------------------
# 2. Load Precomputed Home Range Areas
# -----------------------------
library(sf)
library(dplyr)

result_df <- read.csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Data used for models/home_ranges.csv")

# Factor relevel for reference category
result_df <- result_df %>%
  mutate(landscape = relevel(factor(landscape), ref = "non_mining"))

# Convert area to numeric home range and standardise ID
result_df$home_range <- as.numeric(result_df$Area_m2)
result_df$new_id <- sub(".$", "", result_df$id)

# -----------------------------
# 3. Exploratory Plots - Home Range Distribution
# -----------------------------
ggplot(result_df, aes(x = home_range)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  geom_density(alpha = 0.2, fill = "#FF6666") +
  labs(title = "Distribution of Home Range", x = "Home Range", y = "Density")

hist(result_df$home_range)
plot(density(result_df$home_range))
lines(density(rnorm(1000, mean = mean(result_df$home_range), sd = sd(result_df$home_range))), col = "blue")

# -----------------------------
# 4. GLMM for All Landscapes
# -----------------------------
library(glmmTMB)
hr_model <- glmmTMB(home_range ~ landscape + (1 | new_id), data = result_df, family = Gamma(link = "log"))
summary(hr_model)

# -----------------------------
# 5. Subset and Relevel Mining Models
# -----------------------------
library(forcats)

# Exclude non-mining and set "original" as reference
result_df_mining <- result_df %>%
  filter(landscape != "non_mining") %>%
  mutate(landscape = fct_relevel(landscape, "original"))

table(result_df_mining$landscape)

# Model with only mining landscapes
hr_model_mining <- glmmTMB(home_range ~ landscape + (1 | new_id), data = result_df_mining, family = Gamma(link = "log"))
summary(hr_model_mining)

# -----------------------------
# 6. Coefficient Extraction Function
# -----------------------------
extract_coefficients <- function(model, model_name) {
  model_summary <- summary(model)
  conditional_coefficients <- as.data.frame(model_summary$coeff$cond)
  conditional_coefficients$term <- rownames(conditional_coefficients)
  conditional_coefficients$conf.low <- conditional_coefficients$Estimate - 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$conf.high <- conditional_coefficients$Estimate + 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$color <- ifelse(conditional_coefficients$`Pr(>|z|)` < 0.05, "black", "grey")
  conditional_coefficients$model <- model_name
  return(conditional_coefficients)
}

# -----------------------------
# 7. Coefficient Plotting Function
# -----------------------------
plot_coefficients <- function(coefficients) {
  coefficients <- coefficients[coefficients$term != "(Intercept)", ]
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#FF9F5F", "gray"))
  
  ggplot(coefficients, aes(x = Estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(title = "Non-mining vs. mining analysis", x = "Home range coefficient (β)", y = "Landscape") +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5)
    ) +
    if ("model" %in% names(coefficients)) facet_grid(model ~ ., scales = "free_y", space = "free_y") else NULL
}

# -----------------------------
# 8. Prepare and Plot: All Landscapes
# -----------------------------
hr_coeff <- extract_coefficients(hr_model, "landscape")

hr_coeff <- hr_coeff %>%
  mutate(term = case_when(
    term == "landscapeoriginal" ~ "Current mining",
    term == "landscapedense" ~ "Consolidated mining",
    term == "landscapefragmented" ~ "Fragmented mining",
    term == "(Intercept)" ~ "Non-mining",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Non-mining", "Consolidated mining", "Current mining", "Fragmented mining"))) %>%
  filter(term != "Non-mining")

# Plot
hr_plot <- plot_coefficients(hr_coeff)
print(hr_plot)

# -----------------------------
# 9. Prepare and Plot: Mining Only
# -----------------------------
hr_mining_coeff <- extract_coefficients(hr_model_mining, "landscape")

hr_mining_coeff <- hr_mining_coeff %>%
  mutate(term = case_when(
    term == "landscapeoriginal" ~ "Current mining",
    term == "landscapedense" ~ "Consolidated mining",
    term == "landscapefragmented" ~ "Fragmented mining",
    term == "(Intercept)" ~ "Non-mining",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Non-mining", "Consolidated mining", "Current mining", "Fragmented mining"))) %>%
  filter(term != "Non-mining")

# Plot
hr_mining_plot <- plot_coefficients(hr_mining_coeff)
print(hr_mining_plot)




# =========================================================
# Movement Cost Analysis
# =========================================================

# -----------------------------
# Load Libraries and Set Environment
# -----------------------------
library(tidyverse)
library(terra)
library(glmmTMB)
library(forcats)
library(ggplot2)
library(plyr)

# Set working directory
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Data used for models")

# -----------------------------
# Load and Prepare Data
# -----------------------------
data <- read.csv("movement_costs.csv")

data <- data %>%
  mutate(landscape = relevel(factor(landscape), ref = "non_mining"))
data$cost_mean <- as.numeric(data$cost_mean)
data$new_id <- sub(".$", "", data$id)  # Standardise ID for grouping

# -----------------------------
# Visualise Cost Distribution
# -----------------------------
ggplot(data, aes(x = cost_mean)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
  geom_density(alpha = .2, fill = "#FF6666") +
  labs(title = "Distribution of Home Range", x = "Home Range", y = "Density")

hist(data$cost_mean)
plot(density(data$cost_mean))
lines(density(rnorm(1000, mean=mean(data$cost_mean), sd=sd(data$cost_mean))), col="blue")

# -----------------------------
# Model: All Landscapes
# -----------------------------
cost_model <- glmmTMB(cost_mean ~ landscape + (1|new_id), data = data, family = ordbeta())
summary(cost_model)

# -----------------------------
# Coefficient Extraction Function
# -----------------------------
extract_coefficients <- function(model, model_name) {
  model_summary <- summary(model)
  conditional_coefficients <- as.data.frame(model_summary$coeff$cond)
  conditional_coefficients$term <- rownames(conditional_coefficients)
  conditional_coefficients$conf.low <- conditional_coefficients$Estimate - 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$conf.high <- conditional_coefficients$Estimate + 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$color <- ifelse(conditional_coefficients$`Pr(>|z|)` < 0.05, "black", "grey")
  conditional_coefficients$model <- model_name
  return(conditional_coefficients)
}

# -----------------------------
# Coefficient Plot Function
# -----------------------------
plot_coefficients <- function(coefficients) {
  coefficients <- coefficients[coefficients$term != "(Intercept)", ]
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#FF9F5F", "gray"))
  
  plot <- ggplot(coefficients, aes(x = Estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(title = "Non-mining vs. mining analysis", x = "Energetic cost", y = "Landscape") +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5)
    )
  
  if ("model" %in% names(coefficients)) {
    plot <- plot + facet_grid(model ~ ., scales = "free_y", space = "free_y")
  }
  
  return(plot)
}

# -----------------------------
# Format Coefficients and Plot (All Landscapes)
# -----------------------------
cost_coeff <- extract_coefficients(cost_model, "landscape")

cost_coeff <- cost_coeff %>%
  mutate(term = case_when(
    term == "landscapeoriginal" ~ "Current mining",
    term == "landscapedense" ~ "Consolidated mining",
    term == "landscapefragmented" ~ "Fragmented mining",
    term == "(Intercept)" ~ "Non-mining",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Non-mining", "Consolidated mining", "Current mining", "Fragmented mining"))) %>%
  filter(term != "Non-mining")

cost_plot <- plot_coefficients(cost_coeff)
print(cost_plot)

# -----------------------------
# Model: Mining Landscapes Only
# -----------------------------
data_mining <- data %>% filter(landscape != "non_mining")
data_mining$landscape <- fct_relevel(data_mining$landscape, "original")

cost_model_mining <- glmmTMB(cost_mean ~ landscape + (1|new_id), data = data_mining, family = ordbeta())
summary(cost_model_mining)

# -----------------------------
# Format Coefficients and Plot (Mining Only)
# -----------------------------
cost_mining_coeff <- extract_coefficients(cost_model_mining, "landscape")

# Replace intercept values to 0 for reference
cost_mining_coeff$Estimate[cost_mining_coeff$term == "(Intercept)"] <- 0
cost_mining_coeff$conf.low[cost_mining_coeff$term == "(Intercept)"] <- 0
cost_mining_coeff$conf.high[cost_mining_coeff$term == "(Intercept)"] <- 0

cost_mining_coeff <- cost_mining_coeff %>%
  mutate(term = case_when(
    term == "landscapedense" ~ "Consolidated mining",
    term == "landscapefragmented" ~ "Fragmented mining",
    term == "(Intercept)" ~ "Current mining",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Current mining", "Consolidated mining", "Fragmented mining")))

cost_mining_plot <- plot_coefficients(cost_mining_coeff)
print(cost_mining_plot)


######################################
######## Calculate Revisitations #####
######################################
library(dplyr)
library(move)
library(recurse)
library(terra)
library(glmmTMB)
library(ggplot2)
library(viridis)
library(tictoc)

rm(list = ls())

# Read combined simulation data
setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Simulation trajectories")
data <- read.csv("all_simulation_trajectories_with_memory.csv")

# Define root path for NDVI rasters
ndvi_paths <- list(
  original = "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/current mining/Current_habitat_classifications.tif",
  fragmented = "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/dispersed mining/Dispersed_habitat_classifications.tif",
  dense = "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/aggregated mining/Aggregated habitat classification.tif",
  non_mining = "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Habitat maps/non_mining/Non_mining_habitat_classifications.tif"
)

# Output directory
output_dir <- "C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Revisitations"

# Radius for revisits
revisit_radius <- 20

# Loop through each landscape
for (land in names(ndvi_paths)) {
  cat("\nProcessing:", land, "\n")
  tic()
  
  # Filter data for current landscape
  land_data <- data %>% filter(landscape == land)
  
  # Format timestamps and factors
  land_data$t_ <- ifelse(nchar(as.character(land_data$t_)) <= 10,
                         paste0(land_data$t_, " 00:00:00"),
                         as.character(land_data$t_))
  land_data$t_ <- as.POSIXct(land_data$t_, format="%d/%m/%Y %H:%M")
  land_data$id <- as.factor(land_data$id)
  land_data <- land_data %>%
    arrange(id, t_) %>%
    mutate(t_ = as.POSIXct(t_, format="%Y-%m-%d %H:%M:%S"))
  
  # Load NDVI raster
  ndvi <- rast(ndvi_paths[[land]])
  proj_info <- as.character(crs(ndvi))
  
  # Prepare results container
  results_df <- data.frame()
  unique_ids <- unique(land_data$id)
  
  for (id in unique_ids) {
    cat("  ID:", id, "\n")
    
    # Subset and create move object
    subset_data <- land_data[land_data$id == id, ]
    subset_move <- move(x = subset_data$x_, y = subset_data$y_,
                        time = subset_data$t_, animal = subset_data$id,
                        proj = proj_info)
    
    # Calculate revisits and extract NDVI
    revisits <- getRecursions(subset_move, revisit_radius)
    revisits$habitat <- terra::extract(ndvi, as.data.frame(subset_move)[, c("x", "y")])[, 2]
    
    # Merge data
    subset_data$revisits <- revisits$revisits
    subset_data$habitat <- revisits$habitat
    
    # Classify habitat
    subset_data <- subset_data %>%
      mutate(habitat_cat = case_when(
        habitat == 1 ~ "dense_veg",
        habitat == 2 ~ "grassland",
        habitat == 3 ~ "rocky",
        habitat == 4 ~ "other_disturbed",
        habitat == 5 ~ "pits_dumps",
        TRUE         ~ as.character(habitat)
      ))
    
    # Set habitat factor levels
    subset_data$habitat_cat <- factor(subset_data$habitat_cat,
                                      levels = c("rocky", "grassland", "dense_veg", "pits_dumps", "other_disturbed"))
    
    if (land == "non_mining") {
      subset_data$habitat_cat <- factor(subset_data$habitat_cat,
                                        levels = c("rocky", "grassland", "dense_veg"))
    }
    
    subset_data$landscape <- land
    results_df <- rbind(results_df, subset_data)
    
    rm(subset_data, subset_move, revisits)
    gc()
  }
  
  # Write output
  # write.csv(results_df,
  #           file.path(output_dir, paste0(gsub("-", "_", land), "_revisitations.csv")),
  #           row.names = FALSE)
  toc()
}


# library(dplyr)
# library(readr)
# 
# # Set directory where all the revisit CSVs are stored
# setwd("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Data used for models")
# # List of landscape types
# landscapes <- c("original", "fragmented", "dense", "non_mining")
# # Read and combine
# all_revisits <- lapply(landscapes, function(land) {
#   file <- paste0(land, "_revisitations.csv")
#   read_csv(file)
# }) %>% bind_rows()
# # Write to a single CSV
# write_csv(all_revisits, "habitat_revisitations.csv")


####################################
######## Revisitation analysis #####
####################################
rm(list = ls())

data <- read.csv("C:/Users/00108247/OneDrive - UWA/Desktop/Mitch stuff/Simulation paper/Data and code/Figshare data/Data used for models/habitat_revisitations.csv")
head(data)
str(data)
unique(data$landscape)
unique(data$id)

data$revisits <- as.numeric(data$revisits)
data$new_id <- data$id
#remove last digit of id to make sure all iterations of that location are the same id
data$new_id <- sub(".$", "", data$new_id)
# Assuming you have data with revisits column

# 
###rocky
# Filter for "rocky" habitat and fit GLMM
rocky_data <- data %>% filter(habitat_cat == "rocky")
str(rocky_data)
# Load necessary library
library(ggplot2)
hist(rocky_data$revisits)
plot(density(rocky_data$revisits))
lines(density(rnorm(1000, mean=mean(rocky_data$revisits), sd=sd(rocky_data$revisits))), col="blue")
rocky_data$landscape <-factor(rocky_data$landscape, levels = c("non_mining", "fragmented", "dense","original"))
summary(rocky_data$revisits)
hist(rocky_data$revisits, breaks=50, main="Histogram of Revisitations", xlab="Number of Revisits")
mean_revisits <- mean(rocky_data$revisits)
var_revisits <- var(rocky_data$revisits)
print(paste("Mean: ", mean_revisits))
print(paste("Variance: ", var_revisits))
if(var_revisits > mean_revisits) {
  print("Data shows overdispersion; consider Negative Binomial.")
} else {
  print("Data may fit a Poisson distribution.")
}
library(glmmTMB)
rocky_data <- rocky_data %>%
  mutate(obs = row_number(rocky_data))
rocky_model <- glmmTMB(revisits ~ landscape + (1|new_id)+ (1|obs), data = rocky_data, family=poisson)
summary(rocky_model)

# Define the function to extract coefficients and confidence intervals
extract_coefficients <- function(model, model_name) {
  model_summary <- summary(model)
  conditional_coefficients <- as.data.frame(model_summary$coeff$cond)
  conditional_coefficients$term <- rownames(conditional_coefficients)
  
  # Calculate confidence intervals
  conditional_coefficients$conf.low <- conditional_coefficients$Estimate - 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$conf.high <- conditional_coefficients$Estimate + 1.96 * conditional_coefficients$`Std. Error`
  
  # Create a new column for color based on significance
  conditional_coefficients$color <- ifelse(conditional_coefficients$`Pr(>|z|)` < 0.05, "black", "grey")
  
  # Add model name as a new column
  conditional_coefficients$model <- model_name
  
  return(conditional_coefficients)
}



library(ggplot2)

plot_coefficients <- function(coefficients) {
  # Exclude the intercept from the plot
  coefficients <- coefficients[coefficients$term != "(Intercept)", ]
  
  # Set color based on confidence intervals and estimate
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#FF9F5F", "gray"))
  
  # Plot setup
  plot <- ggplot(coefficients, aes(x = Estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(
      title = "Non-mining vs. mining analysis", # Add your plot title here
      x = "Revisitations coefficient (β)", 
      y = "Landscape"
    ) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5) # Center the title
    )
  
  # If comparing models, add facets
  if ("model" %in% names(coefficients)) {
    plot <- plot + facet_grid(model ~ ., scales = "free_y", space = "free_y")
  }
  
  return(plot)
}

rocky_coeff <- extract_coefficients(rocky_model, "landscape")

library(dplyr)
library(forcats) # For fct_relevel

# Renaming and reordering the terms
rocky_coeff <- rocky_coeff %>%
  mutate(term = case_when(
    term == "landscapeoriginal" ~ "Current landscape",
    term == "landscapedense" ~ "Aggregated landscape",
    term == "landscapefragmented" ~ "Dispersed landscape",
    term == "(Intercept)" ~ "Non-mining landscape",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Non-mining landscape", "Aggregated landscape", "Current landscape", "Dispersed landscape")))

rocky_coeff <- rocky_coeff[rocky_coeff$term != "Non-mining landscape",]

# Assuming you've already extracted coefficients into `all_coeff`
rocky_plot <- plot_coefficients(rocky_coeff)

# Display the plot
print(rocky_plot)



###grassland
# Filter for "grassland" habitat and fit GLMM
grassland_data <- data %>% filter(habitat_cat == "grassland")
str(grassland_data)
# Load necessary library
library(ggplot2)
hist(grassland_data$revisits)
plot(density(grassland_data$revisits))
lines(density(rnorm(1000, mean=mean(grassland_data$revisits), sd=sd(grassland_data$revisits))), col="blue")
grassland_data$landscape <-factor(grassland_data$landscape, levels = c("non_mining", "fragmented", "dense","original"))
summary(grassland_data$revisits)
hist(grassland_data$revisits, breaks=50, main="Histogram of Revisitations", xlab="Number of Revisits")
mean_revisits <- mean(grassland_data$revisits)
var_revisits <- var(grassland_data$revisits)
print(paste("Mean: ", mean_revisits))
print(paste("Variance: ", var_revisits))
if(var_revisits > mean_revisits) {
  print("Data shows overdispersion; consider Negative Binomial.")
} else {
  print("Data may fit a Poisson distribution.")
}
# Assuming your data frame is named grassland_data
grassland_data <- grassland_data %>%
  mutate(obs = row_number(grassland_data))
library(glmmTMB)
grassland_model <- glmmTMB(revisits ~ landscape + (1|new_id) + (1|obs), data = grassland_data, family=poisson)
summary(grassland_model)


# Define the function to extract coefficients and confidence intervals
extract_coefficients <- function(model, model_name) {
  model_summary <- summary(model)
  conditional_coefficients <- as.data.frame(model_summary$coeff$cond)
  conditional_coefficients$term <- rownames(conditional_coefficients)
  
  # Calculate confidence intervals
  conditional_coefficients$conf.low <- conditional_coefficients$Estimate - 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$conf.high <- conditional_coefficients$Estimate + 1.96 * conditional_coefficients$`Std. Error`
  
  # Create a new column for color based on significance
  conditional_coefficients$color <- ifelse(conditional_coefficients$`Pr(>|z|)` < 0.05, "black", "grey")
  
  # Add model name as a new column
  conditional_coefficients$model <- model_name
  
  return(conditional_coefficients)
}



library(ggplot2)

plot_coefficients <- function(coefficients) {
  # Exclude the intercept from the plot
  coefficients <- coefficients[coefficients$term != "(Intercept)", ]
  
  # Set color based on confidence intervals and estimate
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#FF9F5F", "gray"))
  
  # Plot setup
  plot <- ggplot(coefficients, aes(x = Estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(
      title = "Non-mining vs. mining analysis", # Add your plot title here
      x = "Revisitations coefficient (β)", 
      y = "Landscape"
    ) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5) # Center the title
    )
  
  # If comparing models, add facets
  if ("model" %in% names(coefficients)) {
    plot <- plot + facet_grid(model ~ ., scales = "free_y", space = "free_y")
  }
  
  return(plot)
}

grassland_coeff <- extract_coefficients(grassland_model, "landscape")

library(dplyr)
library(forcats) # For fct_relevel

# Renaming and reordering the terms
grassland_coeff <- grassland_coeff %>%
  mutate(term = case_when(
    term == "landscapeoriginal" ~ "Current landscape",
    term == "landscapedense" ~ "Aggregated landscape",
    term == "landscapefragmented" ~ "Dispersed landscape",
    term == "(Intercept)" ~ "Non-mining landscape",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Non-mining landscape", "Aggregated landscape", "Current landscape", "Dispersed landscape")))

grassland_coeff <- grassland_coeff[grassland_coeff$term != "Non-mining landscape",]

# Assuming you've already extracted coefficients into `all_coeff`
grassland_plot <- plot_coefficients(grassland_coeff)

# Display the plot
print(grassland_plot)






###dense_veg
# Filter for "dense_veg" habitat and fit GLMM
dense_veg_data <- data %>% filter(habitat_cat == "dense_veg")
str(dense_veg_data)
# Load necessary library
library(ggplot2)
hist(dense_veg_data$revisits)
plot(density(dense_veg_data$revisits))
lines(density(rnorm(1000, mean=mean(dense_veg_data$revisits), sd=sd(dense_veg_data$revisits))), col="blue")
dense_veg_data$landscape <-factor(dense_veg_data$landscape, levels = c("non_mining", "fragmented", "dense","original"))
summary(dense_veg_data$revisits)
hist(dense_veg_data$revisits, breaks=50, main="Histogram of Revisitations", xlab="Number of Revisits")
mean_revisits <- mean(dense_veg_data$revisits)
var_revisits <- var(dense_veg_data$revisits)
print(paste("Mean: ", mean_revisits))
print(paste("Variance: ", var_revisits))
if(var_revisits > mean_revisits) {
  print("Data shows overdispersion; consider Negative Binomial.")
} else {
  print("Data may fit a Poisson distribution.")
}
# Assuming your data frame is named dense_veg_data
dense_veg_data <- dense_veg_data %>%
  mutate(obs = row_number(dense_veg_data))
library(glmmTMB)
dense_veg_model <- glmmTMB(revisits ~ landscape + (1|new_id) + (1|obs), data = dense_veg_data, family=poisson)
summary(dense_veg_model)


# Define the function to extract coefficients and confidence intervals
extract_coefficients <- function(model, model_name) {
  model_summary <- summary(model)
  conditional_coefficients <- as.data.frame(model_summary$coeff$cond)
  conditional_coefficients$term <- rownames(conditional_coefficients)
  
  # Calculate confidence intervals
  conditional_coefficients$conf.low <- conditional_coefficients$Estimate - 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$conf.high <- conditional_coefficients$Estimate + 1.96 * conditional_coefficients$`Std. Error`
  
  # Create a new column for color based on significance
  conditional_coefficients$color <- ifelse(conditional_coefficients$`Pr(>|z|)` < 0.05, "black", "grey")
  
  # Add model name as a new column
  conditional_coefficients$model <- model_name
  
  return(conditional_coefficients)
}



library(ggplot2)

plot_coefficients <- function(coefficients) {
  # Exclude the intercept from the plot
  coefficients <- coefficients[coefficients$term != "(Intercept)", ]
  
  # Set color based on confidence intervals and estimate
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#FF9F5F", "gray"))
  
  # Plot setup
  plot <- ggplot(coefficients, aes(x = Estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(
      title = "Non-mining vs. mining analysis", # Add your plot title here
      x = "Revisitations coefficient (β)", 
      y = "Landscape"
    ) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5) # Center the title
    )
  
  # If comparing models, add facets
  if ("model" %in% names(coefficients)) {
    plot <- plot + facet_grid(model ~ ., scales = "free_y", space = "free_y")
  }
  
  return(plot)
}

dense_veg_coeff <- extract_coefficients(dense_veg_model, "landscape")

library(dplyr)
library(forcats) # For fct_relevel

# Renaming and reordering the terms
dense_veg_coeff <- dense_veg_coeff %>%
  mutate(term = case_when(
    term == "landscapeoriginal" ~ "Current landscape",
    term == "landscapedense" ~ "Aggregated landscape",
    term == "landscapefragmented" ~ "Dispersed landscape",
    term == "(Intercept)" ~ "Non-mining landscape",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Non-mining landscape", "Aggregated landscape", "Current landscape", "Dispersed landscape")))

dense_veg_coeff <- dense_veg_coeff[dense_veg_coeff$term != "Non-mining landscape",]

# Assuming you've already extracted coefficients into `all_coeff`
dense_veg_plot <- plot_coefficients(dense_veg_coeff)

# Display the plot
print(dense_veg_plot)






###PLOT all
library(dplyr)
library(ggplot2)
library(forcats) # For fct_relevel

# Assuming dense_veg_coeff, grassland_coeff, and rocky_coeff are already loaded and have the correct structure

# Step 1: Add habitat type and combine data frames
dense_veg_coeff$habitat_type <- 'Riparian habitat'
grassland_coeff$habitat_type <- 'Spinifex grassland'
rocky_coeff$habitat_type <- 'Rocky habitat'

# Combine into one data frame
combined_coeff <- rbind(dense_veg_coeff, grassland_coeff, rocky_coeff)

# Assuming combined_coeff is already loaded and contains the correct structure with the 'term' column

# Set the order of landscapes
combined_coeff <- combined_coeff %>%
  mutate(landscape = fct_relevel(term,
                                 "Current landscape",
                                 "Dispersed landscape",
                                 "Aggregated landscape")) %>%
  mutate(landscape = fct_recode(landscape,
                                "Dispersed mining" = "Dispersed landscape",
                                "Current mining" = "Current landscape",
                                "Aggregated mining" = "Aggregated landscape"))

# Set the order of habitat types
combined_coeff <- combined_coeff %>%
  mutate(habitat_type = fct_relevel(habitat_type,
                                    "Rocky habitat",
                                    "Spinifex grassland",
                                    "Riparian habitat"))

# Continue with the rest of the plotting code


# Check the updated data frame
print(combined_coeff)


# Updated function to format the faceted plots
plot_coefficients_horizontal_facets <- function(coefficients) {
  # Set color based on confidence intervals and estimate
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#440154FF", "gray"))
  
  # Plot setup
  plot <- ggplot(coefficients, aes(x = Estimate, y = habitat_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(
      title = "Non-mining vs. mining model",
      x = "Coefficient (β)",
      y = "Revisitations"
    ) +
    facet_wrap(~landscape, ncol = 3) + # Place facets side by side
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5), # Center the title
      strip.background = element_rect(fill = "white"),
      strip.text.x = element_text(size = 14, color = "black", face = "bold")
    )
  
  return(plot)
}

# Create the formatted horizontally faceted plot
formatted_horizontal_faceted_plot <- plot_coefficients_horizontal_facets(combined_coeff)

# Print the formatted plot
print(formatted_horizontal_faceted_plot)
























###################
# Mining version revisitations
########################
#rocky
#Filter for "rocky" habitat and fit GLMM
rocky_data <- data %>% filter(habitat_cat == "rocky")
str(rocky_data)
# Remove 'non_mining' from the landscape column
rocky_mining_data <- rocky_data %>%
  filter(landscape != "non_mining")
# Set 'original' as the reference category for the landscape variable
rocky_mining_data$landscape <- fct_relevel(rocky_mining_data$landscape, "original")

# Load necessary library
library(ggplot2)
hist(rocky_mining_data$revisits)
plot(density(rocky_mining_data$revisits))
lines(density(rnorm(1000, mean=mean(rocky_mining_data$revisits), sd=sd(rocky_mining_data$revisits))), col="blue")
rocky_mining_data$landscape <-factor(rocky_mining_data$landscape, levels = c("original", "fragmented", "dense"))
summary(rocky_mining_data$revisits)
hist(rocky_mining_data$revisits, breaks=50, main="Histogram of Revisitations", xlab="Number of Revisits")
mean_revisits <- mean(rocky_mining_data$revisits)
var_revisits <- var(rocky_mining_data$revisits)
print(paste("Mean: ", mean_revisits))
print(paste("Variance: ", var_revisits))
if(var_revisits > mean_revisits) {
  print("Data shows overdispersion; consider Negative Binomial.")
} else {
  print("Data may fit a Poisson distribution.")
}
library(glmmTMB)
rocky_mining_data <- rocky_mining_data %>%
  mutate(obs = row_number(rocky_mining_data))
rocky_mining_model <- glmmTMB(revisits ~ landscape + (1|new_id)+ (1|obs), data = rocky_mining_data, family=poisson)
summary(rocky_mining_model)

# Define the function to extract coefficients and confidence intervals
extract_coefficients <- function(model, model_name) {
  model_summary <- summary(model)
  conditional_coefficients <- as.data.frame(model_summary$coeff$cond)
  conditional_coefficients$term <- rownames(conditional_coefficients)
  
  # Calculate confidence intervals
  conditional_coefficients$conf.low <- conditional_coefficients$Estimate - 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$conf.high <- conditional_coefficients$Estimate + 1.96 * conditional_coefficients$`Std. Error`
  
  # Create a new column for color based on significance
  conditional_coefficients$color <- ifelse(conditional_coefficients$`Pr(>|z|)` < 0.05, "black", "grey")
  
  # Add model name as a new column
  conditional_coefficients$model <- model_name
  
  return(conditional_coefficients)
}



library(ggplot2)

plot_coefficients <- function(coefficients) {
  # Exclude the intercept from the plot
  coefficients <- coefficients[coefficients$term != "(Intercept)", ]
  
  # Set color based on confidence intervals and estimate
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#FF9F5F", "gray"))
  
  # Plot setup
  plot <- ggplot(coefficients, aes(x = Estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(
      title = "Mining composition analysis", # Add your plot title here
      x = "Coefficient (β)", 
      y = "Revisitations"
    ) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5) # Center the title
    )
  
  # If comparing models, add facets
  if ("model" %in% names(coefficients)) {
    plot <- plot + facet_grid(model ~ ., scales = "free_y", space = "free_y")
  }
  
  return(plot)
}

rocky_mining_coeff <- extract_coefficients(rocky_mining_model, "landscape")



# Identify the row with the intercept
intercept_row <- which(rocky_mining_coeff$term == "(Intercept)")
# Set the Estimate of the intercept to 0
rocky_mining_coeff$Estimate[intercept_row] <- 0
# Identify the row with the intercept
intercept_row <- which(rocky_mining_coeff$term == "(Intercept)")
# Set the Estimate of the intercept to 0
rocky_mining_coeff$conf.low[intercept_row] <- 0
# Identify the row with the intercept
intercept_row <- which(rocky_mining_coeff$term == "(Intercept)")
# Set the Estimate of the intercept to 0
rocky_mining_coeff$conf.high[intercept_row] <- 0


library(dplyr)
library(forcats) # For fct_relevel

# Renaming and reordering the terms
rocky_mining_coeff <- rocky_mining_coeff %>%
  mutate(term = case_when(
    term == "landscapedense" ~ "Aggregated mining",
    term == "landscapefragmented" ~ "Dispersed mining",
    term == "(Intercept)" ~ "Current mining",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Current mining","Aggregated mining", "Dispersed mining")))

# rocky_mining_coeff <- rocky_mining_coeff[rocky_mining_coeff$term != "Current mining",]

# Assuming you've already extracted coefficients into `all_coeff`
rocky_mining_plot <- plot_coefficients(rocky_mining_coeff)

# Display the plot
print(rocky_mining_plot)








###grassland
# Filter for "grassland" habitat and fit GLMM
grassland_data <- data %>% filter(habitat_cat == "grassland")
str(grassland_data)
# Remove 'non_mining' from the landscape column
grassland_mining_data <- grassland_data %>%
  filter(landscape != "non_mining")
# Set 'original' as the reference category for the landscape variable
grassland_mining_data$landscape <- fct_relevel(grassland_mining_data$landscape, "original")
# Load necessary library
library(ggplot2)
hist(grassland_mining_data$revisits)
plot(density(grassland_mining_data$revisits))
lines(density(rnorm(1000, mean=mean(grassland_mining_data$revisits), sd=sd(grassland_mining_data$revisits))), col="blue")
grassland_mining_data$landscape <-factor(grassland_mining_data$landscape, levels = c("original", "fragmented", "dense"))
summary(grassland_mining_data$revisits)
hist(grassland_mining_data$revisits, breaks=50, main="Histogram of Revisitations", xlab="Number of Revisits")
mean_revisits <- mean(grassland_mining_data$revisits)
var_revisits <- var(grassland_mining_data$revisits)
print(paste("Mean: ", mean_revisits))
print(paste("Variance: ", var_revisits))
if(var_revisits > mean_revisits) {
  print("Data shows overdispersion; consider Negative Binomial.")
} else {
  print("Data may fit a Poisson distribution.")
}
# Assuming your data frame is named grassland_mining_data
grassland_mining_data <- grassland_mining_data %>%
  mutate(obs = row_number(grassland_mining_data))
library(glmmTMB)
grassland_mining_model <- glmmTMB(revisits ~ landscape + (1|new_id) + (1|obs), data = grassland_mining_data, family=poisson)
summary(grassland_mining_model)


# Define the function to extract coefficients and confidence intervals
extract_coefficients <- function(model, model_name) {
  model_summary <- summary(model)
  conditional_coefficients <- as.data.frame(model_summary$coeff$cond)
  conditional_coefficients$term <- rownames(conditional_coefficients)
  
  # Calculate confidence intervals
  conditional_coefficients$conf.low <- conditional_coefficients$Estimate - 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$conf.high <- conditional_coefficients$Estimate + 1.96 * conditional_coefficients$`Std. Error`
  
  # Create a new column for color based on significance
  conditional_coefficients$color <- ifelse(conditional_coefficients$`Pr(>|z|)` < 0.05, "black", "grey")
  
  # Add model name as a new column
  conditional_coefficients$model <- model_name
  
  return(conditional_coefficients)
}



library(ggplot2)

plot_coefficients <- function(coefficients) {
  # Exclude the intercept from the plot
  coefficients <- coefficients[coefficients$term != "(Intercept)", ]
  
  # Set color based on confidence intervals and estimate
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#FF9F5F", "gray"))
  
  # Plot setup
  plot <- ggplot(coefficients, aes(x = Estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(
      title = "Non-mining vs. mining analysis", # Add your plot title here
      x = "Revisitations coefficient (β)", 
      y = "Landscape"
    ) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5) # Center the title
    )
  
  # If comparing models, add facets
  if ("model" %in% names(coefficients)) {
    plot <- plot + facet_grid(model ~ ., scales = "free_y", space = "free_y")
  }
  
  return(plot)
}

grassland_mining_coeff <- extract_coefficients(grassland_mining_model, "landscape")

library(dplyr)
library(forcats) # For fct_relevel

# Renaming and reordering the terms
grassland_mining_coeff <- grassland_mining_coeff %>%
  mutate(term = case_when(
    term == "landscapeoriginal" ~ "Current mining",
    term == "landscapedense" ~ "Aggregated mining",
    term == "landscapefragmented" ~ "Dispersed mining",
    term == "(Intercept)" ~ "Non-mining",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Non-mining", "Aggregated mining", "Current mining", "Dispersed mining")))

grassland_mining_coeff <- grassland_mining_coeff[grassland_mining_coeff$term != "Non-mining",]

# Assuming you've already extracted coefficients into `all_coeff`
grassland_mining_plot <- plot_coefficients(grassland_mining_coeff)

# Display the plot
print(grassland_mining_plot)










##dense_veg
#Filter for "dense_veg" habitat and fit GLMM
dense_veg_data <- data %>% filter(habitat_cat == "dense_veg")
str(dense_veg_data)
# Remove 'non_mining' from the landscape column
dense_veg_mining_data <- dense_veg_data %>%
  filter(landscape != "non_mining")
# Set 'original' as the reference category for the landscape variable
dense_veg_mining_data$landscape <- fct_relevel(dense_veg_mining_data$landscape, "original")
# Load necessary library
library(ggplot2)
hist(dense_veg_mining_data$revisits)
plot(density(dense_veg_mining_data$revisits))
lines(density(rnorm(1000, mean=mean(dense_veg_mining_data$revisits), sd=sd(dense_veg_mining_data$revisits))), col="blue")
dense_veg_mining_data$landscape <-factor(dense_veg_mining_data$landscape, levels = c("original", "fragmented", "dense"))
summary(dense_veg_mining_data$revisits)
hist(dense_veg_mining_data$revisits, breaks=50, main="Histogram of Revisitations", xlab="Number of Revisits")
mean_revisits <- mean(dense_veg_mining_data$revisits)
var_revisits <- var(dense_veg_mining_data$revisits)
print(paste("Mean: ", mean_revisits))
print(paste("Variance: ", var_revisits))
if(var_revisits > mean_revisits) {
  print("Data shows overdispersion; consider Negative Binomial.")
} else {
  print("Data may fit a Poisson distribution.")
}
# Assuming your data frame is named dense_veg_mining_data
dense_veg_mining_data <- dense_veg_mining_data %>%
  mutate(obs = row_number(dense_veg_mining_data))
library(glmmTMB)
dense_veg_mining_model <- glmmTMB(revisits ~ landscape + (1|new_id) + (1|obs), data = dense_veg_mining_data, family=poisson)
summary(dense_veg_mining_model)


# Define the function to extract coefficients and confidence intervals
extract_coefficients <- function(model, model_name) {
  model_summary <- summary(model)
  conditional_coefficients <- as.data.frame(model_summary$coeff$cond)
  conditional_coefficients$term <- rownames(conditional_coefficients)
  
  # Calculate confidence intervals
  conditional_coefficients$conf.low <- conditional_coefficients$Estimate - 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$conf.high <- conditional_coefficients$Estimate + 1.96 * conditional_coefficients$`Std. Error`
  
  # Create a new column for color based on significance
  conditional_coefficients$color <- ifelse(conditional_coefficients$`Pr(>|z|)` < 0.05, "black", "grey")
  
  # Add model name as a new column
  conditional_coefficients$model <- model_name
  
  return(conditional_coefficients)
}



library(ggplot2)

plot_coefficients <- function(coefficients) {
  # Exclude the intercept from the plot
  coefficients <- coefficients[coefficients$term != "(Intercept)", ]
  
  # Set color based on confidence intervals and estimate
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#FF9F5F", "gray"))
  
  # Plot setup
  plot <- ggplot(coefficients, aes(x = Estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(
      title = "Non-mining vs. mining analysis", # Add your plot title here
      x = "Revisitations coefficient (β)", 
      y = "Landscape"
    ) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5) # Center the title
    )
  
  # If comparing models, add facets
  if ("model" %in% names(coefficients)) {
    plot <- plot + facet_grid(model ~ ., scales = "free_y", space = "free_y")
  }
  
  return(plot)
}

dense_veg_mining_coeff <- extract_coefficients(dense_veg_mining_model, "landscape")

library(dplyr)
library(forcats) # For fct_relevel

# Renaming and reordering the terms
dense_veg_mining_coeff <- dense_veg_mining_coeff %>%
  mutate(term = case_when(
    term == "landscapeoriginal" ~ "Current mining",
    term == "landscapedense" ~ "Aggregated mining",
    term == "landscapefragmented" ~ "Dispersed mining",
    term == "(Intercept)" ~ "Non-mining",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Non-mining", "Aggregated mining", "Current mining", "Dispersed mining")))

dense_veg_mining_coeff <- dense_veg_mining_coeff[dense_veg_mining_coeff$term != "Non-mining",]

# Assuming you've already extracted coefficients into `all_coeff`
dense_veg_mining_plot <- plot_coefficients(dense_veg_mining_coeff)

# Display the plot
print(dense_veg_mining_plot)









###pits_dumps
# Filter for "pits_dumps" habitat and fit GLMM
pits_dumps_data <- data %>% filter(habitat_cat == "pits_dumps")
str(pits_dumps_data)
# Load necessary library
library(ggplot2)
hist(pits_dumps_data$revisits)
plot(density(pits_dumps_data$revisits))
lines(density(rnorm(1000, mean=mean(pits_dumps_data$revisits), sd=sd(pits_dumps_data$revisits))), col="blue")
pits_dumps_data$landscape <-factor(pits_dumps_data$landscape, levels = c("original", "fragmented", "dense"))
summary(pits_dumps_data$revisits)
hist(pits_dumps_data$revisits, breaks=50, main="Histogram of Revisitations", xlab="Number of Revisits")
mean_revisits <- mean(pits_dumps_data$revisits)
var_revisits <- var(pits_dumps_data$revisits)
print(paste("Mean: ", mean_revisits))
print(paste("Variance: ", var_revisits))
if(var_revisits > mean_revisits) {
  print("Data shows overdispersion; consider Negative Binomial.")
} else {
  print("Data may fit a Poisson distribution.")
}
# Assuming your data frame is named pits_dumps_data
pits_dumps_data <- pits_dumps_data %>%
  mutate(obs = row_number(pits_dumps_data))
library(glmmTMB)
pits_dumps_model <- glmmTMB(revisits ~ landscape + (1|new_id) + (1|obs), data = pits_dumps_data, family=poisson)
summary(pits_dumps_model)


# Define the function to extract coefficients and confidence intervals
extract_coefficients <- function(model, model_name) {
  model_summary <- summary(model)
  conditional_coefficients <- as.data.frame(model_summary$coeff$cond)
  conditional_coefficients$term <- rownames(conditional_coefficients)
  
  # Calculate confidence intervals
  conditional_coefficients$conf.low <- conditional_coefficients$Estimate - 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$conf.high <- conditional_coefficients$Estimate + 1.96 * conditional_coefficients$`Std. Error`
  
  # Create a new column for color based on significance
  conditional_coefficients$color <- ifelse(conditional_coefficients$`Pr(>|z|)` < 0.05, "black", "grey")
  
  # Add model name as a new column
  conditional_coefficients$model <- model_name
  
  return(conditional_coefficients)
}



library(ggplot2)

plot_coefficients <- function(coefficients) {
  # Exclude the intercept from the plot
  coefficients <- coefficients[coefficients$term != "(Intercept)", ]
  
  # Set color based on confidence intervals and estimate
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#FF9F5F", "gray"))
  
  # Plot setup
  plot <- ggplot(coefficients, aes(x = Estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(
      title = "Mining composition analysis", # Add your plot title here
      x = "Revisitations coefficient (β)", 
      y = "Landscape"
    ) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5) # Center the title
    )
  
  # If comparing models, add facets
  if ("model" %in% names(coefficients)) {
    plot <- plot + facet_grid(model ~ ., scales = "free_y", space = "free_y")
  }
  
  return(plot)
}

pits_dumps_coeff <- extract_coefficients(pits_dumps_model, "landscape")

library(dplyr)
library(forcats) # For fct_relevel

# Renaming and reordering the terms
pits_dumps_coeff <- pits_dumps_coeff %>%
  mutate(term = case_when(
    term == "landscapeoriginal" ~ "Current mining",
    term == "landscapedense" ~ "Aggregated mining",
    term == "landscapefragmented" ~ "Dispersed mining",
    term == "(Intercept)" ~ "Non-mining",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Non-mining", "Aggregated mining", "Current mining", "Dispersed mining")))

pits_dumps_coeff <- pits_dumps_coeff[pits_dumps_coeff$term != "Non-mining",]

# Assuming you've already extracted coefficients into `all_coeff`
pits_dumps_plot <- plot_coefficients(pits_dumps_coeff)

# Display the plot
print(pits_dumps_plot)












###other_disturbed
#Filter for "other_disturbed" habitat and fit GLMM
other_disturbed_data <- data %>% filter(habitat_cat == "other_disturbed")
str(other_disturbed_data)
# Load necessary library
library(ggplot2)
hist(other_disturbed_data$revisits)
plot(density(other_disturbed_data$revisits))
lines(density(rnorm(1000, mean=mean(other_disturbed_data$revisits), sd=sd(other_disturbed_data$revisits))), col="blue")
other_disturbed_data$landscape <-factor(other_disturbed_data$landscape, levels = c("original", "fragmented", "dense"))
summary(other_disturbed_data$revisits)
hist(other_disturbed_data$revisits, breaks=50, main="Histogram of Revisitations", xlab="Number of Revisits")
mean_revisits <- mean(other_disturbed_data$revisits)
var_revisits <- var(other_disturbed_data$revisits)
print(paste("Mean: ", mean_revisits))
print(paste("Variance: ", var_revisits))
if(var_revisits > mean_revisits) {
  print("Data shows overdispersion; consider Negative Binomial.")
} else {
  print("Data may fit a Poisson distribution.")
}
# Assuming your data frame is named other_disturbed_data
other_disturbed_data <- other_disturbed_data %>%
  mutate(obs = row_number(other_disturbed_data))
library(glmmTMB)
other_disturbed_model <- glmmTMB(revisits ~ landscape + (1|new_id) + (1|obs), data = other_disturbed_data, family=poisson)
summary(other_disturbed_model)


# Define the function to extract coefficients and confidence intervals
extract_coefficients <- function(model, model_name) {
  model_summary <- summary(model)
  conditional_coefficients <- as.data.frame(model_summary$coeff$cond)
  conditional_coefficients$term <- rownames(conditional_coefficients)
  
  # Calculate confidence intervals
  conditional_coefficients$conf.low <- conditional_coefficients$Estimate - 1.96 * conditional_coefficients$`Std. Error`
  conditional_coefficients$conf.high <- conditional_coefficients$Estimate + 1.96 * conditional_coefficients$`Std. Error`
  
  # Create a new column for color based on significance
  conditional_coefficients$color <- ifelse(conditional_coefficients$`Pr(>|z|)` < 0.05, "black", "grey")
  
  # Add model name as a new column
  conditional_coefficients$model <- model_name
  
  return(conditional_coefficients)
}



library(ggplot2)

plot_coefficients <- function(coefficients) {
  # Exclude the intercept from the plot
  coefficients <- coefficients[coefficients$term != "(Intercept)", ]
  
  # Set color based on confidence intervals and estimate
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#FF9F5F", "gray"))
  
  # Plot setup
  plot <- ggplot(coefficients, aes(x = Estimate, y = term)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(
      title = "Mining composition analysis", # Add your plot title here
      x = "Revisitations coefficient (β)", 
      y = "Landscape"
    ) +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5) # Center the title
    )
  
  # If comparing models, add facets
  if ("model" %in% names(coefficients)) {
    plot <- plot + facet_grid(model ~ ., scales = "free_y", space = "free_y")
  }
  
  return(plot)
}

other_disturbed_coeff <- extract_coefficients(other_disturbed_model, "landscape")

library(dplyr)
library(forcats) # For fct_relevel

# Renaming and reordering the terms
other_disturbed_coeff <- other_disturbed_coeff %>%
  mutate(term = case_when(
    term == "landscapeoriginal" ~ "Current mining",
    term == "landscapedense" ~ "Aggregated mining",
    term == "landscapefragmented" ~ "Dispersed mining",
    term == "(Intercept)" ~ "Non-mining",
    TRUE ~ term
  )) %>%
  mutate(term = factor(term, levels = c("Non-mining", "Aggregated mining", "Current mining", "Dispersed mining")))

other_disturbed_coeff <- other_disturbed_coeff[other_disturbed_coeff$term != "Non-mining",]

# Assuming you've already extracted coefficients into `all_coeff`
other_disturbed_plot <- plot_coefficients(other_disturbed_coeff)

# Display the plot
print(other_disturbed_plot)







###PLOT all
library(dplyr)
library(ggplot2)
library(forcats) # For fct_relevel

# Assuming dense_veg_coeff, grassland_coeff, and rocky_coeff are already loaded and have the correct structure

# Step 1: Add habitat type and combine data frames
dense_veg_mining_coeff$habitat_type <- 'Riparian habitat'
grassland_mining_coeff$habitat_type <- 'Spinifex grassland'
rocky_mining_coeff$habitat_type <- 'Rocky habitat'
pits_dumps_coeff$habitat_type <- 'Mine pits and waste dumps'
other_disturbed_coeff$habitat_type <- 'Other disturbed land'




# Combine into one data frame
combined_coeff <- rbind(dense_veg_mining_coeff, grassland_mining_coeff, rocky_mining_coeff,pits_dumps_coeff,other_disturbed_coeff)

# Assuming combined_coeff is already loaded and contains the correct structure with the 'term' column

# Set the order of landscapes
combined_coeff <- combined_coeff %>%
  mutate(landscape = fct_relevel(term,
                                 
                                 "Dispersed mining",
                                 "Aggregated mining",
                                 "Current mining"))

# Set the order of habitat types
combined_coeff <- combined_coeff %>%
  mutate(habitat_type = fct_relevel(habitat_type,
                                    "Rocky habitat",
                                    "Spinifex grassland",
                                    "Riparian habitat",
                                    "Mine pits and waste dumps",
                                    "Other disturbed land"))

# Continue with the rest of the plotting code


# Check the updated data frame
print(combined_coeff)


# Updated function to format the faceted plots
library(scales) # For additional formatting options

plot_coefficients_horizontal_facets <- function(coefficients) {
  # Set color based on confidence intervals and estimate
  coefficients$color <- ifelse(coefficients$conf.low > 0, "#9CC065",
                               ifelse(coefficients$conf.high < 0, "#440154FF", "gray"))
  
  # Plot setup
  plot <- ggplot(coefficients, aes(x = Estimate, y = habitat_type)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "red", linewidth = 2) +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high, color = color), height = 0.3, linewidth = 2) +
    geom_point(size = 5, aes(color = color)) +
    scale_color_identity() +
    theme_minimal() +
    labs(
      title = "Mining composition analysis (Intercept: Current mining landscape)",
      x = "Coefficient (β)",
      y = "Revisitations"
    ) +
    facet_wrap(~landscape, ncol = 3) + # Place facets side by side
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(color = "black", linewidth = 1.5),
      text = element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14),
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20)),
      plot.title = element_text(hjust = 0.5), # Center the title
      strip.background = element_rect(fill = "white"),
      strip.text.x = element_text(size = 14, color = "black", face = "bold")
    ) +
    scale_x_continuous(labels = function(x) format(x, nsmall = 2)) # Format x-axis labels to 2 decimal places
  
  return(plot)
}

# Assuming combined_coeff is prepared and valid for plotting as per previous steps
formatted_horizontal_faceted_plot <- plot_coefficients_horizontal_facets(combined_coeff)

# Print the formatted plot
print(formatted_horizontal_faceted_plot)


