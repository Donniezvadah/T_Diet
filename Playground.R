

# Load required libraries with tidyverse for data manipulation and visualization
library(tidyverse)    # Includes dplyr, tidyr, ggplot2, etc.
library(hmmTMB)       # For fitting hidden Markov models
library(readxl)       # For reading Excel files

# Clear workspace
rm(list = ls())

# ===============================================================================
# DATA LOADING AND PREPROCESSING
# ===============================================================================

# Load and inspect the dataset
cat("Loading and preparing nutrition data...\n")
nutrition_data <- read_excel("sampled_data_nutrition.xlsx")

# Display dataset information
cat("\n=== Dataset Information ===\n")
cat("Dimensions:", nrow(nutrition_data), "rows x", ncol(nutrition_data), "columns\n")
cat("\nFirst 3 rows:\n")
print(utils::head(nutrition_data, 3))

# Helper function to find columns by pattern
find_columns <- function(patterns, data = nutrition_data) {
  matches <- lapply(patterns, function(p) {
    grep(p, names(data), ignore.case = TRUE, value = TRUE)
  })
  unlist(Filter(length, matches))
}

# Set ID and Visit columns directly
id_col <- "IDNumber"
time_col <- "Visit"

# Create IDNumber and Visit columns if they don't exist
if (!id_col %in% names(nutrition_data)) {
  nutrition_data[[id_col]] <- seq_len(nrow(nutrition_data))
  cat("Created", id_col, "column as row identifier\n")
}

if (!time_col %in% names(nutrition_data)) {
  nutrition_data[[time_col]] <- seq_len(nrow(nutrition_data))
  cat("Created", time_col, "column as time indicator\n")
}

# Identify food group columns (variables 56-73 or by pattern)
food_group_cols <- find_columns(c("FG.*Percentage", "Percentage.*FG", "Food.*Group", "Energy.*Contribution"))

# Fallback to position-based selection if no columns found by pattern
if (length(food_group_cols) == 0) {
  food_group_cols <- if (ncol(nutrition_data) >= 73) {
    names(nutrition_data)[56:73]
  } else {
    tail(names(nutrition_data), 18)
  }
}

# Create final dataset with complete cases only
hmm_data <- nutrition_data %>% 
  select(all_of(c(id_col, time_col, food_group_cols))) %>%
  filter(complete.cases(.))

# Convert food group columns to numeric and scale
hmm_data <- hmm_data %>%
  mutate(across(all_of(food_group_cols), ~ as.numeric(as.character(.)))) %>%
  filter(complete.cases(.)) %>%
  mutate(across(all_of(food_group_cols), ~ scale(.)[,1]))

# Display final dataset information
cat("\n=== Processed Dataset ===\n")
cat("Final dimensions:", nrow(hmm_data), "rows x", ncol(hmm_data), "columns\n")
cat("\nSelected columns:", paste(names(hmm_data), collapse = ", "), "\n")

# Summary of food group data
cat("\n=== Food Group Summary ===\n")
print(summary(hmm_data[food_group_cols]))

# ===============================================================================
# HMM MODEL FITTING
# ===============================================================================

# Prepare data for hmmTMB
cat("\n=== Preparing HMM Model ===\n")

# Create observation matrix
obs_matrix <- as.matrix(hmm_data[food_group_cols])

# Basic HMM with 2 states (can be adjusted)
cat("Fitting HMM with 2 states...\n")
basic_hmm <- hmmTMB(
  obs = obs_matrix,
  data = hmm_data,
  n_states = 2,
  formula = ~ 1,
  family = list(gaussian())
)

# Display model summary
if (!is.null(basic_hmm)) {
  cat("\n=== HMM Model Summary ===\n")
  print(summary(basic_hmm))
  
  # Extract and display state probabilities
  states <- viterbi(basic_hmm)
  cat("\nState assignments:")
  print(table(states$state))
  
  # Plot state probabilities
  plot_states <- ggplot(states, aes(x = .data[[time_col]], y = .data[[id_col]], fill = factor(state))) +
    geom_tile() +
    labs(title = "HMM State Assignments",
         x = "Time",
         y = "Subject ID",
         fill = "State") +
    theme_minimal()
  
  print(plot_states)
  
} else {
  cat("Failed to fit HMM model\n")
}

cat("\nAnalysis complete.\n")
