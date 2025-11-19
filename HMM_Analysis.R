# ===============================================================================
# HIDDEN MARKOV MODEL ANALYSIS USING hmmTMB PACKAGE
# ===============================================================================
# Purpose: Analyze nutrition data using food group percentage energy contributions
# Author: Analysis of T_Diet dataset
# Date: November 2025
# ===============================================================================

# Load required libraries --------------------------------------------------------
library(dplyr)         # Data manipulation and filtering
library(hmmTMB)        # Main package for fitting hidden Markov models
library(ggplot2)       # Data visualization
library(readxl)        # Reading Excel files
library(tidyr)         # Data tidying
library(magrittr)      # Pipe operators

# Clear workspace to ensure clean environment -----------------------------------
rm(list = ls())

# ===============================================================================
# DATA LOADING SECTION
# ===============================================================================

# Read the main nutrition dataset -----------------------------------------------
# This dataset contains nutritional information with food group energy contributions
# Variables 56-73 represent percentage energy contributions from 18 food groups
cat("Loading nutrition dataset...\n")
nutrition_data <- read_excel("sampled_data_nutrition.xlsx")

# Display basic information about the dataset ------------------------------------
cat("Dataset dimensions:", dim(nutrition_data), "\n")
cat("Column names:\n")
print(names(nutrition_data))

# View first few rows to understand data structure -------------------------------
cat("First 6 rows of nutrition data:\n")
print(head(nutrition_data))

# Check data structure -----------------------------------------------------------
cat("Data structure:\n")
str(nutrition_data)

# ===============================================================================
# FOOD GROUP VARIABLES SELECTION
# ===============================================================================

# According to the documentation, variables 56-73 represent:
# Food groups Percentage Energy Contribution: 
# [56] FG1 Percentage Energy Contribution through [73] FG18 Percentage Energy Contribution
# These percentages are computed as ratio between each food group energy (variables 38-55) 
# and Total Energy (variable 9)

# Extract the food group percentage variables (columns 56-73) --------------------
# First, let's identify which columns correspond to variables 56-73
# We'll assume column names follow a pattern or we'll use position-based selection

# Method 1: If column names are descriptive
food_group_cols <- names(nutrition_data) %>% 
  grep("FG.*Percentage|Percentage.*FG|Food.*Group", ., value = TRUE, ignore.case = TRUE)

# Method 2: If no clear pattern, use position-based selection (columns 56-73)
if (length(food_group_cols) == 0) {
  food_group_cols <- names(nutrition_data)[56:73]
}

cat("Selected food group percentage columns:\n")
print(food_group_cols)

# Extract the food group percentage data ----------------------------------------
food_group_data <- nutrition_data %>% 
  select(all_of(food_group_cols))

# Display summary statistics for food group percentages -------------------------
cat("Summary statistics for food group percentages:\n")
print(summary(food_group_data))

# Check for missing values -------------------------------------------------------
cat("Missing values per food group:\n")
print(colSums(is.na(food_group_data)))

# ===============================================================================
# DATA PREPROCESSING FOR hmmTMB
# ===============================================================================

# First, let's examine the actual column names to identify ID and time variables
cat("All column names in dataset:\n")
print(names(nutrition_data))

# Look for potential ID columns (subject/participant identifiers)
id_patterns <- c("ID", "id", "subject", "participant", "person", "user", "respondent")
id_col <- NULL
for (pattern in id_patterns) {
  matches <- grep(pattern, names(nutrition_data), ignore.case = TRUE, value = TRUE)
  if (length(matches) > 0) {
    id_col <- matches[1]  # Take first match
    cat("Found ID column:", id_col, "\n")
    break
  }
}

# Look for potential time/visit columns
time_patterns <- c("time", "visit", "date", "day", "week", "month", "occasion", "measurement")
time_col <- NULL
for (pattern in time_patterns) {
  matches <- grep(pattern, names(nutrition_data), ignore.case = TRUE, value = TRUE)
  if (length(matches) > 0) {
    time_col <- matches[1]  # Take first match
    cat("Found time column:", time_col, "\n")
    break
  }
}

# If no ID column found, create a simple row identifier
if (is.null(id_col)) {
  nutrition_data$subject_id <- seq_len(nrow(nutrition_data))
  id_col <- "subject_id"
  cat("Created subject_id column as row identifier\n")
}

# If no time column found, create a simple visit number
if (is.null(time_col)) {
  nutrition_data$visit_number <- seq_len(nrow(nutrition_data))
  time_col <- "visit_number"
  cat("Created visit_number column as time indicator\n")
}

# Extract the food group percentage variables (columns 56-73) --------------------
# Method 1: If column names are descriptive
food_group_cols <- names(nutrition_data) %>% 
  grep("FG.*Percentage|Percentage.*FG|Food.*Group|Energy.*Contribution", ., value = TRUE, ignore.case = TRUE)

# Method 2: If no clear pattern, use position-based selection (columns 56-73)
if (length(food_group_cols) == 0) {
  if (ncol(nutrition_data) >= 73) {
    food_group_cols <- names(nutrition_data)[56:73]
  } else {
    # Use last 18 columns as food groups if dataset has fewer columns
    food_group_cols <- tail(names(nutrition_data), 18)
  }
}

cat("Selected food group percentage columns:\n")
print(food_group_cols)

# Create dataset with ID, time, and food group variables ------------------------
hmm_data <- nutrition_data %>% 
  select(all_of(c(id_col, time_col, food_group_cols)))

cat("HMM dataset dimensions:", dim(hmm_data), "\n")
cat("HMM dataset columns:\n")
print(names(hmm_data))

# Check for missing values in key variables -----------------------------------
cat("Missing values summary:\n")
cat("ID column missing:", sum(is.na(hmm_data[[id_col]])), "\n")
cat("Time column missing:", sum(is.na(hmm_data[[time_col]])), "\n")
cat("Food group columns missing:\n")
food_group_missing <- colSums(is.na(hmm_data[, food_group_cols]))
print(food_group_missing)

# Remove rows with missing values in food group variables -----------------------
# Keep ID and time columns for interpretation
complete_indices <- complete.cases(hmm_data[, food_group_cols])
complete_data <- hmm_data[complete_indices, ]

cat("Original observations:", nrow(hmm_data), "\n")
cat("Complete cases:", nrow(complete_data), "\n")
cat("Removed due to missing values:", nrow(hmm_data) - nrow(complete_data), "\n")

# Extract only numeric food group columns for scaling ---------------------------
food_group_numeric <- complete_data[, food_group_cols, drop = FALSE]

# Convert all columns to numeric (in case some are character/factor)
for (col in names(food_group_numeric)) {
  if (!is.numeric(food_group_numeric[[col]])) {
    food_group_numeric[[col]] <- as.numeric(as.character(food_group_numeric[[col]]))
  }
}

# Check if all data is numeric now
cat("Data types after conversion:\n")
print(sapply(food_group_numeric, class))

# Remove any remaining non-numeric rows
numeric_indices <- complete.cases(food_group_numeric)
food_group_numeric <- food_group_numeric[numeric_indices, ]

# Also filter the ID and time columns to match
complete_data <- complete_data[numeric_indices, ]

cat("Final dataset after numeric conversion:", nrow(complete_data), "rows\n")

# Standardize the food group percentages -----------------------------------------
# This helps with model convergence and interpretation
# We'll scale each food group to have mean=0 and sd=1
if (nrow(food_group_numeric) > 0 && ncol(food_group_numeric) > 0) {
  scaled_food_groups <- scale(food_group_numeric)
  
  # Convert back to data frame with proper column names
  scaled_food_groups <- as.data.frame(scaled_food_groups)
  names(scaled_food_groups) <- food_group_cols
  
  # Add back ID and time columns
  scaled_food_groups[[id_col]] <- complete_data[[id_col]]
  scaled_food_groups[[time_col]] <- complete_data[[time_col]]
  
  cat("Data after scaling - summary:\n")
  print(summary(scaled_food_groups[, food_group_cols]))
  
  cat("Scaling completed successfully!\n")
} else {
  cat("No valid numeric data for scaling\n")
  scaled_food_groups <- NULL
}

# ===============================================================================
# HIDDEN MARKOV MODEL SETUP
# ===============================================================================

# Define the number of hidden states -------------------------------------------
# For nutritional patterns, we might expect 2-4 states (e.g., low, medium, high intake)
# Let's start with 3 states as a reasonable assumption
n_states <- 3

# Create observation formula for hmmTMB ----------------------------------------
# We'll model the food group percentages as multivariate normal observations
# Each state will have different mean patterns for the 18 food groups

# Formula specification for hmmTMB:
# ~ 1 indicates intercept-only model (no covariates on observation parameters)
# The response variables are the food group percentages
obs_formula <- as.formula(paste("~ 1"))

# ===============================================================================
# MODEL FITTING - BASIC HMM
# ===============================================================================

# Check if we have valid data for modeling
if (is.null(scaled_food_groups) || nrow(scaled_food_groups) == 0) {
  cat("ERROR: No valid data available for HMM modeling\n")
  stop("Cannot proceed with HMM analysis - check data preprocessing")
}

# Extract only the food group columns for modeling (exclude ID and time)
model_data <- scaled_food_groups[, food_group_cols, drop = FALSE]

cat("Fitting basic Hidden Markov Model with", n_states, "states...\n")
cat("Using", nrow(model_data), "observations and", ncol(model_data), "food group variables\n")

# Fit the hmmTMB model ---------------------------------------------------------
# Parameters explanation:
# - par: initial parameter values (NULL for automatic initialization)
# - formula: observation model formula
# - data: observation data (food group percentages)
# - family: distribution family (gaussian for continuous percentages)
# - n_states: number of hidden states
# - optimizer: optimization method ("NLopt" is recommended)

# Note: This is a basic model. We may need to adjust based on data characteristics
basic_hmm <- tryCatch({
  hmmTMB::fit_hmm(
    par = NULL,                    # Let hmmTMB choose initial parameters
    formula = obs_formula,         # Observation model
    data = model_data,             # Scaled food group percentages only
    family = gaussian(),           # Gaussian distribution for continuous data
    n_states = n_states,           # Number of hidden states
    optimizer = "NLopt"            # Use NLopt optimizer
  )
}, error = function(e) {
  cat("Error in basic model fitting:", e$message, "\n")
  cat("This might be due to:\n")
  cat("- Insufficient data for the number of states\n")
  cat("- Poor initial parameter values\n")
  cat("- Numerical issues with the optimization\n")
  cat("- Data structure problems\n")
  return(NULL)
})

# ===============================================================================
# MODEL FITTING - WITH COVARIATES ON TRANSITION PROBABILITIES
# ===============================================================================

# If we have a time column, we could model its effect on transition probabilities
# For time series analysis, we need to ensure data is properly ordered
time_covariate <- scaled_food_groups[[time_col]]

# Check if time variable is suitable for modeling
if (is.numeric(time_covariate) && length(unique(time_covariate)) > 1) {
  cat("Fitting HMM with time-varying transition probabilities...\n")
  cat("Time variable range:", min(time_covariate), "to", max(time_covariate), "\n")
  
  time_hmm <- tryCatch({
    hmmTMB::fit_hmm(
      par = NULL,
      formula = obs_formula,
      data = model_data,
      family = gaussian(),
      n_states = n_states,
      optimizer = "NLopt",
      # Add covariate effect on transition probabilities
      transition_formula = ~ time_covariate
    )
  }, error = function(e) {
    cat("Error in time-varying model:", e$message, "\n")
    return(NULL)
  })
} else {
  cat("Time variable not suitable for modeling - skipping time-varying transitions\n")
  time_hmm <- NULL
}

# ===============================================================================
# MODEL RESULTS AND INTERPRETATION
# ===============================================================================

# Function to display model results ----------------------------------------------
display_model_results <- function(model, model_name) {
  if (is.null(model)) {
    cat("Model", model_name, "failed to fit\n\n")
    return()
  }
  
  cat("=== RESULTS FOR", model_name, "===\n")
  
  # Model convergence information
  cat("Convergence code:", model$optim_info$convergence, "\n")
  cat("Number of function evaluations:", model$optim_info$iterations, "\n")
  cat("Log-likelihood:", model$loglik, "\n")
  cat("AIC:", AIC(model), "\n")
  
  # State-dependent means (for each food group in each state)
  cat("\nState-dependent means (food group percentages):\n")
  print(model$par$obs)
  
  # Transition probability matrix
  cat("\nTransition probability matrix:\n")
  print(model$par$tpm)
  
  # Stationary distribution
  cat("\nStationary distribution:\n")
  print(model$par$delta)
  
  cat("\n")
}

# Display results for both models -----------------------------------------------
display_model_results(basic_hmm, "BASIC HMM")
display_model_results(time_hmm, "TIME-VARYING HMM")

# ===============================================================================
# STATE DECODING AND VISUALIZATION
# ===============================================================================

# Function to decode states and create visualizations ---------------------------
analyze_and_visualize <- function(model, data, model_name, id_col_name, time_col_name, complete_data_full) {
  if (is.null(model)) {
    cat("Cannot analyze", model_name, "- model failed\n\n")
    return()
  }
  
  cat("=== STATE ANALYSIS FOR", model_name, "===\n")
  
  # Decode the most likely state sequence using Viterbi algorithm
  viterbi_states <- hmmTMB::viterbi(model, data)
  
  # Get state probabilities (forward-backward algorithm)
  state_probs <- hmmTMB::state_probs(model, data)
  
  # Add decoded states to analysis data
  analysis_data <- data.frame(
    decoded_state = viterbi_states,
    time_index = seq_len(nrow(data)),
    id = complete_data_full[[id_col_name]],
    time = complete_data_full[[time_col_name]]
  )
  
  # Add food group data for interpretation
  for (col in food_group_cols) {
    analysis_data[[col]] <- data[[col]]
  }
  
  # Summary of state assignments
  cat("State assignment summary:\n")
  print(table(analysis_data$decoded_state))
  
  # Calculate mean food group percentages by state (in original scale)
  # Convert back from standardized to original scale
  original_scale_data <- data.frame()
  for (col in food_group_cols) {
    # Reverse the scaling: original = standardized * sd + mean
    original_values <- data[[col]] * attr(food_group_numeric, "scaled:scale")[col] + 
                       attr(food_group_numeric, "scaled:center")[col]
    original_scale_data[[col]] <- original_values
  }
  original_scale_data$decoded_state <- viterbi_states
  
  state_means <- original_scale_data %>% 
    group_by(decoded_state) %>% 
    summarise(across(all_of(food_group_cols), mean, .names = "mean_{.col}"))
  
  cat("\nMean food group percentages by state (original scale):\n")
  print(round(state_means[, -1], 2))  # Exclude decoded_state column
  
  # Create visualization of state assignments over time
  p1 <- ggplot(analysis_data, aes(x = time_index, y = decoded_state, color = factor(decoded_state))) +
    geom_point(alpha = 0.6) +
    geom_line(alpha = 0.3) +
    labs(title = paste("Decoded States Over Time -", model_name),
         x = "Observation Index",
         y = "Hidden State",
         color = "State") +
    theme_minimal()
  
  # Create heatmap of state means
  state_means_long <- state_means %>% 
    pivot_longer(cols = starts_with("mean_"), 
                 names_to = "food_group", 
                 values_to = "mean_percentage") %>%
    mutate(food_group = gsub("mean_", "", food_group))
  
  p2 <- ggplot(state_means_long, aes(x = food_group, y = decoded_state, fill = mean_percentage)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(title = paste("Mean Food Group Percentages by State -", model_name),
         x = "Food Group",
         y = "Hidden State",
         fill = "Mean Percentage") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Display plots
  print(p1)
  print(p2)
  
  # Return analysis results for further use
  return(list(
    analysis_data = analysis_data,
    state_means = state_means,
    viterbi_states = viterbi_states,
    state_probs = state_probs
  ))
}

# Analyze and visualize both models ---------------------------------------------
basic_results <- analyze_and_visualize(basic_hmm, model_data, "BASIC HMM", id_col, time_col, complete_data)
time_results <- analyze_and_visualize(time_hmm, model_data, "TIME-VARYING HMM", id_col, time_col, complete_data)

# ===============================================================================
# MODEL COMPARISON AND SELECTION
# ===============================================================================

# Compare models using information criteria -------------------------------------
cat("=== MODEL COMPARISON ===\n")

if (!is.null(basic_hmm) && !is.null(time_hmm)) {
  models_list <- list(Basic = basic_hmm, TimeVarying = time_hmm)
  
  # Calculate AIC and BIC for each model
  comparison_table <- data.frame(
    Model = names(models_list),
    AIC = sapply(models_list, AIC),
    BIC = sapply(models_list, BIC),
    LogLik = sapply(models_list, function(x) x$loglik),
    Parameters = sapply(models_list, function(x) length(x$par$obs) + 
                  length(x$par$tpm) + length(x$par$delta) - 1)
  )
  
  print(comparison_table)
  
  # Determine which model is better (lower AIC/BIC is preferred)
  best_aic <- comparison_table$Model[which.min(comparison_table$AIC)]
  best_bic <- comparison_table$Model[which.min(comparison_table$BIC)]
  
  cat("Best model by AIC:", best_aic, "\n")
  cat("Best model by BIC:", best_bic, "\n")
} else {
  cat("Cannot compare models - one or both failed to fit\n")
}

# ===============================================================================
# SAVING RESULTS
# ===============================================================================

# Save the best model (if available) -------------------------------------------
best_model <- if (!is.null(basic_hmm) && !is.null(time_hmm)) {
  if (AIC(basic_hmm) < AIC(time_hmm)) basic_hmm else time_hmm
} else if (!is.null(basic_hmm)) {
  basic_hmm
} else if (!is.null(time_hmm)) {
  time_hmm
} else {
  NULL
}

if (!is.null(best_model)) {
  saveRDS(best_model, "best_hmm_model.rds")
  cat("Best model saved as 'best_hmm_model.rds'\n")
}

# Save analysis results --------------------------------------------------------
results_summary <- list(
  food_group_columns = food_group_cols,
  n_states = n_states,
  n_observations = nrow(scaled_food_groups),
  basic_model_fitted = !is.null(basic_hmm),
  time_model_fitted = !is.null(time_hmm),
  best_model_type = if (!is.null(best_model)) {
    if (identical(best_model, basic_hmm)) "Basic" else "Time-Varying"
  } else {
    "None"
  }
)

saveRDS(results_summary, "hmm_analysis_summary.rds")
cat("Analysis summary saved as 'hmm_analysis_summary.rds'\n")

# ===============================================================================
# END OF ANALYSIS
# ===============================================================================

cat("Hidden Markov Model analysis completed!\n")
cat("Check the generated plots and saved results for detailed interpretation.\n")
