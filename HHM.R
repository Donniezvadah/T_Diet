
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