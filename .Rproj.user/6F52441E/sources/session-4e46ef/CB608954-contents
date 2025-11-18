# Comprehensive Longitudinal Data Analysis Script
# Nutrition Dataset Exploratory Data Analysis

# ===== LIBRARY IMPORTS =====
# Data manipulation and wrangling
library(tidyverse)        # Core data science toolkit (dplyr, ggplot2, etc.)
library(data.table)       # Fast data manipulation for large datasets
library(dplyr)            # Data manipulation grammar
library(tidyr)            # Data tidying
library(lubridate)        # Date and time manipulation
library(stringr)          # String manipulation

# Data import and export
library(readxl)           # Excel file reading
# library(writexl)          # Excel file writing
# library(openxlsx)         # Advanced Excel operations

# Statistical analysis
library(stats)            # Basic statistical functions
library(car)              # Regression diagnostics
library(MASS)             # Modern Applied Statistics with S
library(nlme)             # Nonlinear mixed-effects models
library(lme4)             # Linear mixed-effects models
library(lmerTest)         # p-values for lmer models

# Visualization
library(ggplot2)          # Grammar of graphics
library(ggplot2.extensions) # Additional ggplot2 extensions
library(ggpubr)           # Publication ready plots
library(patchwork)        # Combining ggplot objects
library(cowplot)          # Publication quality themes
library(scales)           # Scale functions for visualization
library(RColorBrewer)     # Color palettes
library(viridis)          # Colorblind-friendly palettes

# Specialized longitudinal analysis
library(ggfortify)        # Autoplotting for statistical models
library(broom)            # Tidying model outputs
library(broom.mixed)      # Tidying mixed model outputs
library(emmeans)          # Estimated marginal means
library(effects)          # Effect displays for statistical models

# Time series analysis
library(zoo)              # Time series utilities
library(xts)              # Extensible time series
library(tsibble)          # Tidy time series data

# Summary statistics and tables
library(knitr)            # Dynamic report generation
library(kableExtra)       # Enhanced table formatting
library(gt)               # Grammar of tables
library(flextable)        # Flexible table creation

# Missing data analysis
library(naniar)           # Exploring missing values
library(VIM)              # Visualization and imputation of missing values
library(mice)             # Multivariate imputation by chained equations

# Correlation and association analysis
library(corrplot)         # Visualization of correlation matrices
library(Hmisc)            # Miscellaneous functions for data analysis
library(psych)            # Psychological research methods

# Outlier detection
library(outliers)         # Tests for outliers
library(DescTools)        # Descriptive statistics and tools

# Data quality checks
library(validate)         # Data validation
library(assertthat)       # Readable assertions for R

# Performance optimization
library(microbenchmark)   # Precise timing
library(profvis)          # Visual profiling for R code

# ===== SET UP =====
# Set working directory and options
setwd("/Users/donnie/Desktop/T_Diet")

# Global options for better output
options(stringsAsFactors = FALSE)
options(scipen = 999)  # Disable scientific notation
options(max.print = 1000)

# Set random seed for reproducibility
set.seed(123)

# Create output directory for plots and results
if (!dir.exists("analysis_results")) {
  dir.create("analysis_results")
  dir.create("analysis_results/plots")
  dir.create("analysis_results/tables")
  dir.create("analysis_results/models")
}

# ===== DATA IMPORT =====
cat("Importing nutrition dataset...\n")

# Read the Excel file
tryCatch({
  nutrition_data <- read_excel("sampled_data_nutrition.xlsx")
  cat("Dataset imported successfully!\n")
}, error = function(e) {
  cat("Error importing dataset:", e$message, "\n")
  stop("Dataset import failed")
})

# ===== BASIC DATA INSPECTION =====
cat("\n===== BASIC DATA INSPECTION =====\n")

# Dataset dimensions
cat("Dataset dimensions:", dim(nutrition_data), "\n")

# Column names
cat("Column names:\n")
print(colnames(nutrition_data))

# First few rows
cat("\nFirst 6 rows of data:\n")
print(head(nutrition_data, 6))

# Last few rows
cat("\nLast 6 rows of data:\n")
print(tail(nutrition_data, 6))

# Data structure
cat("\nData structure:\n")
str(nutrition_data)

# ===== DATA TYPE ANALYSIS =====
cat("\n===== DATA TYPE ANALYSIS =====\n")

# Identify data types
data_types <- sapply(nutrition_data, class)
type_summary <- table(data_types)
print(type_summary)

# Separate columns by type
numeric_cols <- names(nutrition_data)[sapply(nutrition_data, is.numeric)]
character_cols <- names(nutrition_data)[sapply(nutrition_data, is.character)]
factor_cols <- names(nutrition_data)[sapply(nutrition_data, is.factor)]
date_cols <- names(nutrition_data)[sapply(nutrition_data, inherits, "Date") || 
                                   sapply(nutrition_data, inherits, "POSIXt")]

cat("Numeric columns:", length(numeric_cols), "\n")
cat("Character columns:", length(character_cols), "\n")
cat("Factor columns:", length(factor_cols), "\n")
cat("Date columns:", length(date_cols), "\n")

# ===== MISSING DATA ANALYSIS =====
cat("\n===== MISSING DATA ANALYSIS =====\n")

# Overall missing data summary
missing_summary <- naniar::miss_var_summary(nutrition_data)
print(missing_summary)

# Missing data patterns
missing_patterns <- naniar::miss_case_summary(nutrition_data)
print(missing_patterns)

# Visualize missing data
p1 <- naniar::vis_miss(nutrition_data) + 
  ggtitle("Missing Data Pattern") +
  theme_minimal()

ggsave("analysis_results/plots/missing_data_pattern.png", p1, width = 10, height = 6)

# Missing data heatmap if there are missing values
if (any(is.na(nutrition_data))) {
  p2 <- VIM::aggr(nutrition_data, col = c('navyblue', 'red'), 
                  numbers = TRUE, sortVars = TRUE)
  dev.copy(png, "analysis_results/plots/missing_data_heatmap.png", width = 800, height = 600)
  dev.off()
}

# ===== DESCRIPTIVE STATISTICS =====
cat("\n===== DESCRIPTIVE STATISTICS =====\n")

# Basic descriptive statistics for numeric variables
if (length(numeric_cols) > 0) {
  basic_stats <- nutrition_data %>%
    select(all_of(numeric_cols)) %>%
    summarise(
      across(everything(), list(
        n = ~sum(!is.na(.)),
        mean = ~mean(., na.rm = TRUE),
        median = ~median(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        min = ~min(., na.rm = TRUE),
        max = ~max(., na.rm = TRUE),
        q25 = ~quantile(., 0.25, na.rm = TRUE),
        q75 = ~quantile(., 0.75, na.rm = TRUE)
      ))
    )
  
  print(basic_stats)
  
  # Save descriptive statistics table
  write.csv(basic_stats, "analysis_results/tables/descriptive_statistics.csv", row.names = FALSE)
}

# Frequency tables for categorical variables
if (length(character_cols) > 0 || length(factor_cols) > 0) {
  cat_cols <- c(character_cols, factor_cols)
  
  for (col in cat_cols) {
    cat("\nFrequency table for", col, ":\n")
    freq_table <- table(nutrition_data[[col]], useNA = "ifany")
    print(freq_table)
    
    # Save frequency tables
    write.csv(as.data.frame(freq_table), 
              paste0("analysis_results/tables/freq_table_", col, ".csv"), 
              row.names = FALSE)
  }
}

# ===== LONGITUDINAL DATA ANALYSIS =====
cat("\n===== LONGITUDINAL DATA ANALYSIS =====\n")

# Identify potential longitudinal structure
# Look for columns that might represent time or subject ID
potential_time_cols <- grep("time|date|week|month|year|visit|assessment", 
                           colnames(nutrition_data), ignore.case = TRUE, value = TRUE)
potential_id_cols <- grep("id|subject|participant|patient", 
                         colnames(nutrition_data), ignore.case = TRUE, value = TRUE)

cat("Potential time-related columns:", potential_time_cols, "\n")
cat("Potential ID columns:", potential_id_cols, "\n")

# If we have potential longitudinal structure, analyze it
if (length(potential_id_cols) > 0 && length(potential_time_cols) > 0) {
  id_col <- potential_id_cols[1]
  time_col <- potential_time_cols[1]
  
  cat("Analyzing longitudinal structure with ID:", id_col, "and Time:", time_col, "\n")
  
  # Number of unique subjects
  n_subjects <- length(unique(nutrition_data[[id_col]]))
  cat("Number of unique subjects:", n_subjects, "\n")
  
  # Time points per subject
  time_points <- nutrition_data %>%
    group_by(!!sym(id_col)) %>%
    summarise(n_timepoints = n(), .groups = "drop")
  
  cat("Summary of time points per subject:\n")
  print(summary(time_points$n_timepoints))
  
  # Distribution of time points
  p3 <- ggplot(time_points, aes(x = n_timepoints)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
    ggtitle("Distribution of Time Points per Subject") +
    xlab("Number of Time Points") +
    ylab("Number of Subjects") +
    theme_minimal()
  
  ggsave("analysis_results/plots/time_points_distribution.png", p3, width = 8, height = 6)
  
  # Time range
  if (is.numeric(nutrition_data[[time_col]])) {
    time_range <- range(nutrition_data[[time_col]], na.rm = TRUE)
    cat("Time range:", time_range[1], "to", time_range[2], "\n")
  }
}

# ===== DISTRIBUTION ANALYSIS =====
cat("\n===== DISTRIBUTION ANALYSIS =====\n")

# Distribution plots for numeric variables
if (length(numeric_cols) > 0) {
  # Create distribution plots
  n_vars <- length(numeric_cols)
  n_cols_plot <- min(3, n_vars)
  n_rows_plot <- ceiling(n_vars / n_cols_plot)
  
  # Histograms
  p_hist <- nutrition_data %>%
    select(all_of(numeric_cols)) %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 30, fill = "steelblue", color = "black", alpha = 0.7) +
    facet_wrap(~variable, scales = "free", ncol = n_cols_plot) +
    ggtitle("Distribution of Numeric Variables") +
    theme_minimal()
  
  ggsave("analysis_results/plots/distributions_histograms.png", p_hist, 
         width = n_cols_plot * 4, height = n_rows_plot * 3)
  
  # Box plots
  p_box <- nutrition_data %>%
    select(all_of(numeric_cols)) %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
    ggplot(aes(x = variable, y = value)) +
    geom_boxplot(fill = "lightblue", alpha = 0.7) +
    coord_flip() +
    ggtitle("Box Plots of Numeric Variables") +
    theme_minimal()
  
  ggsave("analysis_results/plots/distributions_boxplots.png", p_box, 
         width = 8, height = n_vars * 0.5)
  
  # Density plots
  p_density <- nutrition_data %>%
    select(all_of(numeric_cols)) %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
    ggplot(aes(x = value, fill = variable)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~variable, scales = "free", ncol = n_cols_plot) +
    ggtitle("Density Plots of Numeric Variables") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave("analysis_results/plots/distributions_density.png", p_density, 
         width = n_cols_plot * 4, height = n_rows_plot * 3)
}

# ===== OUTLIER DETECTION =====
cat("\n===== OUTLIER DETECTION =====\n")

if (length(numeric_cols) > 0) {
  outlier_summary <- data.frame(
    variable = numeric_cols,
    n_outliers = NA_integer_,
    outlier_percentage = NA_numeric_
  )
  
  for (i in seq_along(numeric_cols)) {
    var <- numeric_cols[i]
    values <- nutrition_data[[var]]
    
    # Remove NA values
    values_clean <- values[!is.na(values)]
    
    if (length(values_clean) > 0) {
      # IQR method for outlier detection
      Q1 <- quantile(values_clean, 0.25)
      Q3 <- quantile(values_clean, 0.75)
      IQR <- Q3 - Q1
      lower_bound <- Q1 - 1.5 * IQR
      upper_bound <- Q3 + 1.5 * IQR
      
      outliers <- values_clean < lower_bound | values_clean > upper_bound
      n_outliers <- sum(outliers)
      
      outlier_summary$n_outliers[i] <- n_outliers
      outlier_summary$outlier_percentage[i] <- (n_outliers / length(values_clean)) * 100
      
      cat("Variable", var, ":", n_outliers, "outliers (", 
          round(outlier_summary$outlier_percentage[i], 2), "%)\n")
    }
  }
  
  # Save outlier summary
  write.csv(outlier_summary, "analysis_results/tables/outlier_summary.csv", row.names = FALSE)
  
  # Visualize outliers
  p_outliers <- nutrition_data %>%
    select(all_of(numeric_cols)) %>%
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
    ggplot(aes(x = variable, y = value)) +
    geom_boxplot(fill = "lightcoral", alpha = 0.7, outlier.color = "red", outlier.shape = 16) +
    coord_flip() +
    ggtitle("Outlier Detection - Box Plots") +
    theme_minimal()
  
  ggsave("analysis_results/plots/outlier_detection.png", p_outliers, 
         width = 8, height = n_vars * 0.5)
}

# ===== CORRELATION ANALYSIS =====
cat("\n===== CORRELATION ANALYSIS =====\n")

if (length(numeric_cols) > 1) {
  # Calculate correlation matrix
  cor_matrix <- cor(nutrition_data[numeric_cols], use = "pairwise.complete.obs")
  
  # Save correlation matrix
  write.csv(cor_matrix, "analysis_results/tables/correlation_matrix.csv")
  
  # Correlation plot
  p_cor <- corrplot(cor_matrix, method = "color", type = "upper", 
                    order = "hclust", tl.cex = 0.8, tl.col = "black",
                    title = "Correlation Matrix", mar = c(0,0,1,0))
  
  dev.copy(png, "analysis_results/plots/correlation_matrix.png", width = 800, height = 800)
  dev.off()
  
  # Find highly correlated pairs
  high_cor_pairs <- which(abs(cor_matrix) > 0.7 & cor_matrix != 1, arr.ind = TRUE)
  if (nrow(high_cor_pairs) > 0) {
    cat("Highly correlated variable pairs (|r| > 0.7):\n")
    for (i in 1:nrow(high_cor_pairs)) {
      row_idx <- high_cor_pairs[i, 1]
      col_idx <- high_cor_pairs[i, 2]
      if (row_idx < col_idx) {  # Avoid duplicates
        cat(numeric_cols[row_idx], "-", numeric_cols[col_idx], ": ", 
            round(cor_matrix[row_idx, col_idx], 3), "\n")
      }
    }
  }
}

# ===== TIME SERIES ANALYSIS (if applicable) =====
if (length(potential_time_cols) > 0 && length(numeric_cols) > 0) {
  cat("\n===== TIME SERIES ANALYSIS =====\n")
  
  time_col <- potential_time_cols[1]
  
  # If time column is numeric, create time series plots
  if (is.numeric(nutrition_data[[time_col]])) {
    # Sort by time
    data_sorted <- nutrition_data[order(nutrition_data[[time_col]]), ]
    
    # Plot key variables over time
    key_vars <- numeric_cols[1:min(4, length(numeric_cols))]  # Plot up to 4 variables
    
    for (var in key_vars) {
      p_ts <- ggplot(data_sorted, aes(x = !!sym(time_col), y = !!sym(var))) +
        geom_line(alpha = 0.7) +
        geom_smooth(method = "loess", se = TRUE, color = "red") +
        ggtitle(paste("Time Series:", var)) +
        xlab("Time") +
        ylab(var) +
        theme_minimal()
      
      ggsave(paste0("analysis_results/plots/timeseries_", var, ".png"), p_ts, 
             width = 10, height = 6)
    }
  }
}

# ===== SUMMARY REPORT =====
cat("\n===== SUMMARY REPORT =====\n")

# Create a comprehensive summary report
summary_report <- list(
  dataset_info = list(
    dimensions = dim(nutrition_data),
    columns = colnames(nutrition_data),
    numeric_columns = numeric_cols,
    categorical_columns = c(character_cols, factor_cols),
    date_columns = date_cols
  ),
  missing_data = missing_summary,
  longitudinal_structure = list(
    id_columns = potential_id_cols,
    time_columns = potential_time_cols,
    n_subjects = if (length(potential_id_cols) > 0) length(unique(nutrition_data[[potential_id_cols[1]]])) else NA
  ),
  outliers = outlier_summary
)

# Save summary report
saveRDS(summary_report, "analysis_results/summary_report.rds")

# Print final summary
cat("\n" + "="*50 + "\n")
cat("ANALYSIS COMPLETE\n")
cat("="*50 + "\n")
cat("Dataset dimensions:", dim(nutrition_data)[1], "rows,", dim(nutrition_data)[2], "columns\n")
cat("Numeric variables:", length(numeric_cols), "\n")
cat("Categorical variables:", length(c(character_cols, factor_cols)), "\n")
if (any(is.na(nutrition_data))) {
  cat("Missing values detected - see missing_data_pattern.png\n")
} else {
  cat("No missing values detected\n")
}
if (length(potential_id_cols) > 0 && length(potential_time_cols) > 0) {
  cat("Longitudinal structure detected with", summary_report$longitudinal_structure$n_subjects, "subjects\n")
} else {
  cat("No clear longitudinal structure detected\n")
}
cat("All plots saved in analysis_results/plots/\n")
cat("All tables saved in analysis_results/tables/\n")
cat("="*50 + "\n")

cat("\nNext steps for longitudinal analysis:\n")
cat("1. Examine the longitudinal structure plots\n")
cat("2. Review missing data patterns\n")
cat("3. Check correlation matrix for multicollinearity\n")
cat("4. Consider mixed-effects models for longitudinal analysis\n")
cat("5. Perform time series analysis if appropriate\n")