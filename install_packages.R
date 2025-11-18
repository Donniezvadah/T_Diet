# Install all required packages for longitudinal data analysis
# This script reads the requirements.txt file and installs all packages

# Set working directory
setwd("/Users/donnie/Desktop/T_Diet")

# Read requirements file
requirements <- readLines("requirements.txt")
requirements <- requirements[!grepl("^#", requirements)]  # Remove comments
requirements <- requirements[requirements != ""]           # Remove empty lines

cat("Installing", length(requirements), "packages...\n")

# Function to install packages safely
install_package_safely <- function(pkg) {
  cat("Installing package:", pkg, "\n")
  tryCatch({
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cran.rstudio.com/", dependencies = TRUE)
      if (require(pkg, character.only = TRUE, quietly = TRUE)) {
        cat("✓ Successfully installed and loaded", pkg, "\n")
        return(TRUE)
      } else {
        cat("✗ Failed to load", pkg, "after installation\n")
        return(FALSE)
      }
    } else {
      cat("✓ Package", pkg, "already installed\n")
      return(TRUE)
    }
  }, error = function(e) {
    cat("✗ Error installing", pkg, ":", e$message, "\n")
    return(FALSE)
  })
}

# Install all packages
installation_results <- sapply(requirements, install_package_safely)

# Summary
successful_installs <- sum(installation_results)
failed_installs <- length(installation_results) - successful_installs

cat("\n" + "="*50 + "\n")
cat("INSTALLATION SUMMARY\n")
cat("="*50 + "\n")
cat("Successfully installed:", successful_installs, "packages\n")
cat("Failed to install:", failed_installs, "packages\n")

if (failed_installs > 0) {
  failed_packages <- names(installation_results)[!installation_results]
  cat("Failed packages:", paste(failed_packages, collapse = ", "), "\n")
}

cat("\nAll essential packages are now ready for longitudinal data analysis!\n")
