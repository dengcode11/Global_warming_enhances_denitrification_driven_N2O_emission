# ----------------------------------------------------------------
# Analysis of denitrification and denitrification-driven N₂O Emission Rate Increase ΔR% -
# 5-Model Ensemble Analysis (using Q10 values)
# Includes boxplot analysis for 2090-2100 period
# ----------------------------------------------------------------

# Load required libraries
library(ncdf4)      # For reading NetCDF files
library(raster)     # For handling raster (gridded) data
library(tidyverse)  # For data manipulation and visualization
library(lubridate)  # For date-time operations
library(patchwork)  # For combining multiple plots
library(rstatix)    # For statistical tests
library(ggpubr)     # For publication-ready plots
library(ggsignif)   # For adding statistical significance annotations to plots

# ----------------------------------------------------------------
# 1. Define Key Parameters
# ----------------------------------------------------------------

# Base directory path containing climate model data files
# Note: Update this path according to your actual file location
base_path <- "/Users/mac/Downloads/SSP370/"

# Q10 represents the temperature sensitivity of biological processes
# Q10 values for denitrification rate specifically
Q10_values <- list(
  Plastisphere = 2.26,    # Q10 for denitrification in plastisphere communities
  Epiphyton = 2.56        # Q10 for denitrification in epiphyton communities
)

# Q10 values for denitrification-driven N2O emissions specifically
# Note: This redefines the previous list - in practice, choose one set of Q10 values
Q10_values <- list(
  Plastisphere = 2.54,    # Q10 for plastisphere-driven N2O emissions
  Epiphyton = 3.78        # Q10 for epiphyton-driven N2O emissions
)

# Baseline period for calculating temperature anomalies (ΔT)
baseline_period <- c(2015, 2025)  # Years included in baseline

# Period for boxplot analysis (end-of-century analysis)
boxplot_period <- c(2090, 2100)  # Years included in statistical analysis

# ----------------------------------------------------------------
# 2. Define File Paths
# ----------------------------------------------------------------

# MPI-ESM1-2-HR model uses segmented files (different time periods)
mpi_periods <- c(
  "201501-201912", "202001-202412", "202501-202912", "203001-203412",
  "203501-203912", "204001-204412", "204501-204912", "205001-205412",
  "205501-205912", "206001-206412", "206501-206912", "207001-207412",
  "207501-207912", "208001-208412", "208501-208912", "209001-209412",
  "209501-209912", "210001-210012"
)

# Create file paths for MPI-ESM1-2-HR model segments
mpi_files <- paste0(base_path, 
                    "tas_Amon_MPI-ESM1-2-HR_ssp370_r1i1p1f1_gn_", 
                    mpi_periods, ".nc")

# Check existence of MPI-ESM1-2-HR files
mpi_existing <- file.exists(mpi_files)
cat("MPI-ESM1-2-HR file existence check:\n")
for(i in 1:length(mpi_files)) {
  cat(paste0("  ", basename(mpi_files[i]), ": ", 
             ifelse(mpi_existing[i], "Exists", "Missing"), "\n"))
}

# Keep only existing MPI files
mpi_files_exist <- mpi_files[mpi_existing]

# Other climate models (single file each covering entire period)
other_models <- list(
  "ACCESS-ESM1-5" = paste0(base_path, "tas_Amon_ACCESS-ESM1-5_ssp370_r1i1p1f1_gn_201501-210012.nc"),
  "GFDL-ESM4" = paste0(base_path, "tas_Amon_GFDL-ESM4_ssp370_r1i1p1f1_gr1_201501-210012.nc"),
  "IPSL-CM6A-LR" = paste0(base_path, "tas_Amon_IPSL-CM6A-LR_ssp370_r1i1p1f1_gr_201501-210012.nc"),
  "MRI-ESM2-0" = paste0(base_path, "tas_Amon_MRI-ESM2-0_ssp370_r1i1p1f1_gn_201501-210012.nc")
)

# Check existence of other model files
other_existing <- sapply(other_models, file.exists)
cat("\nOther model file existence check:\n")
for(model in names(other_models)) {
  cat(paste0(model, ": ", ifelse(other_existing[model], "Exists", "Missing"), "\n"))
}

# Combine all model paths (only existing files)
all_models <- list()

# Add MPI-ESM1-2-HR if files exist
if(length(mpi_files_exist) > 0) {
  all_models[["MPI-ESM1-2-HR"]] <- mpi_files_exist
}

# Add other existing models
for(model in names(other_models)) {
  if(other_existing[model]) {
    all_models[[model]] <- other_models[[model]]
  }
}

# Check if we have at least 2 models for ensemble analysis
if(length(all_models) < 2) {
  stop(paste("At least 2 models are required for analysis, currently only", 
             length(all_models), "models available"))
}

cat(paste0("\nTotal models to process: ", length(all_models), "\n"))
cat(paste0("  ", paste(names(all_models), collapse = ", "), "\n"))

# ----------------------------------------------------------------
# 3. Define Helper Functions
# ----------------------------------------------------------------

# Simplified function to read single NetCDF file using raster package
read_single_nc_simple <- function(file_path) {
  tryCatch({
    cat(paste0("    Reading: ", basename(file_path), "\n"))
    
    # Read NetCDF file as raster brick (multi-layer raster)
    tas_brick <- brick(file_path, varname = "tas")
    
    # Get time information from raster
    time_days <- getZ(tas_brick)
    
    # If no time information in file, extract from filename
    if(is.null(time_days)) {
      # Extract year range from filename
      file_name <- basename(file_path)
      year_match <- regmatches(file_name, regexpr("\\d{6}-\\d{6}", file_name))
      if(length(year_match) > 0) {
        start_year <- as.numeric(substr(year_match[1], 1, 4))
        start_month <- as.numeric(substr(year_match[1], 5, 6))
        end_year <- as.numeric(substr(year_match[1], 8, 11))
        end_month <- as.numeric(substr(year_match[1], 12, 13))
        
        # Create monthly date sequence
        total_months <- nlayers(tas_brick)
        dates <- seq(from = as.Date(paste(start_year, start_month, "15", sep = "-")),
                     by = "month", length.out = total_months)
      } else {
        # Default to CMIP6 time convention
        dates <- as.Date(time_days, origin = as.Date("1850-01-01"))
      }
    } else {
      # Use standard CMIP6 time convention
      dates <- as.Date(time_days, origin = as.Date("1850-01-01"))
    }
    
    # Calculate global mean temperature (simple spatial average)
    global_tas <- cellStats(tas_brick, mean, na.rm = TRUE)
    
    # Create monthly data table
    monthly_data <- tibble(
      date = dates,
      year = year(dates),
      month = month(dates),
      tas_k = global_tas,
      tas_c = tas_k - 273.15  # Convert Kelvin to Celsius
    )
    
    return(monthly_data)
    
  }, error = function(e) {
    warning(paste("Failed to read file:", basename(file_path), "-", e$message))
    return(NULL)
  })
}

# Simplified function to process MPI-ESM1-2-HR multiple files
process_mpi_multifile_simple <- function(file_paths) {
  cat("    Processing MPI-ESM1-2-HR data (", length(file_paths), "files)...\n")
  
  all_monthly <- list()
  for(file_path in file_paths) {
    monthly_data <- read_single_nc_simple(file_path)
    if(!is.null(monthly_data)) {
      all_monthly <- c(all_monthly, list(monthly_data))
    }
  }
  
  if(length(all_monthly) == 0) {
    stop("No valid MPI-ESM1-2-HR data found")
  }
  
  # Combine all monthly data
  mpi_data <- bind_rows(all_monthly) %>%
    arrange(date) %>%
    distinct(date, .keep_all = TRUE)  # Ensure no duplicate dates
  
  cat(paste0("    MPI-ESM1-2-HR: ", nrow(mpi_data), " monthly data points\n"))
  cat(paste0("    Time range: ", min(mpi_data$date), " to ", max(mpi_data$date), "\n"))
  
  return(mpi_data)
}

# Simplified function to process single model file
process_single_model_simple <- function(model_name, file_path) {
  cat("  Processing", model_name, "...\n")
  
  # Handle MPI-ESM1-2-HR separately due to multiple files
  if(model_name == "MPI-ESM1-2-HR" && length(file_path) > 1) {
    monthly_data <- process_mpi_multifile_simple(file_path)
  } else {
    # Other models have single files
    monthly_data <- read_single_nc_simple(file_path)
  }
  
  if(is.null(monthly_data)) {
    warning(paste("Unable to process model:", model_name))
    return(NULL)
  }
  
  # Calculate annual mean temperature from monthly data
  annual_data <- monthly_data %>%
    group_by(year) %>%
    summarise(
      annual_tas = mean(tas_c, na.rm = TRUE),
      n_months = n(),
      .groups = "drop"
    ) %>%
    filter(n_months >= 10)  # Require at least 10 months of data per year
  
  # Check data completeness
  if(nrow(annual_data) == 0) {
    warning(paste("Model", model_name, "has no valid annual data"))
    return(NULL)
  }
  
  # Add model name identifier
  annual_data$model <- model_name
  
  cat(paste0("    ", model_name, ": ", nrow(annual_data), " annual data points\n"))
  cat(paste0("    Time range: ", min(annual_data$year), "-", max(annual_data$year), "\n"))
  
  return(annual_data)
}

# Function to calculate ΔR% for a single model using Q10 formula
calculate_deltaR_for_model <- function(model_data, Q10_val, Q10_name, baseline_period) {
  if(is.null(model_data)) return(NULL)
  
  # Calculate baseline period mean temperature
  baseline_data <- model_data %>%
    filter(year >= baseline_period[1] & year <= baseline_period[2])
  
  if(nrow(baseline_data) == 0) {
    warning(paste("Model", unique(model_data$model), "has no baseline period data"))
    return(NULL)
  }
  
  baseline_tas <- mean(baseline_data$annual_tas, na.rm = TRUE)
  
  # Calculate ΔT and ΔR% using Q10 formula: ΔR% = [Q₁₀^(ΔT/10) - 1] × 100%
  result <- model_data %>%
    mutate(
      delta_T = annual_tas - baseline_tas,
      delta_R_percent = (Q10_val^(delta_T/10) - 1) * 100,
      baseline_tas = baseline_tas,
      Q10_type = Q10_name
    )
  
  return(result)
}

# ----------------------------------------------------------------
# 4. Process All Model Data
# ----------------------------------------------------------------
cat("\n=== Processing All Model Data ===\n")

# Store annual temperature data for all models
all_annual_data <- list()

for(model_name in names(all_models)) {
  cat(paste0("\nProcessing model: ", model_name, "\n"))
  
  annual_data <- process_single_model_simple(model_name, all_models[[model_name]])
  
  if(!is.null(annual_data)) {
    all_annual_data[[model_name]] <- annual_data
    cat(paste0("  ✓ ", model_name, " processed: ", min(annual_data$year), "-", max(annual_data$year), 
               " (", nrow(annual_data), " years)\n"))
  } else {
    cat(paste0("  ✗ ", model_name, " processing failed\n"))
  }
}

# Check if we have sufficient data
if(length(all_annual_data) == 0) {
  stop("No valid model data found. Please check file paths and formats.")
}

cat(paste0("\n✓ Successfully processed ", length(all_annual_data), " models\n"))

# ----------------------------------------------------------------
# 5. Calculate ΔR% for Each Q10 Value
# ----------------------------------------------------------------
cat("\n=== Calculating ΔR% ===\n")

deltaR_results <- list()

for(Q10_name in names(Q10_values)) {
  cat(paste0("\nCalculating for Q10: ", Q10_name, " (Q₁₀ = ", Q10_values[[Q10_name]], ")\n"))
  
  Q10_val <- Q10_values[[Q10_name]]
  Q10_data_list <- list()
  
  for(model_name in names(all_annual_data)) {
    model_result <- calculate_deltaR_for_model(
      all_annual_data[[model_name]],
      Q10_val,
      Q10_name,
      baseline_period
    )
    
    if(!is.null(model_result)) {
      Q10_data_list[[model_name]] <- model_result
      cat(paste0("  ✓ ", model_name, "\n"))
    } else {
      cat(paste0("  ✗ ", model_name, " (calculation failed)\n"))
    }
  }
  
  if(length(Q10_data_list) > 0) {
    # Combine results from all models
    deltaR_results[[Q10_name]] <- bind_rows(Q10_data_list)
    cat(paste0("  ✓ ", Q10_name, ": ", length(Q10_data_list), " models calculated\n"))
  }
}

# Combine all results for boxplot analysis
all_deltaR_data <- bind_rows(deltaR_results)

# ----------------------------------------------------------------
# 6. Prepare Boxplot Data (2090-2100)
# ----------------------------------------------------------------
cat("\n=== Preparing Boxplot Data (2090-2100) ===\n")

# Filter data for 2090-2100 period
boxplot_data <- all_deltaR_data %>%
  filter(year >= boxplot_period[1] & year <= boxplot_period[2])

# Check data characteristics
cat(paste0("Boxplot data contains:\n"))
cat(paste0("  Time range: ", boxplot_period[1], "-", boxplot_period[2], "\n"))
cat(paste0("  Total data points: ", nrow(boxplot_data), "\n"))
cat(paste0("  Number of models: ", length(unique(boxplot_data$model)), "\n"))
cat(paste0("  Q10 types: ", paste(unique(boxplot_data$Q10_type), collapse = ", "), "\n"))

# Statistics by model and Q10 type
model_counts <- boxplot_data %>%
  group_by(model, Q10_type) %>%
  summarise(n_years = n(), .groups = "drop")

cat("\nNumber of years per model:\n")
print(model_counts)

# ----------------------------------------------------------------
# 7. Perform Statistical Analysis
# ----------------------------------------------------------------
cat("\n=== Performing Statistical Analysis ===\n")

# 7.1 Shapiro-Wilk normality test
cat("\n1. Shapiro-Wilk Normality Test:\n")

normality_test <- boxplot_data %>%
  group_by(Q10_type) %>%
  summarise(
    statistic = round(shapiro.test(delta_R_percent)$statistic, 4),
    p_value = round(shapiro.test(delta_R_percent)$p.value, 4),
    .groups = "drop"
  )

print(normality_test)

# Interpret normality test results
for(i in 1:nrow(normality_test)) {
  Q10_type <- normality_test$Q10_type[i]
  p_val <- normality_test$p_value[i]
  if(p_val > 0.05) {
    cat(paste0("  ", Q10_type, ": Data follows normal distribution (p = ", p_val, ")\n"))
  } else {
    cat(paste0("  ", Q10_type, ": Data does not follow normal distribution (p = ", p_val, 
               "), non-parametric tests recommended\n"))
  }
}

# 7.2 Wilcoxon rank-sum test (comparing two Q10 types)
cat("\n2. Wilcoxon Rank-Sum Test (comparing Plastisphere and Epiphyton):\n")

wilcox_test_result <- wilcox.test(
  delta_R_percent ~ Q10_type,
  data = boxplot_data,
  exact = FALSE,
  conf.int = TRUE
)

cat(paste0("  W = ", round(wilcox_test_result$statistic, 2), "\n"))
cat(paste0("  p-value = ", format.pval(wilcox_test_result$p.value, digits = 4), "\n"))
cat(paste0("  95% Confidence Interval: [", round(wilcox_test_result$conf.int[1], 2), ", ", 
           round(wilcox_test_result$conf.int[2], 2), "]\n"))
cat(paste0("  Median Difference: ", round(wilcox_test_result$estimate, 2), "\n"))

# 7.3 Effect size calculation (Cliff's delta)
cat("\n3. Effect Size Calculation (Cliff's delta):\n")

# Extract data for each group
plastisphere_data <- boxplot_data %>% 
  filter(Q10_type == "Plastisphere") %>% 
  pull(delta_R_percent)
epiphyton_data <- boxplot_data %>% 
  filter(Q10_type == "Epiphyton") %>% 
  pull(delta_R_percent)

# Calculate Cliff's delta effect size
cliffs_delta <- 2 * wilcox_test_result$statistic / (length(plastisphere_data) * length(epiphyton_data)) - 1

cat(paste0("  Cliff's delta = ", round(cliffs_delta, 3), "\n"))

# Interpret effect size magnitude
if(abs(cliffs_delta) < 0.147) {
  cat("  Effect size magnitude: Negligible\n")
} else if(abs(cliffs_delta) < 0.33) {
  cat("  Effect size magnitude: Small\n")
} else if(abs(cliffs_delta) < 0.474) {
  cat("  Effect size magnitude: Medium\n")
} else {
  cat("  Effect size magnitude: Large\n")
}

# 7.4 Descriptive statistics
cat("\n4. Descriptive Statistics (2090-2100):\n")

descriptive_stats <- boxplot_data %>%
  group_by(Q10_type) %>%
  summarise(
    n = n(),
    Mean = round(mean(delta_R_percent, na.rm = TRUE), 2),
    SD = round(sd(delta_R_percent, na.rm = TRUE), 2),
    Median = round(median(delta_R_percent, na.rm = TRUE), 2),
    Q1 = round(quantile(delta_R_percent, 0.25, na.rm = TRUE), 2),
    Q3 = round(quantile(delta_R_percent, 0.75, na.rm = TRUE), 2),
    Min = round(min(delta_R_percent, na.rm = TRUE), 2),
    Max = round(max(delta_R_percent, na.rm = TRUE), 2),
    IQR = round(IQR(delta_R_percent, na.rm = TRUE), 2),
    .groups = "drop"
  )

print(descriptive_stats)

# ----------------------------------------------------------------
# 8. Create Boxplot
# ----------------------------------------------------------------
cat("\n=== Creating Boxplot ===\n")

# Define colors for Q10 types
Q10_colors <- c(
  "Plastisphere" = "#FF6B6B",  # Red color
  "Epiphyton" = "#32CD32"      # Green color
)

# 8.1 Basic boxplot
boxplot_basic <- ggplot(boxplot_data, aes(x = Q10_type, y = delta_R_percent, fill = Q10_type)) +
  geom_boxplot(
    alpha = 0.8,
    outlier.shape = 16,
    outlier.size = 2,
    outlier.alpha = 0.6,
    width = 0.6
  ) +
  # Add individual data points (jittered)
  geom_jitter(
    width = 0.15,
    height = 0,
    size = 3,
    alpha = 0.5,
    color = "gray30"
  ) +
  scale_fill_manual(values = Q10_colors) +
  labs(
    title = "N₂O Emission Rate Increase ΔR% (2090-2100)",
    subtitle = paste("SSP370 Scenario |", length(unique(boxplot_data$model)), "Model Ensemble"),
    x = "Q₁₀ Type",
    y = "ΔR (%)",
    caption = paste(
      "Wilcoxon test: p =", format.pval(wilcox_test_result$p.value, digits = 3),
      "| Cliff's delta =", round(cliffs_delta, 3)
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    plot.caption = element_text(hjust = 0.5, size = 9, color = "gray50"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.2)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = expansion(mult = 0.05),
    limits = c(-5, 75)
  ) +
  # Add statistical significance annotation
  geom_signif(
    comparisons = list(c("Plastisphere", "Epiphyton")),
    map_signif_level = TRUE,
    textsize = 4,
    vjust = -0.2,
    tip_length = 0.01
  ) +
  # Add mean markers (diamonds)
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 5,
    color = "black"
  )

# Display the boxplot
print(boxplot_basic)

# Save the boxplot as a narrow/tall PDF
ggsave(
  filename = "DEN_narrow_boxplot.pdf",  # Output filename
  plot = boxplot_basic,                  # Plot object to save
  width = 3,                             # Width in inches (narrow)
  height = 6,                            # Height in inches (tall)
  units = "in",                          # Units for dimensions
  dpi = 300                              # Resolution for publication
)

# ----------------------------------------------------------------
# 9. Calculate Ensemble Statistics (for main trend plot)
# ----------------------------------------------------------------
cat("\n=== Calculating Ensemble Statistics ===\n")

ensemble_stats <- list()

for(Q10_name in names(deltaR_results)) {
  Q10_data <- deltaR_results[[Q10_name]]
  
  # Calculate ensemble statistics by year
  stats_by_year <- Q10_data %>%
    group_by(year, Q10_type) %>%
    summarise(
      delta_R_mean = mean(delta_R_percent, na.rm = TRUE),
      delta_R_sd = sd(delta_R_percent, na.rm = TRUE),
      delta_R_median = median(delta_R_percent, na.rm = TRUE),
      delta_R_min = min(delta_R_percent, na.rm = TRUE),
      delta_R_max = max(delta_R_percent, na.rm = TRUE),
      delta_R_q025 = quantile(delta_R_percent, 0.025, na.rm = TRUE),
      delta_R_q975 = quantile(delta_R_percent, 0.975, na.rm = TRUE),
      n_models = n_distinct(model),
      .groups = "drop"
    )
  
  ensemble_stats[[Q10_name]] <- stats_by_year
  cat(paste0("  ✓ ", Q10_name, " statistics calculated\n"))
}

# Combine statistics for all Q10 types
all_stats <- bind_rows(ensemble_stats)

# Filter for 2015-2100 period
all_stats <- all_stats %>% filter(year >= 2015 & year <= 2100)

cat(paste0("\n✓ Statistics calculated, data range: ", min(all_stats$year), "-", 
           max(all_stats$year), "\n"))

# ----------------------------------------------------------------
# 10. Create Main Trend Plot
# ----------------------------------------------------------------
cat("\n=== Creating Main Trend Plot ===\n")

# Create main trend plot
main_plot <- ggplot(all_stats, aes(x = year, color = Q10_type, fill = Q10_type)) +
  # Add baseline period background shading
  annotate("rect",
           xmin = baseline_period[1], xmax = baseline_period[2],
           ymin = -Inf, ymax = Inf,
           alpha = 0.1, fill = "gray50") +
  # Add boxplot analysis period background shading
  annotate("rect",
           xmin = boxplot_period[1], xmax = boxplot_period[2],
           ymin = -Inf, ymax = Inf,
           alpha = 0.05, fill = "blue") +
  # Plot 95% confidence interval as ribbon
  geom_ribbon(aes(ymin = delta_R_q025, ymax = delta_R_q975),
              alpha = 0.2, color = NA) +
  # Plot mean line
  geom_line(aes(y = delta_R_mean), linewidth = 1.2) +
  # Set colors and fills
  scale_color_manual(values = Q10_colors, name = "Q₁₀ Type") +
  scale_fill_manual(values = Q10_colors, name = "Q₁₀ Type") +
  # Plot labels
  labs(
    title = "N₂O Emission Rate Increase ΔR% (2015-2100)",
    subtitle = paste(
      "SSP370 Scenario |", length(all_annual_data), "Model Ensemble",
      paste("Models:", paste(names(all_annual_data), collapse = ", ")),
      paste("Baseline:", baseline_period[1], "-", baseline_period[2]),
      paste("Boxplot analysis:", boxplot_period[1], "-", boxplot_period[2]),
      sep = "\n"
    ),
    x = "Year",
    y = "ΔR (%)",
    caption = paste(
      "Formula: ΔR% = [Q₁₀^(ΔT/10) - 1] × 100%",
      "Shaded area: 95% confidence interval",
      paste0("Plastisphere: Q₁₀ = ", Q10_values$Plastisphere),
      paste0("Epiphyton: Q₁₀ = ", Q10_values$Epiphyton),
      sep = "\n"
    )
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
    plot.caption = element_text(hjust = 0, size = 8, color = "gray50"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray92", linewidth = 0.2)
  ) +
  scale_x_continuous(
    breaks = seq(2020, 2100, 10),
    limits = c(2015, 2100)
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(-5, 60)
  ) +
  # Add boxplot analysis period label
  annotate("text",
           x = mean(boxplot_period),
           y = max(all_stats$delta_R_mean) * 0.95,
           label = paste("Boxplot analysis\n", boxplot_period[1], "-", boxplot_period[2]),
           color = "blue",
           size = 3,
           alpha = 0.7)

# Display the main plot
print(main_plot)

# Save the main trend plot
ggsave(
  filename = "DEN_DeltaR_Main_Trend_2015-2100.pdf",  # Output filename
  plot = main_plot,                                   # Plot object to save
  width = 6,                                          # Width in inches
  height = 6,                                         # Height in inches
  dpi = 300                                           # Resolution
)
cat("✓ Main trend plot saved\n")

# ----------------------------------------------------------------
# 11. Save All Data and Results
# ----------------------------------------------------------------
cat("\n=== Saving All Data and Results ===\n")

# Save detailed data
write_csv(all_deltaR_data, "DEN_DeltaR_Detailed_Data_All_Years.csv")
write_csv(boxplot_data, "DEN_DeltaR_Boxplot_Data_2090-2100.csv")
write_csv(all_stats, "DEN_DeltaR_Trend_Statistics.csv")

# Save statistical results
stat_results <- list(
  Normality_Test = normality_test,
  Wilcoxon_Test = data.frame(
    Test = "Wilcoxon rank-sum test",
    W_Statistic = wilcox_test_result$statistic,
    P_Value = wilcox_test_result$p.value,
    Estimate = wilcox_test_result$estimate,
    CI_Lower = wilcox_test_result$conf.int[1],
    CI_Upper = wilcox_test_result$conf.int[2]
  ),
  Effect_Size = data.frame(
    Measure = "Cliff's delta",
    Value = cliffs_delta
  ),
  Descriptive_Stats = descriptive_stats
)

# Save statistical results as CSV files
write_csv(normality_test, "DEN_Statistical_Normality_Test.csv")
write_csv(descriptive_stats, "DEN_Statistical_Descriptive_Stats.csv")
write_csv(data.frame(
  Test = "Wilcoxon rank-sum test",
  W_Statistic = wilcox_test_result$statistic,
  P_Value = wilcox_test_result$p.value,
  Estimate = wilcox_test_result$estimate,
  CI_Lower = wilcox_test_result$conf.int[1],
  CI_Upper = wilcox_test_result$conf.int[2],
  Cliffs_Delta = cliffs_delta
), "DEN_Statistical_Hypothesis_Tests.csv")

cat("✓ All data and results saved successfully\n")
cat("=== Analysis Complete ===\n")