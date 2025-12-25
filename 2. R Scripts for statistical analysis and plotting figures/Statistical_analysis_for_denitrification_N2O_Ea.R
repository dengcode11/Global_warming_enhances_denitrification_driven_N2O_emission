# Load necessary packages
library(tidyr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(effsize)
library(coin)
library(car)

# ==================== PART 1: DATA INPUT ====================
# Users can replace with their own data here
# Format requirement: First column is the grouping variable, subsequent columns are variables to compare
data <- data.frame(
  Group = rep(c("Group1", "Group2"), each = 3),  # Grouping variable
  Variable1 = c(2.281, 2.241, 2.258, 2.536, 2.563, 2.534),  # First variable
  Variable2 = c(2.178, 2.180, 2.166, 2.357, 2.207, 2.240)   # Second variable
)

# Display the data
cat("Raw data:\n")
print(data)

# ==================== PART 2: INDEPENDENT SAMPLES COMPARISON (Differences between groups for the same variable) ====================

# Define a general analysis function for independent comparisons
analyze_independent_comparison <- function(data, var_name, var_label) {
  cat(paste("\n=== ", var_label, " (", var_name, ") Between-Group Comparison ===\n", sep=""))
  
  # Extract data for the current variable
  var_data <- data[, c("Group", var_name)]
  names(var_data) <- c("Group", "Value")
  
  # Descriptive statistics
  desc_stats <- var_data %>%
    group_by(Group) %>%
    summarise(
      n = n(),
      Mean = mean(Value),
      SD = sd(Value),
      SE = sd(Value)/sqrt(n()),
      Median = median(Value),
      Min = min(Value),
      Max = max(Value),
      CV = (sd(Value)/mean(Value))*100  # Coefficient of variation
    )
  
  cat("Descriptive Statistics:\n")
  print(as.data.frame(desc_stats))
  
  # Normality test (Shapiro-Wilk)
  normality_test <- var_data %>%
    group_by(Group) %>%
    summarise(
      Shapiro_p = shapiro.test(Value)$p.value
    )
  
  cat("\nNormality Test (Shapiro-Wilk):\n")
  print(normality_test)
  
  # Homogeneity of variance test (Levene's test)
  var_test <- leveneTest(Value ~ Group, data = var_data)
  cat("\nHomogeneity of Variance Test (Levene's Test):\n")
  print(var_test)
  
  # Select statistical method based on test results
  if(all(normality_test$Shapiro_p > 0.05) && var_test$`Pr(>F)`[1] > 0.05) {
    # If normality and homogeneity of variance are met, use t-test
    cat("\nUsing Independent Samples t-test:\n")
    t_test_result <- t.test(Value ~ Group, data = var_data, var.equal = TRUE)
    print(t_test_result)
    
    # Calculate effect size (Cohen's d)
    cohen_d_result <- cohens_d(Value ~ Group, data = var_data)
    cat("\nEffect Size (Cohen's d):\n")
    print(cohen_d_result)
    
    test_method <- "Independent t-test"
    p_value <- t_test_result$p.value
    effect_size <- cohen_d_result$effsize
  } else {
    # If parametric assumptions are not met, use Mann-Whitney U test
    cat("\nUsing Mann-Whitney U Test (Wilcoxon rank-sum test):\n")
    wilcox_result <- wilcox.test(Value ~ Group, data = var_data, exact = FALSE)
    print(wilcox_result)
    
    # Calculate effect size (Cliff's delta)
    cliff_delta_result <- cliff.delta(Value ~ Group, data = var_data)
    cat("\nEffect Size (Cliff's delta):\n")
    print(cliff_delta_result)
    
    test_method <- "Mann-Whitney U test"
    p_value <- wilcox_result$p.value
    effect_size <- cliff_delta_result$estimate
  }
  
  # Visualization
  p <- ggplot(var_data, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = NA) +
    geom_jitter(size = 3, shape = 21, width = 0.1, height = 0) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "white") +
    labs(title = paste("Between-Group Comparison of", var_label),
         x = "Group", 
         y = var_label) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12),
          legend.position = "none") +
    scale_fill_manual(values = c("Group1" = "#4E84C4", "Group2" = "#D16103"))
  
  print(p)
  
  # Return results summary
  result_summary <- data.frame(
    Variable = var_label,
    Test_Method = test_method,
    P_Value = p_value,
    Effect_Size = effect_size,
    Group1_Mean = desc_stats$Mean[desc_stats$Group == "Group1"],
    Group2_Mean = desc_stats$Mean[desc_stats$Group == "Group2"],
    Group1_SD = desc_stats$SD[desc_stats$Group == "Group1"],
    Group2_SD = desc_stats$SD[desc_stats$Group == "Group2"]
  )
  
  return(result_summary)
}

# Perform independent samples comparison for each variable
independent_results <- data.frame()

# Get all variable names to analyze (exclude the grouping column)
variables_to_analyze <- names(data)[-1]

for (i in 1:length(variables_to_analyze)) {
  var_name <- variables_to_analyze[i]
  var_label <- paste("Variable", i)  # User can modify labels as needed
  independent_results <- rbind(independent_results, 
                               analyze_independent_comparison(data, var_name, var_label))
}

# Output comprehensive results for independent samples comparison
cat("\n\n=== Comprehensive Independent Samples Comparison Results ===\n")
print(independent_results)

# ==================== PART 3: PAIRED SAMPLES COMPARISON (Differences between variables within the same group) ====================

# Define function for paired comparisons
analyze_paired_comparison <- function(data, group_name, var_names, var_labels) {
  cat(paste("\n\n=== Paired Comparison Within ", group_name, " Group ===\n", sep=""))
  
  # Extract data for the current group
  group_data <- data[data$Group == group_name, ]
  
  # Check if there are enough variables for paired comparison
  if (length(var_names) < 2) {
    cat("Error: At least two variables are required for paired comparison\n")
    return(NULL)
  }
  
  # Perform all possible paired comparisons
  paired_results <- data.frame()
  
  for (i in 1:(length(var_names)-1)) {
    for (j in (i+1):length(var_names)) {
      var1 <- var_names[i]
      var2 <- var_names[j]
      label1 <- var_labels[i]
      label2 <- var_labels[j]
      
      cat(paste("\n--- ", label1, " vs ", label2, " ---\n", sep=""))
      
      # Extract paired data
      paired_values <- data.frame(
        Var1 = group_data[, var1],
        Var2 = group_data[, var2]
      )
      
      # Descriptive statistics
      desc_stats <- data.frame(
        Variable = c(label1, label2),
        Mean = c(mean(paired_values$Var1), mean(paired_values$Var2)),
        SD = c(sd(paired_values$Var1), sd(paired_values$Var2)),
        Median = c(median(paired_values$Var1), median(paired_values$Var2))
      )
      
      cat("Descriptive Statistics:\n")
      print(desc_stats)
      
      # Calculate paired differences
      differences <- paired_values$Var1 - paired_values$Var2
      
      # Normality test for paired differences
      shapiro_test <- shapiro.test(differences)
      cat("\nNormality Test for Paired Differences (Shapiro-Wilk):\n")
      print(shapiro_test)
      
      # Select statistical method based on normality test
      if(shapiro_test$p.value > 0.05) {
        # If normality is met, use paired t-test
        cat("\nUsing Paired t-test:\n")
        paired_t_result <- t.test(paired_values$Var1, paired_values$Var2, 
                                  paired = TRUE)
        print(paired_t_result)
        
        # Calculate effect size (Cohen's d for paired samples)
        mean_diff <- mean(differences)
        sd_diff <- sd(differences)
        cohen_d_paired <- mean_diff / sd_diff
        
        test_method <- "Paired t-test"
        p_value <- paired_t_result$p.value
        effect_size <- cohen_d_paired
      } else {
        # If normality is not met, use Wilcoxon signed-rank test
        cat("\nUsing Wilcoxon Signed-Rank Test:\n")
        wilcox_paired_result <- wilcox.test(paired_values$Var1, paired_values$Var2, 
                                            paired = TRUE)
        print(wilcox_paired_result)
        
        # Calculate effect size (matched rank biserial correlation)
        n_pairs <- nrow(paired_values)
        w_stat <- wilcox_paired_result$statistic
        r_effect <- abs(w_stat/(n_pairs*(n_pairs+1)/2) - 0.5) * 2
        
        test_method <- "Wilcoxon signed-rank test"
        p_value <- wilcox_paired_result$p.value
        effect_size <- r_effect
      }
      
      # Create results summary
      result_summary <- data.frame(
        Group = group_name,
        Comparison = paste(label1, "vs", label2),
        Test_Method = test_method,
        P_Value = p_value,
        Effect_Size = effect_size,
        Mean_Var1 = desc_stats$Mean[1],
        Mean_Var2 = desc_stats$Mean[2],
        SD_Var1 = desc_stats$SD[1],
        SD_Var2 = desc_stats$SD[2]
      )
      
      paired_results <- rbind(paired_results, result_summary)
    }
  }
  
  return(paired_results)
}

# Perform paired comparisons for each group
paired_results_all <- data.frame()

# Get all unique groups
groups <- unique(data$Group)

# Define variable labels (user should modify according to their data)
variable_labels <- c("Variable 1", "Variable 2")  # Should correspond to variable names

for (group in groups) {
  group_paired_results <- analyze_paired_comparison(data, group, variables_to_analyze, variable_labels)
  if (!is.null(group_paired_results)) {
    paired_results_all <- rbind(paired_results_all, group_paired_results)
  }
}

# Output comprehensive results for paired comparisons
if (nrow(paired_results_all) > 0) {
  cat("\n\n=== Comprehensive Paired Comparison Results ===\n")
  print(paired_results_all)
}

# ==================== PART 4: RESULTS SUMMARY ====================

cat("\n\n=================== Statistical Analysis Results Summary ===================\n")

# Add significance markers for independent results
if (nrow(independent_results) > 0) {
  independent_results$Significance <- ifelse(independent_results$P_Value < 0.001, "***",
                                             ifelse(independent_results$P_Value < 0.01, "**",
                                                    ifelse(independent_results$P_Value < 0.05, "*", "ns")))
  
  cat("\n1. Independent Samples Comparison Results:\n")
  cat("Note: *** p<0.001, ** p<0.01, * p<0.05, ns not significant\n\n")
  
  summary_independent <- independent_results[, c("Variable", "Test_Method", 
                                                 "Group1_Mean", "Group2_Mean", 
                                                 "Group1_SD", "Group2_SD", 
                                                 "P_Value", "Effect_Size", "Significance")]
  
  names(summary_independent) <- c("Variable", "Test Method", "Group1 Mean", "Group2 Mean", 
                                  "Group1 SD", "Group2 SD", "P Value", "Effect Size", "Significance")
  
  print(summary_independent, row.names = FALSE)
}

# Add significance markers for paired results
if (nrow(paired_results_all) > 0) {
  paired_results_all$Significance <- ifelse(paired_results_all$P_Value < 0.001, "***",
                                            ifelse(paired_results_all$P_Value < 0.01, "**",
                                                   ifelse(paired_results_all$P_Value < 0.05, "*", "ns")))
  
  cat("\n\n2. Paired Comparison Results:\n")
  cat("Note: *** p<0.001, ** p<0.01, * p<0.05, ns not significant\n\n")
  
  summary_paired <- paired_results_all[, c("Group", "Comparison", "Test_Method",
                                           "Mean_Var1", "Mean_Var2",
                                           "SD_Var1", "SD_Var2",
                                           "P_Value", "Effect_Size", "Significance")]
  
  names(summary_paired) <- c("Group", "Comparison", "Test Method", 
                             "Variable 1 Mean", "Variable 2 Mean", 
                             "Variable 1 SD", "Variable 2 SD", 
                             "P Value", "Effect Size", "Significance")
  
  print(summary_paired, row.names = FALSE)
}

# ==================== PART 5: MULTIVARIATE VISUALIZATION ====================

cat("\n\n=================== Multivariate Data Visualization ===================\n")

# Convert data to long format for visualization
data_long <- pivot_longer(data, 
                          cols = all_of(variables_to_analyze),
                          names_to = "Variable",
                          values_to = "Value")

# Add variable labels (if needed)
data_long$VariableLabel <- factor(data_long$Variable,
                                  levels = variables_to_analyze,
                                  labels = variable_labels)

# Create grouped boxplots
p_group <- ggplot(data_long, aes(x = VariableLabel, y = Value, fill = Group)) +
  geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = NA) +
  geom_jitter(size = 2, shape = 21, width = 0.1, height = 0, alpha = 0.7) +
  facet_wrap(~ VariableLabel, scales = "free") +
  labs(title = "Multivariate Group Comparison",
       x = "Variable", 
       y = "Measurement Value") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 10)) +
  scale_fill_manual(values = c("Group1" = "#4E84C4", "Group2" = "#D16103"))

print(p_group)

# Create paired data visualizations (for within-group paired comparisons)
for (group in groups) {
  group_data_long <- data_long[data_long$Group == group, ]
  
  p_paired <- ggplot(group_data_long, aes(x = VariableLabel, y = Value, group = 1)) +
    geom_line(aes(group = interaction(rep(1:nrow(data)/2, each=2))), 
              color = "gray", alpha = 0.5) +
    geom_point(size = 3, aes(color = VariableLabel)) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 5, color = "red") +
    labs(title = paste("Paired Comparison Within", group, "Group"),
         x = "Variable", 
         y = "Measurement Value") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_text(size = 12))
  
  print(p_paired)
}

# ==================== PART 6: SAVE RESULTS ====================

# Create timestamp for file naming
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Save independent comparison results
write.csv(independent_results, 
          paste0("independent_comparison_results_", timestamp, ".csv"), 
          row.names = FALSE)

# Save paired comparison results
if (nrow(paired_results_all) > 0) {
  write.csv(paired_results_all, 
            paste0("paired_comparison_results_", timestamp, ".csv"), 
            row.names = FALSE)
}

# Save summary results
if (exists("summary_independent") && exists("summary_paired")) {
  all_summary <- list(
    Independent_Comparison = summary_independent,
    Paired_Comparison = summary_paired
  )
  
  # Save as Excel format (install writexl package if needed)
  # library(writexl)
  # write_xlsx(all_summary, paste0("all_results_summary_", timestamp, ".xlsx"))
  
  # Save as CSV
  write.csv(summary_independent, 
            paste0("summary_independent_", timestamp, ".csv"), 
            row.names = FALSE)
  
  if (nrow(summary_paired) > 0) {
    write.csv(summary_paired, 
              paste0("summary_paired_", timestamp, ".csv"), 
              row.names = FALSE)
  }
}

cat(paste("\n\nAnalysis completed! Results saved to the current working directory.\nTimestamp:", timestamp, "\n"))