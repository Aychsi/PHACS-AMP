
# This file has helper functions for AMP codes.


library("skimr")
library("dplyr")
library("RcmdrMisc")
library("haven")
library("sjPlot")
library(ggplot2)
library(purrr)
library(tidyr)
library(ggrepel)
library(lme4)
library(LDATS)
library(dplyr)
library(BioAge)
library(Metrics)
library(ggpubr)
library(broom)
library(performance)
library(nnet)
library(pscl)
library(e1071)
library(gridExtra)


# Standardize a column or a list of columns in a dataframe
# Used to standardize poa and kdm
# Standardize selected columns grouped by "sex"
standardize_columns <- function(data, column_names, group_var) {
  if (!is.vector(column_names)) {
    column_names <- list(column_names)
  }
  
  # Group the data by the specified grouping variable
  data <- data %>%
    group_by({{ group_var }})
  
  for (col in column_names) {
    # Standardize within each group
    data <- data %>%
      mutate(!!col := (.data[[col]] - mean(.data[[col]], na.rm = TRUE)) / 
               sd(.data[[col]], na.rm = TRUE))
  }
  
  # Remove the grouping information
  data <- ungroup(data)
  
  return(data)
}


# This function filters the dataframes by the specified week numbers, pivots them into wide format, 
# and renames the columns based on the timepoints. The result is a list of wide-format dataframes, 
# one for each input dataframe.

wideFormatByWeek <- function(df, week_nums, timepoint1, timepoint2, values_to_pivot) {
  tryCatch({
    # Filter the dataframe by the specified week numbers
    filtered_df <- df %>%
      filter(week %in% week_nums)
    
    # Pivot the specified columns from long to wide format
    wide_df <- filtered_df %>%
      pivot_wider(names_from = week, values_from = {{values_to_pivot}})
    
    # Unlist the pivoted columns and fill with NA if NULL
    wide_df <- wide_df %>%
      unnest(everything())
    
    return(wide_df)
  }, error = function(e) {
    # Print an error message and return NULL if an error occurs
    cat("Error message:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Example usage:
# Assuming you have a dataframe df, specified week numbers week_nums,
# and timepoints timepoint1 and timepoint2, and values_to_pivot is a column name or a vector of column names to pivot
#wide_df <- wideFormatByWeek(df, week_nums = c(1, 2), timepoint1 = "Baseline", timepoint2 = "FollowUp", values_to_pivot = c("Column1", "Column2"))



# This function takes two arguments:
#   
# df: The dataframe you want to operate on.
# string_to_match: The string to look for at the end of column names.
# result_col_name: Column name of the resulting sum scores.
# 
# It then uses grep to find column names containing the specified string, and c_across from the 
#dplyr package to sum those columns for each row. The result is added as a new column with a name 
#like "Sum_string_to_match".

rowSumColumnsWithSuffix <- function(df, string_to_match, result_col_name) {
  # Find column names that end with the specified string
  matching_columns <- grep(paste0(string_to_match, "$"), names(df), value = TRUE)
  
  # Extract the matching columns and calculate row sums
  df <- df %>%
    rowwise() %>%
    mutate(!!result_col_name := if(all(is.na(c_across(all_of(matching_columns))))) NA_real_ else
      sum(c_across(all_of(matching_columns)), na.rm = TRUE))
  
  return(df)
}


# This function takes two arguments:
#   
# df: The dataframe you want to operate on.
# columns_to_match: list of column names to look for at the end of column names.
# result_col_name: Column name of the resulting sum scores.
# 
# It then uses grep to find column names containing the specified string, and c_across from the 
#dplyr package to sum those columns for each row. The result is added as a new column with a name 
#like "Sum_string_to_match".

rowSumColumnsWithNames <- function(df, columns_to_match, result_col_name) {
  # Calculate row sums for the specified columns
  df <- df %>%
    rowwise() %>%
    mutate(!!result_col_name := if(all(is.na(c_across(all_of(columns_to_match))))) NA_real_ else
      sum(c_across(all_of(columns_to_match)), na.rm = TRUE))
  
  return(df)
}


# Function Name: getOutcomeNames
# 
# Input: A wide-format dataframe with columns ending in "_96" or "_48".
# 
# Output: A vector of unique outcome names without the week numbers.
# 
# Description: This function extracts and returns unique outcome names from the input dataframe's 
# column names, specifically those ending with "_96" or "_48".

getOutcomeNames <- function(wide_df) {
  # Extract column names
  col_names <- names(wide_df)
  
  # Extract outcome names by looking for columns ending with "_96" or "_48"
  outcome_names <- col_names[grep("_(96|48)$", col_names)]
  
  # Remove duplicates and return the unique outcome names
  unique_outcome_names <- unique(sub("_(96|48)$", "", outcome_names))
  return(unique_outcome_names)
}

# Same as above, but specified weeks
getOutcomeNamesCustom <- function(wide_df, suffixes) {
  # Convert the suffixes to a single regex pattern
  suffix_pattern <- paste0("(_(", paste(suffixes, collapse="|"), ")$)")
  
  # Extract column names
  col_names <- names(wide_df)
  
  # Extract outcome names using the suffix pattern
  outcome_names <- col_names[grep(suffix_pattern, col_names)]
  
  # Remove duplicates and return the unique outcome names
  unique_outcome_names <- unique(sub(suffix_pattern, "", outcome_names))
  return(unique_outcome_names)
}


# Function Summary:
# Analyzes any outcome data for a specific eventname to check which outcomes are zero-inflated.

# Input:
#   - data: The dataframe containing biomarker and outcome data.
#   - outcome_list: A list of outcome variable names to analyze.
#   - formula: The formula for the regression analysis.

# Output:
#   - Returns a list with two sets of results:
#     1. Fits for variables with low zero-inflation (fits_reg).
#     2. Fits for variables with high zero-inflation (fits_zero).
#     Also returns corresponding variable names (fits_reg_name and fits_zero_name).

# Usage:
#   result <- analyze_cbcl_data(data, cbcl_list, eventname, formula)

check_zero_inflation <- function(my_data, outcome_list, my_formula) {
  fits_reg <- list()
  fits_reg_name <- list()
  fits_zero <- list()
  fits_zero_name <- list()
  
  for (outcome in outcome_list) { 
    print(outcome)
    # Check if there are any zero values in the column
    zero_check <- all(my_data[[outcome]] != 0, na.rm = TRUE)
    
    # If there are no zero values, subtract the smallest value
    if (zero_check) {
      print("no Zeroes")
      min_value <- min(my_data[[outcome]], na.rm = TRUE)
      my_data[[outcome]] <- my_data[[outcome]] - min_value
    }
    
    outcome_model <- glm(my_data, formula = as.formula(paste(outcome, "~", my_formula)), family = poisson)
    
    
    if (check_zeroinflation(outcome_model)$ratio < 0.05) {
      fits_zero <- rbind(fits_zero, list(outcome, check_zeroinflation(outcome_model)$ratio))
      fits_zero_name <- append(fits_zero_name, outcome)
    } else {
      fits_reg <- rbind(fits_reg, list(outcome, check_zeroinflation(outcome_model)$ratio))
      fits_reg_name <- append(fits_reg_name, outcome)
    }
  }
  
  return(list(fits_reg = fits_reg, fits_reg_name = fits_reg_name, fits_zero = fits_zero, fits_zero_name = fits_zero_name))
}




# This function fit_and_save_lmer takes the dataframe df, the formula, and the list of outcome names 
# outcome_list as inputs. It fits linear mixed-effects models using lmer for each outcome in the 
# list and stores the coefficient estimates and p-values in a list called results. 
# You can then access the results for each outcome and print or process them as needed.
# last_time: usually 96 weeks (should be a string)
# first_time: Either 0 or 48 weeks (should be a string)



fit_lm_getOutcomeMetrics <- function(dataframe, outcome_list, formula_name, last_time, first_time) {
  results <- list()
  
  for (name in outcome_list) { 
    tryCatch({
      name_last <- paste0(name, "_", last_time)
      
      # Check the proportion of NA values in the outcome variable
      na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
      
      if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
        # Check if all values in the column are 0
        if (all(dataframe[[name_last]] == 0)) {
          print(paste(name_last, "contains all 0's. Skipping."))
          next  # Go to the next iteration
        }
        
        name_base <- paste0(name, "_", first_time)
        
        print(paste(name_last, "~", formula_name, "+", name_base))
        
        mod <- lm(data = dataframe, 
                  formula = as.formula(paste(name_last, "~", formula_name, "+", name_base)),
                  na.action = na.exclude)  # This will exclude NA values
        
        # Predict the outcome variable using the model
        predicted <- predict(mod, newdata = dataframe)
        
        rmse1 <- function(o, p, m = TRUE) {
          sqrt(mean((o - p)^2, na.rm = m))
        }
        
        mae1 <- function(o, p, m = TRUE) {
          mean(abs(o - p), na.rm = m)
        }
        
        
        # Calculate MAE1 and RMSE1, removing NA values
        mae_value <- mae1(dataframe[[name_last]], predicted)
        rmse_value <- rmse1(dataframe[[name_last]], predicted)
        
        
        # Extract the R-squared value
        r_squared <- summary(mod)$r.squared
        
        # Get coefficients summary
        coeff_summary <- summary(mod)$coefficients
        
        # Extract p-values for independent variables (excluding intercept)
        p_values <- coeff_summary[-1, 4]  # Excludes the intercept
        
        # Create a data frame with R-squared, MAE1, RMSE1, and p-values
        fit_summary <- data.frame(
          R_Squared = r_squared, 
          MAE1 = mae_value, 
          RMSE1 = rmse_value,
          P_Values = I(list(p_values))  # Store p-values as a list column in the data frame
        )
        
        results[[name_last]] <- fit_summary
      } else {
        print(paste("Skipping", name_last, "due to more than 75% NA values"))
      }
    }, error = function(e) {
      print(paste("Error in modeling", name_last))
      print(e)
    })
  }
  
  return(results)
}

# Function 1: Run the model for each name in outcome_list
run_lin_models <- function(dataframe, outcome_list, formula_name, last_time, first_time) {
  model_list <- list()
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Check the proportion of NA values in the outcome variable
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
      # Check if all values in the column are 0
      if (all(dataframe[[name_last]] == 0)) {
        print(paste(name_last, "contains all 0's. Skipping."))
        next  # Go to the next iteration
      }
      
      mod <- lm(data = dataframe, 
                formula = as.formula(paste(name_last, "~", formula_name, "+", name_base)),
                na.action = na.exclude)  # This will exclude NA values
      model_list[[name_last]] <- mod
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(model_list)
}



# Function 2: Extract RMSE, MAE, and R-squared from each model
get_lin_model_metrics <- function(dataframe, model_list) {
  metric_list <- list()
  
  # Helper functions for RMSE and MAE
  rmse1 <- function(o, p, m = TRUE) {
    sqrt(mean((o - p)^2, na.rm = m))
  }
  
  mae1 <- function(o, p, m = TRUE) {
    mean(abs(o - p), na.rm = m)
  }
  
  for (name_last in names(model_list)) {
    mod <- model_list[[name_last]]
    
    # Predict using the model
    predicted <- predict(mod, newdata = dataframe)
    
    # Actual outcomes from the data
    actual_outcome <- dataframe[[name_last]]
    
    # Calculate MAE and RMSE
    mae_value <- mae1(actual_outcome, predicted)
    rmse_value <- rmse1(actual_outcome, predicted)
    
    # Extract the R-squared value
    r_squared <- summary(mod)$r.squared
    
    # Combine the metrics into a data frame
    metrics <- data.frame(R_Squared = r_squared, MAE1 = mae_value, RMSE1 = rmse_value)
    metric_list[[name_last]] <- metrics
  }
  
  return(metric_list)
}


get_lin_model_partial_rsquared <- function(dataframe, model_list, variable_of_interest) {
  metric_list <- list()
  
  # Helper functions for RMSE and MAE
  rmse1 <- function(o, p, m = TRUE) {
    sqrt(mean((o - p)^2, na.rm = m))
  }
  
  mae1 <- function(o, p, m = TRUE) {
    mean(abs(o - p), na.rm = m)
  }
  
  for (name_last in names(model_list)) {
    mod <- model_list[[name_last]]
    
    # Calculate SSR for the full model
    ssr_with_variable <- sum(residuals(mod)^2)
    
    # Fit a model without the variable of interest and calculate its SSR
    mod_without_var <- lm(update(formula(mod), paste(". ~ . -", variable_of_interest)), data = dataframe)
    ssr_without_variable <- sum(residuals(mod_without_var)^2)
    
    # Calculate partial R-squared
    partial_r2 <- (ssr_without_variable - ssr_with_variable) / ssr_without_variable
    
    # Predict using the model
    predicted <- predict(mod, newdata = dataframe)
    
    # Actual outcomes from the data
    actual_outcome <- dataframe[[name_last]]
    
    # Calculate MAE and RMSE
    mae_value <- mae1(actual_outcome, predicted)
    rmse_value <- rmse1(actual_outcome, predicted)
    
    # Extract the R-squared value for the whole model
    r_squared <- summary(mod)$r.squared
    
    # Combine the metrics into a data frame
    metrics <- data.frame(R_Squared = r_squared, Partial_R_Squared = partial_r2, MAE1 = mae_value, RMSE1 = rmse_value)
    metric_list[[name_last]] <- metrics
  }
  
  return(metric_list)
}



# Function 3: Extract betas and p-values of independent variables from each model
get_betas_and_pvalues <- function(model_list) {
  results_list <- list()
  
  for (name_last in names(model_list)) {
    mod <- model_list[[name_last]]
    
    # Extract the coefficient estimates and their p-values
    betas <- summary(mod)$coefficients[, "Estimate"]
    p_values <- summary(mod)$coefficients[, "Pr(>|t|)"]
    
    # Exclude the intercept (if you want to keep it, remove these lines)
    betas <- betas[-1]
    p_values <- p_values[-1]
    
    # Store results in a dataframe and then in the list
    model_results <- data.frame(Betas = betas, P_Values = p_values)
    results_list[[name_last]] <- model_results
  }
  
  return(results_list)
}

get_betas_and_pvalues_zeroinfl <- function(model_list) {
  results_list <- list()
  
  for (name_last in names(model_list)) {
    mod <- model_list[[name_last]]
    summary_mod <- summary(mod)
    
    # Extract the coefficient estimates and their p-values for the count part
    betas_count <- summary_mod$count$coefficients[, "Estimate"]
    p_values_count <- summary_mod$count$coefficients[, "Pr(>|z|)"]
    
    # Extract the coefficient estimates and their p-values for the zero-inflation part
    betas_zero <- summary_mod$zero$coefficients[, "Estimate"]
    p_values_zero <- summary_mod$zero$coefficients[, "Pr(>|z|)"]
    
    # Exclude the intercept for both parts (if you want to keep it, remove these lines)
    betas_count <- betas_count[-1]
    p_values_count <- p_values_count[-1]
    betas_zero <- betas_zero[-1]
    p_values_zero <- p_values_zero[-1]
    
    # Store results in dataframes
    model_results_count <- data.frame(Betas_Count = betas_count, P_Values_Count = p_values_count)
    model_results_zero <- data.frame(Betas_Zero = betas_zero, P_Values_Zero = p_values_zero)
    
    # Combine the results
    combined_results <- cbind(model_results_count, model_results_zero)
    
    results_list[[name_last]] <- combined_results
  }
  
  return(results_list)
}



# Use sensemakr package to get partial r-squared
get_lin_model_metrics_sensemakr <- function(dataframe, model_list, variable_of_interest) {
  library(sensemakr)
  metric_list <- list()
  
  # Helper functions for RMSE and MAE
  rmse1 <- function(o, p, m = TRUE) {
    sqrt(mean((o - p)^2, na.rm = m))
  }
  
  mae1 <- function(o, p, m = TRUE) {
    mean(abs(o - p), na.rm = m)
  }
  
  for (name_last in names(model_list)) {
    mod <- model_list[[name_last]]
    
    # Conduct sensitivity analysis using sensemakr
    sensitivity_result <- sensemakr(model = mod, 
                                    treatment = variable_of_interest, 
                                    data = dataframe)
    
    # Extract partial R^2 from sensemakr result
    partial_r2 <- sensitivity_result$sensitivity_stats$r2yd.x
    
    # Predict using the model
    predicted <- predict(mod, newdata = dataframe)
    
    # Actual outcomes from the data
    actual_outcome <- dataframe[[name_last]]
    
    # Calculate MAE and RMSE
    mae_value <- mae1(actual_outcome, predicted)
    rmse_value <- rmse1(actual_outcome, predicted)
    
    # Extract the R-squared value for the whole model
    r_squared <- summary(mod)$r.squared
    
    # Combine the metrics into a data frame
    metrics <- data.frame(R_Squared = r_squared, 
                          Partial_R_Squared_sensemakr = partial_r2, 
                          MAE1 = mae_value, 
                          RMSE1 = rmse_value)
    metric_list[[name_last]] <- metrics
  }
  
  return(metric_list)
}


get_zip_model_metrics <- function(dataframe, model_list, variable_of_interest) {
  library(pscl)
  
  metric_list <- list()
  
  # Helper functions for RMSE and MAE
  rmse1 <- function(o, p, m = TRUE) {
    sqrt(mean((o - p)^2, na.rm = m))
  }
  
  mae1 <- function(o, p, m = TRUE) {
    mean(abs(o - p), na.rm = m)
  }
  
  for (name_last in names(model_list)) {
    mod <- model_list[[name_last]]
    
    # Predict using the model to get predicted counts
    predicted <- predict(mod, newdata = dataframe, type = "response")
    
    # Actual outcomes from the data
    actual_outcome <- dataframe[[name_last]]
    
    # Calculate MAE and RMSE
    mae_value <- mae1(actual_outcome, predicted)
    rmse_value <- rmse1(actual_outcome, predicted)
    
    # Combine the metrics into a data frame
    metrics <- data.frame(MAE1 = mae_value, 
                          RMSE1 = rmse_value)
    metric_list[[name_last]] <- metrics
  }
  
  return(metric_list)
}


# This function runs stepAIC for the three different methods of pace of aging
run_stepwise_model <- function(dataframe, outcome_list, last_time, first_time) {
  library(MASS)
  
  optimal_models <- list()
  
  formula_base <- " ~ poa_chang*sex + poa_belsky*sex + kdm_advance_0_kdm*sex + entryage + income + prace"
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Check the proportion of NA values in the outcome variable
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
      # Check if all values in the column are 0
      if (all(dataframe[[name_last]] == 0)) {
        print(paste(name_last, "contains all 0's. Skipping."))
        next  # Go to the next iteration
      }
      
      # Drop NA values for the entire dataset based on this model
      complete_data <- dataframe[complete.cases(dataframe[, all.vars(as.formula(paste(name_last, formula_base, "+", name_base)))]), ]
      
      # Print the number of rows removed
      rows_removed <- nrow(dataframe) - nrow(complete_data)
      print(paste("For model", name_last, ",", rows_removed, "rows were removed due to missing values."))
      
      initial_model <- lm(data = complete_data, 
                          formula = as.formula(paste(name_last, formula_base, "+", name_base)))
      
      # Use stepAIC to get the optimal model based on AIC but ensure that certain predictors always remain
      scope_list <- list(lower = as.formula(paste("~ sex + entryage + income + prace +", name_base)),
                         upper = as.formula(paste(name_last, formula_base, "+", name_base)))
      
      optimal_model <- MASS::stepAIC(initial_model, direction = "both", scope = scope_list, trace = FALSE)  # trace=FALSE suppresses output
      
      # Extract the right-hand side of the model's formula
      selected_vars <- optimal_model$coefficients   # Exclude the outcome variable (it's the 1st item)
      optimal_models[[name_last]] <- selected_vars
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(optimal_models)
}

run_stepwise_binlog_model <- function(dataframe, outcome_list, last_time, first_time) {
  library(MASS)
  
  optimal_models <- list()
  
  formula_base <- " ~ poa_chang*sex + poa_belsky*sex + kdm_advance_0_kdm*sex + entryage + income + prace"
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Check the proportion of NA values in the outcome variable
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
      
      # Ensure the outcome variable is numeric
      dataframe[[name_last]] <- as.numeric(as.character(dataframe[[name_last]]))
      
      tryCatch({
        
        # Drop NA values for the entire dataset based on this model
        complete_data <- dataframe[complete.cases(dataframe[, all.vars(as.formula(paste(name_last, formula_base, "+", name_base)))]), ]
        
        # Print the number of rows removed
        rows_removed <- nrow(dataframe) - nrow(complete_data)
        print(paste("For model", name_last, ",", rows_removed, "rows were removed due to missing values."))
        
        # Oversample the minority class
        oversampled_data <- oversample_minority(complete_data, name_last)
        
        initial_model <- glm(data = oversampled_data, 
                             formula = as.formula(paste(name_last, formula_base, "+", name_base)),
                             family = "binomial")
        
        # Use stepAIC to get the optimal model based on AIC but ensure that certain predictors always remain
        scope_list <- list(lower = as.formula(paste("~ sex + entryage + income + prace +", name_base)),
                           upper = as.formula(paste(name_last, formula_base, "+", name_base)))
        
        optimal_model <- MASS::stepAIC(initial_model, direction = "both", scope = scope_list, trace = FALSE)  # trace=FALSE suppresses output
        
        # Extract the right-hand side of the model's formula
        selected_vars <- optimal_model$coefficients
        optimal_models[[name_last]] <- selected_vars
        
      }, 
      error = function(e) {
        if(all(dataframe[[name_last]] == 0)) {
          print(paste(name_last, "contains all 0's. Skipping."))
        } else {
          print(paste("Error encountered with", name_last, ":", e$message))
        }
      })
      
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(optimal_models)
}




run_stepwise_zip_model <- function(dataframe, outcome_list, last_time, first_time) {
  library(MASS)
  library(pscl)
  
  optimal_models <- list()
  
  formula_base <- " ~ poa_chang*sex + poa_belsky*sex + kdm_advance_0_kdm*sex + entryage + income + prace"
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Check the proportion of NA values in the outcome variable
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
      
      # Adjust the outcome variable to have a minimum value of 0
      dataframe[[name_last]] <- dataframe[[name_last]] - min(dataframe[[name_last]], na.rm = TRUE)
      
      # Check if all values in the column are 0
      if (all(dataframe[[name_last]] == 0)) {
        print(paste(name_last, "contains all 0's. Skipping."))
        next  # Go to the next iteration
      }
      
      # Drop NA values for the entire dataset based on this model
      complete_data <- dataframe[complete.cases(dataframe[, all.vars(as.formula(paste(name_last, formula_base, "+", name_base)))]), ]
      
      # Print the number of rows removed
      rows_removed <- nrow(dataframe) - nrow(complete_data)
      print(paste("For model", name_last, ",", rows_removed, "rows were removed due to missing values."))
      
      tryCatch({
        initial_model <- zeroinfl(data = complete_data, 
                                  formula = as.formula(paste(name_last, formula_base, "+", name_base)))
        
        # Use stepAIC to get the optimal model based on AIC but ensure that certain predictors always remain
        scope_list <- list(lower = as.formula(paste("~ sex + entryage + income + prace +", name_base)),
                           upper = as.formula(paste(name_last, formula_base, "+", name_base)))
        
        optimal_model <- MASS::stepAIC(initial_model, direction = "both", scope = scope_list, trace = FALSE)  # trace=FALSE suppresses output
        
        # Extract the right-hand side of the model's formula
        selected_vars <- all.vars(optimal_model$formula)[-1]  # Exclude the outcome variable (it's the 1st item)
        optimal_models[[name_last]] <- selected_vars
        
      }, error = function(e) {
        print(paste("Error encountered for", name_last, ":", e$message))
      })
      
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(optimal_models)
}








# Function to run initial linear models for each outcome
run_initial_models <- function(dataframe, outcome_list, last_time, first_time) {
  model_list <- list()
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Check the proportion of NA values
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {
      if (all(dataframe[[name_last]] == 0)) {
        print(paste(name_last, "contains all 0's. Skipping."))
        next
      }
      
      # Fit initial model
      model <- lm(data = dataframe, 
                  formula = as.formula(paste(name_last, "~ poa_chang + poa_belsky + kdm_advance_0_kdm + sex + entryage + income + prace +", name_base)),
                  na.action = na.exclude)
      model_list[[name_last]] <- model
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(model_list)
}

run_stepwise_binlog_model <- function(dataframe, outcome_list, last_time, first_time) {
  library(MASS)
  
  optimal_models <- list()
  
  formula_base <- " ~ poa_chang*sex + poa_belsky*sex + kdm_advance_0_kdm*sex + entryage + income + prace"
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Check the proportion of NA values in the outcome variable
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
      
      # Ensure the outcome variable is numeric
      dataframe[[name_last]] <- as.numeric(as.character(dataframe[[name_last]]))
      
      tryCatch({
        
        # Drop NA values for the entire dataset based on this model
        complete_data <- dataframe[complete.cases(dataframe[, all.vars(as.formula(paste(name_last, formula_base, "+", name_base)))]), ]
        
        # Print the number of rows removed
        rows_removed <- nrow(dataframe) - nrow(complete_data)
        print(paste("For model", name_last, ",", rows_removed, "rows were removed due to missing values."))
        
        # Check if all values in the column are 0
        if (all(dataframe[[name_last]] == 0)) {
          stop(paste(name_last, "contains all 0's."))
        }
        
        # Oversample the minority class
        oversampled_data <- oversample_minority(complete_data, name_last)
        
        initial_model <- glm(data = oversampled_data, 
                             formula = as.formula(paste(name_last, formula_base, "+", name_base)),
                             family = "binomial")
        
        # Use stepAIC to get the optimal model based on AIC but ensure that certain predictors always remain
        scope_list <- list(lower = as.formula(paste("~ sex + entryage + income + prace +", name_base)),
                           upper = as.formula(paste(name_last, formula_base, "+", name_base)))
        
        optimal_model <- MASS::stepAIC(initial_model, direction = "both", scope = scope_list, trace = FALSE)  # trace=FALSE suppresses output
        
        # Extract the right-hand side of the model's formula
        selected_vars <- optimal_model$coefficients
        optimal_models[[name_last]] <- selected_vars
        
      }, 
      error = function(e) {
        print(paste("Error encountered with", name_last, ":", e$message))
      })
      
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(optimal_models)
}



run_initial_binlog_models <- function(dataframe, outcome_list, last_time, first_time) {
  model_list <- list()
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Check the proportion of NA values
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {
      
      # Ensure the outcome variable is numeric
      dataframe[[name_last]] <- as.numeric(as.character(dataframe[[name_last]]))
      
      tryCatch({
        # Drop NA values for the entire dataset based on this model
        complete_data <- dataframe[complete.cases(dataframe[, all.vars(as.formula(paste(name_last, "~ poa_chang + poa_belsky + kdm_advance_0_kdm + sex + entryage + income + prace +", name_base)))]), ]
        
        # Check if all values in the column are 0
        if (all(dataframe[[name_last]] == 0)) {
          stop(paste(name_last, "contains all 0's."))
        }
        
        # Oversample the minority class
        oversampled_data <- oversample_minority(complete_data, name_last)
        
        # Fit initial binomial logistic regression model
        model <- glm(data = oversampled_data, 
                     formula = as.formula(paste(name_last, "~ poa_chang + poa_belsky + kdm_advance_0_kdm + sex + entryage + income + prace +", name_base)),
                     family = "binomial", 
                     na.action = na.exclude)
        
        model_list[[name_last]] <- model
      },
      error = function(e) {
        print(paste("Error encountered with", name_last, ":", e$message))
      })
      
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(model_list)
}





# Runs models but checks the VIF of each predictor. Only returns VIF
run_initial_models_with_vif <- function(dataframe, outcome_list, last_time, first_time) {
  library(car) # for vif function
  
  results_list <- list()
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Check the proportion of NA values
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {
      if (all(dataframe[[name_last]] == 0)) {
        print(paste(name_last, "contains all 0's. Skipping."))
        next
      }
      
      # Fit initial model
      model <- lm(data = dataframe, 
                  formula = as.formula(paste(name_last, "~ poa_chang + poa_belsky + kdm_advance_0_kdm + sex + entryage + income + prace +", name_base)),
                  na.action = na.exclude)
      
      # Compute VIF for the model
      model_vif <- vif(model)
      
      results_list[[name_last]] <- list(vif = model_vif)
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(results_list)
}


run_initial_binlog_models_with_vif <- function(dataframe, outcome_list, last_time, first_time) { 
  
  library(car)  # for vif function
  
  results_list <- list()
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Check the proportion of NA values
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {
      
      # Ensure the outcome variable is numeric
      dataframe[[name_last]] <- as.numeric(as.character(dataframe[[name_last]]))
      
      # Oversample the minority class for this outcome variable
      dataframe_oversampled <- oversample_minority(dataframe, name_last)
      
      tryCatch({
        # Fit initial binomial logistic regression model
        model <- glm(data = dataframe_oversampled, 
                     formula = as.formula(paste(name_last, "~ poa_chang + poa_belsky + kdm_advance_0_kdm + sex + entryage + income + prace +", name_base)),
                     family = binomial, 
                     na.action = na.exclude)
        
        model_vif <- vif(model)
        
        results_list[[name_last]] <- list(vif = model_vif)
      },
      error = function(e) {
        if(all(dataframe[[name_last]] == 0)) {
          print(paste(name_last, "contains all 0's. Skipping."))
        } else {
          print(paste("Error encountered with", name_last, ":", e$message))
        }
      })
      
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(results_list)
}


run_initial_models_with_vif_zi <- function(dataframe, outcome_list, last_time, first_time) {
  library(car)     # for vif function
  library(pscl)    # for zeroinfl function
  
  results_list <- list()
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Subtract 1 from the outcome variable
    dataframe[[name_last]] <- dataframe[[name_last]] - min(dataframe[[name_last]], na.rm = T)
    # Subtract 1 from the baseline variable as well
    dataframe[[name_base]] <- dataframe[[name_base]] - min(dataframe[[name_base]], na.rm = T)
    
    # Check the proportion of NA values
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {
      if (all(dataframe[[name_last]] == 0)) {
        print(paste(name_last, "contains all 0's. Skipping."))
        next
      }
      
      # Check if there are any usable rows for modeling after handling NAs
      usable_rows <- complete.cases(dataframe[, c("poa_chang", "poa_belsky", "kdm_advance_0_kdm", "sex", "entryage", "income", "prace", name_base, name_last)])
      if (sum(usable_rows) == 0) {
        print(paste("No usable rows for", name_last, ". Skipping."))
        next
      }
      
      tryCatch({
        # Fit zero-inflated model
        zi_model <- zeroinfl(formula = as.formula(paste(name_last, "~ poa_chang + poa_belsky + kdm_advance_0_kdm + sex + entryage + income + prace +", name_base, 
                                                        "| poa_chang + poa_belsky + kdm_advance_0_kdm + sex + entryage + income + prace")), 
                             data = dataframe, dist = "poisson", link = "logit")
        
        # Compute VIF for the count model
        count_model <- lm(formula = as.formula(paste(name_last, "~ poa_chang + poa_belsky + kdm_advance_0_kdm + sex + entryage + income + prace +", name_base)), 
                          data = dataframe)
        count_vif <- vif(count_model)
        
        # Compute VIF for the zero-inflation model (assuming the same predictors as count model for simplicity)
        zi_model_logistic <- glm((dataframe[[name_last]] == 0) ~ poa_chang + poa_belsky + kdm_advance_0_kdm + sex + entryage + income + prace, family = "binomial", data = dataframe)
        zi_vif <- vif(zi_model_logistic)
        
        results_list[[name_last]] <- list(count_vif = count_vif, zi_vif = zi_vif)
      }, error = function(e) {
        print(paste("Error encountered for", name_last, ":", e$message))
      })
      
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(results_list)
}






run_initial_models_zip <- function(dataframe, outcome_list, last_time, first_time) {
  library(pscl)
  model_list <- list()
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Adjust the outcome so its minimum is 0
    dataframe[[name_last]] <- dataframe[[name_last]] - min(dataframe[[name_last]], na.rm = TRUE)
    
    # Check the proportion of NA values
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {
      if (all(dataframe[[name_last]] == 0)) {
        print(paste(name_last, "contains all 0's. Skipping."))
        next
      }
      
      tryCatch({
        # Fit initial zero-inflated model (with both count and zero-inflation parts)
        model_formula <- as.formula(paste(name_last, "~ poa_chang + poa_belsky + kdm_advance_0_kdm + sex + entryage + income + prace +", name_base, 
                                          "|", 
                                          "poa_chang + poa_belsky + kdm_advance_0_kdm + sex + entryage + income + prace +", name_base))
        
        model <- zeroinfl(formula = model_formula, data = dataframe, 
                          dist = "poisson", link = "logit", 
                          na.action = na.exclude)
        model_list[[name_last]] <- model
        
      }, error = function(e) {
        print(paste("Error encountered for", name_last, ":", e$message))
      })
      
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(model_list)
}




# Function to optimize a given list of models by ANOVA comparison
optimize_models_by_anova <- function(model_list) {
  optimized_list <- list()
  
  for (name in names(model_list)) {
    model <- model_list[[name]]
    
    # Extract outcome variable
    outcome_var <- all.vars(formula(model))[1]  # The first variable is the outcome
    
    # DEBUGGING: Check that the outcome variable exists in the dataset
    if (!(outcome_var %in% names(model$model))) {
      stop(paste("Error: Outcome variable", outcome_var, "not found in the dataset for model", name))
    }
    
    # Variables to consider for removal
    potential_remove_vars <- c("poa_chang", "poa_belsky", "kdm_advance_0_kdm")
    current_vars <- setdiff(all.vars(formula(model)), outcome_var)
    
    changed <- TRUE
    
    while(changed) {
      changed <- FALSE
      
      # Try removing variables
      for (var in potential_remove_vars) {
        if (var %in% current_vars) {
          # Construct formula without the variable to be removed
          reduced_formula_str <- paste(outcome_var, "~", paste(setdiff(current_vars, var), collapse = " + "))
          reduced_model <- lm(as.formula(reduced_formula_str), data = model$model)
          
          # Compare models using anova
          comparison <- anova(reduced_model, model)
          p_value <- comparison$"Pr(>F)"[2]
          
          # If p-value is large, then the reduced model is not significantly worse
          if (p_value > 0.05) {
            model <- reduced_model
            current_vars <- setdiff(all.vars(formula(model)), outcome_var)
            changed <- TRUE
          }
        }
      }
    }
    
    optimized_list[[name]] <- model
  }
  
  return(optimized_list)
}


optimize_models_by_chisq <- function(model_list) {
  library(pscl)
  library(lmtest)
  
  optimized_list <- list()
  
  for (name in names(model_list)) {
    
    tryCatch({
      model <- model_list[[name]]
      
      # Check if model$model has less than 50 rows
      if (nrow(model$model) < 50) {
        stop(paste(name, "has less than 50 rows."))
      }
      
      print(model)
      
      # Extract outcome variable
      outcome_var <- all.vars(formula(model))[1]
      
      # DEBUGGING: Check that the outcome variable exists in the dataset
      if (!(outcome_var %in% names(model$model))) {
        stop(paste("Error: Outcome variable", outcome_var, "not found in the dataset for model", name))
      }
      
      # Variables to consider for removal
      potential_remove_vars <- c("poa_chang", "poa_belsky", "kdm_advance_0_kdm")
      current_vars <- setdiff(all.vars(formula(model)), outcome_var)
      
      changed <- TRUE
      
      while(changed) {
        changed <- FALSE
        
        # Try removing variables
        for (var in potential_remove_vars) {
          if (var %in% current_vars) {
            # Construct formula without the variable to be removed
            reduced_formula_str <- paste(outcome_var, "~", paste(setdiff(current_vars, var), collapse = " + "))
            reduced_model <- glm(as.formula(reduced_formula_str), data = model$model, family = "binomial")
            
            # Compare models using anova for GLMs (Likelihood ratio test)
            comparison <- anova(reduced_model, model, test = "Chisq")
            p_value <- comparison$`Pr(>Chi)`[2]
            
            # If p-value is large, then the reduced model is not significantly worse
            if (p_value > 0.05) {
              model <- reduced_model
              current_vars <- setdiff(all.vars(formula(model)), outcome_var)
              changed <- TRUE
            }
          }
        }
      }
      
      optimized_list[[name]] <- model
    },
    error = function(e) {
      print(e$message)
    })
    
  }
  
  return(optimized_list)
}






optimize_models_by_anova_zip <- function(initial_models) {
  
  library(pscl)
  library(lmtest)
  
  optimized_list <- list()
  
  for (name in names(initial_models)) {
    tryCatch({
      model <- initial_models[[name]]
      
      # Extract outcome variable
      outcome_var <- all.vars(formula(model))[1]
      
      # DEBUGGING: Check that the outcome variable exists in the dataset
      if (!(outcome_var %in% names(model$model))) {
        stop(paste("Error: Outcome variable", outcome_var, "not found in the dataset for model", name))
      }
      
      # Variables to consider for removal
      potential_remove_vars <- c("poa_chang", "poa_belsky", "kdm_advance_0_kdm")
      current_vars <- setdiff(all.vars(formula(model)), outcome_var)
      
      changed <- TRUE
      
      while(changed) {
        changed <- FALSE
        
        # Try removing variables
        for (var in potential_remove_vars) {
          if (var %in% current_vars) {
            
            # Construct formula without the variable to be removed
            reduced_formula_str <- paste(outcome_var, "~", paste(setdiff(current_vars, var),
                                                                 collapse = " + "))
            reduced_model <- zeroinfl(as.formula(reduced_formula_str), data = model$model, 
                                      dist = "poisson", link = "logit")
            
            # Compare models using likelihood ratio test
            comparison <- lrtest(reduced_model, model)
            p_value <- comparison$"Pr(>Chisq)"[2]
            
            # If p-value is large, then the reduced model is not significantly worse
            if (p_value > 0.05) {
              model <- reduced_model
              current_vars <- setdiff(all.vars(formula(model)), outcome_var)
              changed <- TRUE
            }
          }
        }
      }
      
      optimized_list[[name]] <- model
      
    }, error = function(e) {
      print(paste("Error encountered for model", name, ":", e$message))
    })
  }
  
  return(optimized_list)
}


# Oversample the minority class in a dataframe
oversample_minority <- function(data, outcome_var) {
  
  # Remove rows where outcome_var is NA
  data <- data[!is.na(data[[outcome_var]]), ]
  
  table_counts <- table(data[[outcome_var]])
  
  # Determine majority and minority classes
  if (table_counts["0"] > table_counts["1"]) {
    majority_class <- 0
    minority_class <- 1
  } else {
    majority_class <- 1
    minority_class <- 0
  }
  
  # Split the data into majority and minority
  majority_data <- data[data[[outcome_var]] == majority_class, ]
  minority_data <- data[data[[outcome_var]] == minority_class, ]
  
  # Sample from the minority class
  oversampled_minority <- minority_data[sample(1:nrow(minority_data), 
                                               nrow(majority_data), 
                                               replace = TRUE), ]
  
  # Combine the majority class with the oversampled minority class
  new_data <- rbind(majority_data, oversampled_minority)
  
  return(new_data)
}


oversample_with_ovun <- function(data, outcome_var) {
  library(ROSE)
  # Remove rows where outcome_var is NA
  data <- data[!is.na(data[[outcome_var]]), ]
  
  # Use ovun.sample to adjust the proportion of the minority class
  oversampled_data <- ovun.sample(as.formula(paste(outcome_var, "~ .")), 
                                  data = data, 
                                  method = "over", 
                                  p = 0.25, 
                                  seed = 12345)$data  # optional: set seed for reproducibility
  
  return(oversampled_data)
}


run_logit_models <- function(dataframe, outcome_list, formula_name, last_time, first_time) {
  library(glmnet)
  
  model_list <- list()
  
  for (name in outcome_list) {
    name_last <- paste0(name, "_", last_time)
    name_base <- paste0(name, "_", first_time)
    
    # Check the proportion of NA values in the outcome variable
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
      
      # Ensure the outcome variable is numeric
      dataframe[[name_last]] <- as.numeric(as.character(dataframe[[name_last]]))
      
      # Oversample the minority class for this outcome variable
      dataframe_oversampled <- oversample_minority(dataframe, name_last)
      
      tryCatch({
        # Fit the logistic regression model on the oversampled data
        mod <- glm(data = dataframe_oversampled, 
                   formula = as.formula(paste(name_last, "~", formula_name, "+", name_base)),
                   family = "binomial", 
                   na.action = na.exclude)  # This will exclude NA values
        
        model_list[[name_last]] <- mod
      },
      error = function(e) {
        if(all(dataframe[[name_last]] == 0)) {
          print(paste(name_last, "contains all 0's. Skipping."))
        } else {
          print(paste("Error encountered with", name_last, ":", e$message))
        }
      })
      
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(model_list)
}




# Extract betas and p-values from logistics regressions
get_betas_and_pvalues_logit <- function(model_list) {
  results_list <- list()
  
  for (name_last in names(model_list)) {
    mod <- model_list[[name_last]]
    
    # Extract the coefficient estimates and their p-values
    betas <- summary(mod)$coefficients[, "Estimate"]
    p_values <- summary(mod)$coefficients[, "Pr(>|z|)"]
    
    # Exclude the intercept (if you want to keep it, remove these lines)
    betas <- betas[-1]
    p_values <- p_values[-1]
    
    # Store results in a dataframe and then in the list
    model_results <- data.frame(Betas = betas, P_Values = p_values)
    results_list[[name_last]] <- model_results
  }
  
  return(results_list)
}

# Gets the accuracy, sensitivity, and specificy from logit models
evaluate_logit_models <- function(dataframe, model_list) {
  library(caret)  # for confusionMatrix function
  results <- list()
  
  for (name_last in names(model_list)) { 
    tryCatch({
      # Check if the column exists in dataframe
      if (name_last %in% colnames(dataframe)) {
        # Predict the outcome variable using the model
        predicted_probs <- predict(model_list[[name_last]], newdata = dataframe, type = "response")
        predicted <- ifelse(predicted_probs >= 0.5, 1, 0)
        
        conf_matrix <- confusionMatrix(as.factor(dataframe[[name_last]]), as.factor(predicted))
        
        # Create a data frame with accuracy, sensitivity, and specificity
        fit_summary <- data.frame(Accuracy = as.numeric(conf_matrix$overall["Accuracy"]), 
                                  Sensitivity = as.numeric(conf_matrix$byClass["Sensitivity"]), 
                                  Specificity = as.numeric(conf_matrix$byClass["Specificity"]))
        
        results[[name_last]] <- fit_summary
      }
    }, error = function(e) {
      print(paste("Error in evaluating", name_last))
      print(e)
    })
  }
  
  return(results)
}




# fit_binomial_logistic <- function(dataframe, outcome_list, formula_name, last_time, first_time) {
#   results <- list()
#   
#   for (name in outcome_list) { 
#     tryCatch({
#       name_last <- paste0(name, "_", last_time)
#       
#       # Check the proportion of NA values in the outcome variable
#       na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
#       
#       if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
#         # Check if all values in the column are 0
#         if (all(dataframe[[name_last]] == 0)) {
#           print(paste(name_last, "contains all 0's. Skipping."))
#           next  # Go to the next iteration
#         }
#         
#         name_base <- paste0(name, "_", first_time)
#         
#         # Create a binomial logistic regression model
#         mod <- glm(as.formula(paste(name_last, "~", formula_name, "+", name_base)), 
#                    data = dataframe, family = binomial)
#         
#         # Predict the outcome variable using the model
#         predicted_probs <- predict(mod, newdata = dataframe, type = "response")
#         predicted <- ifelse(predicted_probs >= 0.5, 1, 0)
#         
#         conf_matrix <- confusionMatrix(as.factor(dataframe[[name_last]]), as.factor(predicted))
#         
#         # Create a data frame with accuracy, sensitivity, and specificity
#         fit_summary <- data.frame(Accuracy = as.numeric(conf_matrix$overall["Accuracy"]), 
#                                   Sensitivity = as.numeric(conf_matrix$byClass["Sensitivity"]), 
#                                   Specificity = as.numeric(conf_matrix$byClass["Specificity"]))
#         
#         results[[name_last]] <- fit_summary
#       } else {
#         print(paste("Skipping", name_last, "due to more than 75% NA values"))
#       }
#     }, error = function(e) {
#       print(paste("Error in modeling", name_last))
#       print(e)
#     })
#   }
#   
#   return(results)
# }







fit_multinomial_logistic <- function(dataframe, outcome_list, formula_name, last_time, first_time) {
  
  library(nnet)
  results <- list()
  
  for (name in outcome_list) { 
    tryCatch({
      name_last <- paste0(name, "_", last_time)
      
      # Check the proportion of NA values in the outcome variable
      na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
      
      if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
        # Check if all values in the column are 0
        if (all(dataframe[[name_last]] == 0)) {
          print(paste(name_last, "contains all 0's. Skipping."))
          next  # Go to the next iteration
        }
        
        name_base <- paste0(name, "_", first_time)
        
        # Create a multinomial logistic regression model
        mod <- multinom(as.formula(paste(name_last, "~", formula_name, "+", name_base)), data = dataframe)
        
        # Predict the outcome variable using the model
        predicted <- predict(mod, newdata = dataframe, type = "class")
        
        conf_matrix <- table(as.factor(dataframe[[name_last]]), predicted)
        
        # Calculate True Positives (TP), True Negatives (TN),
        # False Positives (FP), and False Negatives (FN)
        TP <- diag(conf_matrix)
        TN <- sum(conf_matrix) - rowSums(conf_matrix) - colSums(conf_matrix) + TP
        FP <- colSums(conf_matrix) - TP
        FN <- rowSums(conf_matrix) - TP
        
        # Calculate accuracy
        accuracy <- (TP + TN) / sum(conf_matrix)
        # Calculate sensitivity (True Positive Rate)
        sensitivity <- TP / rowSums(conf_matrix)
        # Calculate specificity
        specificity <- TN / (TN + FP)
        
        # Create a data frame with accuracy, sensitivity, and specificity
        fit_summary <- data.frame(Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity)
        
        results[[name_last]] <- fit_summary
      } else {
        print(paste("Skipping", name_last, "due to more than 75% NA values"))
      }
    }, error = function(e) {
      print(paste("Error in modeling", name_last))
      print(e)
    })
  }
  
  return(results)
}

# Function 1: Fit multinomial logistic models
fit_multinomial_logistic_models <- function(dataframe, outcome_list, formula_name, last_time, first_time) {
  
  library(nnet)
  model_results <- list()
  
  for (name in outcome_list) {
    tryCatch({
      name_last <- paste0(name, "_", last_time)
      
      # Check the proportion of NA values in the outcome variable
      na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
      
      if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
        # Check if all values in the column are 0
        if (all(dataframe[[name_last]] == 0)) {
          print(paste(name_last, "contains all 0's. Skipping."))
          next  # Go to the next iteration
        }
        
        name_base <- paste0(name, "_", first_time)
        
        # Create a multinomial logistic regression model
        mod <- multinom(as.formula(paste(name_last, "~", formula_name, "+", name_base)), data = dataframe)
        
        model_results[[name_last]] <- mod
        
      } else {
        print(paste("Skipping", name_last, "due to more than 75% NA values"))
      }
    }, error = function(e) {
      print(paste("Error in modeling", name_last))
      print(e)
    })
  }
  
  return(model_results)
}

# Function 2: Evaluate multinomial logistic models
evaluate_multinomial_models <- function(dataframe, model_list) {
  eval_results <- list()
  
  for (name_last in names(model_list)) {
    mod <- model_list[[name_last]]
    
    # Predict the outcome variable using the model
    predicted <- predict(mod, newdata = dataframe, type = "class")
    
    conf_matrix <- table(as.factor(dataframe[[name_last]]), predicted)
    
    # Calculate True Positives (TP), True Negatives (TN),
    # False Positives (FP), and False Negatives (FN)
    TP <- diag(conf_matrix)
    TN <- sum(conf_matrix) - rowSums(conf_matrix) - colSums(conf_matrix) + TP
    FP <- colSums(conf_matrix) - TP
    FN <- rowSums(conf_matrix) - TP
    
    # Calculate accuracy
    accuracy <- (TP + TN) / sum(conf_matrix)
    # Calculate sensitivity (True Positive Rate)
    sensitivity <- TP / rowSums(conf_matrix)
    # Calculate specificity
    specificity <- TN / (TN + FP)
    
    # Create a data frame with accuracy, sensitivity, and specificity
    eval_summary <- data.frame(Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity)
    
    eval_results[[name_last]] <- eval_summary
  }
  
  return(eval_results)
}


# This function gets the significant predictors of the binomial multinomial function
get_significant_predictors_multinom <- function(model_list, threshold = 0.05) {
  results_list <- list()
  
  pValue_extract <- function(x){
    z <- summary(x)$coefficients/summary(x)$standard.errors
    # 2-tailed Wald z tests to test significance of coefficients
    p <- (1 - pnorm(abs(z), 0, 1)) * 2
    p
  }
  
  for (name_last in names(model_list)) {
    mod <- model_list[[name_last]]
    
    # Extracting coefficients for each level of the outcome
    coefficients_list <- summary(mod)$coefficients
    names(coefficients_list)
    
    significance_list <- list()
    for (level in names(coefficients_list)) {
      betas <- coefficients_list
      # Use the custom function for p-values
      p_values <- pValue_extract(mod)
      
      # Exclude the intercept and focus only on coefficients containing "kdm" or "poa"
      mask <- grepl("kdm|poa", names(betas))
      betas <- betas[mask]
      p_values <- p_values[mask]
      
      # Filter for significance
      significant_vars <- names(betas)[p_values < 0.05]
      
      if (length(significant_vars) > 0) {
        significance_list[[level]] <- data.frame(Variable = significant_vars, Betas = betas[significant_vars], P_Values = p_values[significant_vars])
      }
    }
    
    if (length(significance_list) > 0) {
      results_list[[name_last]] <- significance_list
    }
  }
  
  return(results_list)
}



# Description:
#   This R function fits a Zero-Inflated Poisson (ZIP) regression model to a list of specified 
# outcome variables in a given dataset. It calculates model performance metrics including Mean 
# Absolute Error (MAE), Root Mean Squared Error (RMSE), and the Vuong test p-value for zero-inflation. 
# The function is designed to handle datasets with potential zero-inflation issues and provides 
# insights into the quality of the ZIP regression models for each specified outcome.

# Parameters:
#   
# dataframe: The dataset containing the outcome variables and predictors.
# outcome_list: A list of outcome variable names to be analyzed.
# formula_name: The formula specifying the model structure, including predictors.
# last_time: A character string specifying the time point for the outcome variable.
# first_time: A character string specifying the initial time point for predictor variables.
# 
# Output:
#   A list of performance metrics for each specified outcome variable, including MAE, RMSE, and the 
# Vuong test p-value for zero-inflation. The results help evaluate the quality of the ZIP 
# regression models.

fit_zip_getOutcomeMetrics <- function(dataframe, outcome_list, formula_name, last_time, first_time) {
  results <- list()
  
  for (name in outcome_list) { 
    tryCatch({
      name_last <- paste0(name, "_", last_time)
      
      # Check the proportion of NA values in the outcome variable
      na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
      
      if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
        # Check if all values in the column are 0
        if (all(dataframe[[name_last]] == 0)) {
          print(paste(name_last, "contains all 0's. Skipping."))
          next  # Go to the next iteration
        }
        
        name_base <- paste0(name, "_", first_time)
        
        print(paste(name_last, "~", formula_name, "+", name_base))
        
        # Fit a Zero-Inflated Poisson Regression model
        mod <- zeroinfl(formula = as.formula(paste(name_last, "~", formula_name, "+", name_base)), 
                        data = dataframe, na.action = na.exclude)  # This will exclude NA values
        
        # Predict the outcome variable using the model
        predicted <- predict(mod, newdata = dataframe, type = "response")
        
        # Calculate MAE and RMSE
        mae_value <- mean(abs(dataframe[[name_last]] - predicted), na.rm = TRUE)
        rmse_value <- sqrt(mean((dataframe[[name_last]] - predicted)^2, na.rm = TRUE))
        
        r_squared <- as.numeric(r2_zeroinflated(mod)$R2_adjusted)
        
        # Create a data frame with the metrics
        fit_summary <- data.frame(MAE = mae_value, RMSE = rmse_value, R_Squared = r_squared)
        
        results[[name_last]] <- fit_summary
      } else {
        print(paste("Skipping", name_last, "due to more than 75% NA values"))
      }
    }, error = function(e) {
      print(paste("Error in modeling", name_last))
      print(e)
    })
  }
  
  return(results)
}






# 
# The compare_lm_metrics function is designed to help you compare the R-squared (R²), 
# Mean Absolute Error (MAE), and Root Mean Squared Error (RMSE) metrics from the results 
# obtained by running the fit_lm_getOutcomeMetrics function on multiple dataframes. 
# Here's a summary of what this function does:
# 
# Summary:
# Input: It takes a list of result data frames as its input, where each data frame contains R², 
# MAE, and RMSE values for different outcomes obtained from running fit_lm_getOutcomeMetrics on 
# various dataframes.
# Output: The function returns a list of combined data frames for each metric (R², MAE, RMSE), 
# making it easier to compare these metrics across different dataframes.

compare_lm_metrics <- function(chang, belsky, kdm) {
  # Initialize data frames to store the combined metrics
  combined_r_squared <- data.frame()
  combined_mae <- data.frame()
  combined_rmse <- data.frame()
  
  # Iterate through the list of results
  for (i in seq_along(result_list)) {
    # Extract R-squared, MAE, and RMSE from the results
    r_squared <- result_list[[i]]$R_Squared
    mae <- result_list[[i]]$MAE1
    rmse <- result_list[[i]]$RMSE1
    
    # Create a data frame with the extracted metrics
    metrics_df <- data.frame(Dataframe = i, R_Squared = r_squared, MAE = mae, RMSE = rmse)
    
    # Append to the combined data frames
    combined_r_squared <- rbind(combined_r_squared, metrics_df)
    combined_mae <- rbind(combined_mae, metrics_df)
    combined_rmse <- rbind(combined_rmse, metrics_df)
  }
  
  # Return the combined metrics data frames
  return(list(R_Squared = combined_r_squared, MAE = combined_mae, RMSE = combined_rmse))
}



run_zero_inflated_model <- function(dataframe, outcome_list, formula_name, last_time, first_time) {
  library(pscl)
  
  models <- list()
  
  for (name in outcome_list) { 
    skip_iteration <- FALSE
    name_last <- paste0(name, "_", last_time)
    
    # Subtract 1 from the outcome variable
    dataframe[[name_last]] <- dataframe[[name_last]] - min(dataframe[[name_last]], na.rm = T)
    
    # Check the proportion of NA values in the outcome variable
    na_proportion <- sum(is.na(dataframe[[name_last]])) / length(dataframe[[name_last]])
    
    if (na_proportion <= 0.75) {  # Skip if more than 75% NAs
      # Check if all values in the column are 0
      if (all(dataframe[[name_last]] == 0)) {
        print(paste(name_last, "contains all 0's. Skipping."))
        next  # Go to the next iteration
      }
      
      name_base <- paste0(name, "_", first_time)
      
      # Subtract 1 from the baseline variable as well
      dataframe[[name_base]] <- dataframe[[name_base]] - min(dataframe[[name_base]], na.rm = T)
      
      # Check if more than 96% of the values in dataframe[[name_base]] are zeros
      if (mean(dataframe[[name_base]] == 0, na.rm = TRUE) > 0.96) {
        print(paste(name_base, "contains more than 96% zeros. Skipping."))
        next  # Go to the next iteration
      }
      
      tryCatch({
        # Fit a Zero-Inflated Poisson Regression model
        mod <- zeroinfl(formula = as.formula(paste(name_last, "~", formula_name, "+", name_base)), 
                        data = dataframe, na.action = na.exclude)
        
        models[[name_last]] <- mod
        
      }, warning = function(w) {
        # Handle the singularity warning
        if (grepl("system is computationally singular", w$message)) {
          print(paste("Model for", name_last, "is singular. Skipping this iteration."))
          skip_iteration <- TRUE
        } else {
          # If not the singularity warning, re-raise the warning
          warning(w)
        }
      }, error = function(e) {
        print(paste("Error in modeling", name_last))
        print(e)
      })
      
      if (skip_iteration) {
        next
      }
      
    } else {
      print(paste("Skipping", name_last, "due to more than 75% NA values"))
    }
  }
  
  return(models)
}








extract_coefficients_zip <- function(model_list, threshold = 0.05) {
  results <- list()
  
  for (name in names(model_list)) {
    model <- model_list[[name]]
    
    # Extracting coefficients and p-values for both count and zero-inflation components
    count_coefs <- coef(summary(model))$count[, "Estimate"]
    count_pvalues <- coef(summary(model))$count[, "Pr(>|z|)"]
    
    zero_coefs <- coef(summary(model))$zero[, "Estimate"]
    zero_pvalues <- coef(summary(model))$zero[, "Pr(>|z|)"]
    
    # Replace NaN values in zero_pvalues with 1.0
    count_pvalues[is.nan(count_pvalues)] <- 1.0
    zero_pvalues[is.nan(zero_pvalues)] <- 1.0
    
    # Filtering for significance and variable names containing "poa" or "kdm"
    sig_count_idx <- (count_pvalues < threshold) & (grepl("poa|kdm", names(count_pvalues)))
    sig_zero_idx <- (zero_pvalues < threshold) & (grepl("poa|kdm", names(zero_pvalues)))
    
    significant_df <- data.frame(
      Variable = c(names(count_coefs)[sig_count_idx], names(zero_coefs)[sig_zero_idx]),
      Betas_Count = c(count_coefs[sig_count_idx], rep(NA, sum(sig_zero_idx))),
      P_Values_Count = c(count_pvalues[sig_count_idx], rep(NA, sum(sig_zero_idx))),
      Betas_Zero = c(rep(NA, sum(sig_count_idx)), zero_coefs[sig_zero_idx]),
      P_Values_Zero = c(rep(NA, sum(sig_count_idx)), zero_pvalues[sig_zero_idx])
    )
    
    # Only add to results if there are any significant coefficients
    if (nrow(significant_df) > 0) {
      results[[name]] <- significant_df
    }
  }
  
  return(results)
}




get_fit_metrics_zip <- function(model_list, dataframe) {
  results <- list()
  
  rmse1 <- function(o, p, m = TRUE) {
    sqrt(mean((o - p)^2, na.rm = m))
  }
  
  mae1 <- function(o, p, m = TRUE) {
    mean(abs(o - p), na.rm = m)
  }
  
  for (name in names(model_list)) {
    model <- model_list[[name]]
    
    # Predict using the count component
    predicted <- predict(model, newdata = dataframe, type = "count")
    true_values <- dataframe[, gsub("formula", "", formula(model)[2])]
    
    # Calculate MAE and RMSE
    mae_value <- mae1(true_values, predicted)
    rmse_value <- rmse1(true_values, predicted)
    
    # Storing results in a list
    results[[name]] <- list(MAE = mae_value, RMSE = rmse_value)
  }
  
  return(results)
}





## Create a list of result data frames from fit_lm_getOutcomeMetrics
# results_list <- list(result_df1, result_df2, result_df3)  # Replace with your actual result data frames

## Compare the metrics
# metrics_comparison <- compare_lm_metrics(results_list)

## Access the combined metrics for R-squared, MAE, and RMSE
# r_squared_combined <- metrics_comparison$R_Squared
# mae_combined <- metrics_comparison$MAE
# rmse_combined <- metrics_comparison$RMSE



# The combine_and_name_dataframes function takes three lists of data frames (hold, hold1, and hold2) 
# and combines them into a list of data frames with specific column names. It adds a new column with 
# predefined values ("Chang", "Belsky", "KDM") to each data frame and assigns names to these data 
# frames based on the names in the hold list. The resulting list contains named data frames ready 
# for further analysis or processing.

combine_and_name_dataframes <- function(hold, hold1, hold2) {
  # Create an empty list to store the data frames
  df1_list <- list()
  
  for (i in seq_along(hold)) {
    df1 <- rbind(rbind(hold[[i]], hold1[[i]]), hold2[[i]])
    
    # Create a vector with the desired values
    names_col <- c("Chang", "Belsky", "KDM")
    
    # Add this vector as the first column to df1
    df1 <- cbind(Names = names_col, df1)
    
    # Get the name of the ith element in hold
    df_name <- names(hold)[i]
    
    # Assign the name to the dataframe
    names(df1) <- c("Names", names(df1)[-1])
    
    # Store df1 in the list with the same name
    df1_list[[df_name]] <- df1
  }
  
  return(df1_list)
}



# Function to get significant p-values
get_significant_vars <- function(data_list, var_string, threshold = 0.05) {
  significant_vars <- list()
  
  # Extract significant p-values
  for (name in names(data_list)) {
    df <- data_list[[name]]
    # Identify rows containing the var_string
    match_rows <- grep(var_string, rownames(df))
    significant <- df[match_rows, ]
    significant <- significant[significant$P_Values <= threshold, ]
    if (nrow(significant) > 0) {
      significant_vars[[name]] <- significant
    }
  }
  
  return(significant_vars)
}






# This R function, compare_multinomial_logistic_metrics, takes three lists of outcome metrics 
# (a, b, c), where a is chang, b is belsky, and c is kdm, obtained from multinomial logistic regression models. 
# It iterates through these lists, extracts accuracy, sensitivity, and specificity metrics for 
# each outcome if they exist, and compiles them into a combined data frame for each of the three 
# datasets: "Chang," "Belsky," and "KDM." It then calculates the mean accuracy, sensitivity, and 
# specificity for each dataset while ignoring NA values and returns a summary data frame indicating
# which dataset performs best in terms of these metrics.

compare_multinomial_logistic_metrics <- function(a, b, c) {
  # Initialize data frames to store the combined metrics
  combined_accuracy <- data.frame()
  combined_sensitivity <- data.frame()
  combined_specificity <- data.frame()
  
  # Iterate through the list of outcomes
  for (i in seq_along(a)) {
    # Initialize variables to store metrics
    accuracy1 <- sensitivity1 <- specificity1 <- NA
    accuracy2 <- sensitivity2 <- specificity2 <- NA
    accuracy3 <- sensitivity3 <- specificity3 <- NA
    
    # Extract metrics from outcomes if they exist
    tryCatch({
      accuracy1 <- a[[i]]$Accuracy
      sensitivity1 <- a[[i]]$Sensitivity
      specificity1 <- a[[i]]$Specificity
    }, error = function(e) {})
    
    tryCatch({
      accuracy2 <- b[[i]]$Accuracy
      sensitivity2 <- b[[i]]$Sensitivity
      specificity2 <- b[[i]]$Specificity
    }, error = function(e) {})
    
    tryCatch({
      accuracy3 <- c[[i]]$Accuracy
      sensitivity3 <- c[[i]]$Sensitivity
      specificity3 <- c[[i]]$Specificity
    }, error = function(e) {})
    
    # Create data frames with the extracted metrics
    metrics_chang <- data.frame(Dataframe = "Chang", Outcome = names(a)[i], Accuracy = accuracy1, 
                                Sensitivity = sensitivity1, Specificity = specificity1)
    metrics_belsky <- data.frame(Dataframe = "Belsky", Outcome = names(b)[i], Accuracy = accuracy2, 
                                 Sensitivity = sensitivity2, Specificity = specificity2)
    metrics_kdm <- data.frame(Dataframe = "KDM", Outcome = names(c)[i], Accuracy = accuracy3, 
                              Sensitivity = sensitivity3, Specificity = specificity3)
    
    # Append to the combined data frames
    combined_accuracy <- rbind(combined_accuracy, metrics_chang, metrics_belsky, metrics_kdm)
    combined_sensitivity <- rbind(combined_sensitivity, metrics_chang, metrics_belsky, metrics_kdm)
    combined_specificity <- rbind(combined_specificity, metrics_chang, metrics_belsky, metrics_kdm)
    
  }
  
  
  # Determine the best dataframe based on accuracy, sensitivity, and specificity
  best_accuracy_df <- combined_accuracy %>%
    group_by(Dataframe) %>%
    filter(all(is.finite(Accuracy))) %>%  # Filter out NaN and Inf
    summarize(Mean_Accuracy = mean(Accuracy, na.rm = TRUE))
  
  best_sensitivity_df <- combined_sensitivity %>%
    group_by(Dataframe) %>%
    filter(all(is.finite(Sensitivity))) %>%  # Filter out NaN and Inf
    summarize(Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE))
  
  best_specificity_df <- combined_specificity %>%
    group_by(Dataframe) %>%
    filter(all(is.finite(Specificity))) %>%  # Filter out NaN and Inf
    summarize(Mean_Specificity = mean(Specificity, na.rm = TRUE))
  
  # Combine the metrics for accuracy, sensitivity, and specificity
  combined_metrics <- merge(best_accuracy_df, best_sensitivity_df, by = "Dataframe")
  combined_metrics <- merge(combined_metrics, best_specificity_df, by = "Dataframe")
  
  
  return(combined_metrics)
}


# The "transform_columns_zero_scale" function takes a dataframe and a list of column names as input. 
# It flips the distribution of specified columns so that the left tail becomes the heavy tail and sets 
# the lowest value to zero while adjusting all other values accordingly. This transformation is applied 
# to each specified column, and the modified dataframe is returned as output.

# Example usage:
# Assuming you have a dataframe df with columns you want to transform: col1, col2, col3
# columns_to_transform <- c("col1", "col2", "col3")
# transformed_df <- transform_columns(df, columns_to_transform)

transform_columns_zero_scale <- function(dataframe, column_names) {
  if (!is.vector(column_names)) {
    column_names <- list(column_names)
  }
  
  for (col in column_names) {
    if (skewness(dataframe[[col]], na.rm = TRUE) < 0) {
      # If left tail heavy, flip the distribution
      min_val <- min(dataframe[[col]], na.rm = TRUE)
      dataframe[[col]] <- max(dataframe[[col]], na.rm = TRUE) - dataframe[[col]]
    } else {
      # If right tail heavy, set the lowest value to zero
      min_val <- min(dataframe[[col]], na.rm = TRUE)
      dataframe[[col]] <- dataframe[[col]] - min_val
    }
  }
  
  return(dataframe)
}


# Function Name: binarize_columns
# 
# Description: This function allows you to binarize selected columns in a DataFrame by specifying a 
# threshold value. Values greater than the threshold are set to 1, while values less than or equal 
# to the threshold are set to 0. It takes the DataFrame, a list of column names, and the threshold as
# input and returns the modified DataFrame.

binarize_columns <- function(dataframe, column_names, threshold) {
  if (!is.vector(column_names)) {
    column_names <- list(column_names)
  }
  
  for (col in column_names) {
    dataframe[[col]] <- ifelse(dataframe[[col]] >= threshold, "1", "0")
    dataframe[[col]] <- as.factor(dataframe[[col]])
  }
  
  return(dataframe)
}


create_combined_plots <- function(dataframe, outcome_names, last_time = "_96") {
  # Ensure outcome_names is a character vector
  if (!is.character(outcome_names)) {
    stop("outcome_names must be a character vector.")
  }
  
  # Add last_time to each outcome name
  outcome_names <- paste0(outcome_names, last_time)
  
  # Create an empty list to store the plots
  plot_lists <- list()
  
  # Iterate through each outcome variable
  for (outcome_name in outcome_names) {
    # Check if the outcome variable is character and convert to factor
    if (is.character(dataframe[[outcome_name]])) {
      dataframe[[outcome_name]] <- as.factor(dataframe[[outcome_name]])
    }
    
    # Check if the outcome variable is categorical (a factor)
    if (is.factor(dataframe[[outcome_name]])) {
      # Create boxplots for the predictor variables
      plots <- lapply(c("poa_chang", "poa_belsky", "kdm_advance_0_kdm"), function(predictor_name) {
        ggplot(data = dataframe, aes_string(x = outcome_name, y = predictor_name)) +
          geom_boxplot() +
          labs(title = paste("Boxplot of", outcome_name),
               x = outcome_name, y = predictor_name) +
          theme_minimal()
      })
    } else {
      # Create scatterplots for the predictor variables
      plots <- lapply(c("poa_chang", "poa_belsky", "kdm_advance_0_kdm"), function(predictor_name) {
        ggplot(data = dataframe, aes_string(x = predictor_name, y = outcome_name)) +
          geom_point() + 
          geom_smooth(method = "lm", se = TRUE) +
          labs(title = paste("Scatterplot of", predictor_name),
               x = predictor_name, y = outcome_name) +
          theme_minimal()
      })
    }
    
    # Store the plots in the list
    plot_lists[[outcome_name]] <- plots
  }
  
  # Return the list of plots
  return(plot_lists)
}


print_histograms_and_lowest_proportions <- function(df, column_names) {
  proportions_list <- list()
  
  for (column_name in column_names) {
    # Check if column exists in dataframe
    if (!column_name %in% names(df)) {
      warning(paste("Column", column_name, "does not exist in the dataframe."))
      next
    }
    
    # Draw histogram
    hist(df[[column_name]], main = paste("Histogram of", column_name), xlab = column_name)
    
    # Calculate the proportion of the lowest value
    lowest_val <- min(df[[column_name]], na.rm = TRUE)
    proportion <- sum(df[[column_name]] == lowest_val, na.rm = TRUE) / sum(!is.na(df[[column_name]]))
    proportions_list[[column_name]] <- proportion
    
    # Print the proportion of the lowest value
    cat(paste("Proportion of the lowest value (", lowest_val, ") in column", column_name, ":", proportion, "\n"))
  }
  
  return(proportions_list)
}

# Example usage:
# Assuming 'your_dataframe' is your dataframe and you want to check columns 'Age' and 'Salary'
result <- print_histograms_and_lowest_proportions(your_dataframe, c("Age", "Salary"))









