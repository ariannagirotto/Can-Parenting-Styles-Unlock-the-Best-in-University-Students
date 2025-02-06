########## -------------------------------------------------------- ###########

################################## LIBRARIES ##################################

library(boot)
library(broom)
library(car)  
library(caret)
library(conflicted)
library(ContaminatedMixt)
library(corpcor)
library(dplyr)
library(flexmix)
library(ggplot2)
library(haven)
library(insight)
library(LMest)
library(lme4)
library(longclust)
library(knitr)
library(MASS)
library(mclust)
library(pROC)
library(ranger)
library(reshape2)
library(robust)
library(robustbase)
library(rpart)
library(rpart.plot)
library(stats)
library(tidyverse)
library(umap)
library(tclust)
library(cluster)
library(factoextra)
library(Metrics)
library(fmsb)
library(tibble)
library(tidyr)
library(summarytools)



########## -------------------------------------------------------- ###########

############################### DATA PREPARATION ###############################

data <- read_sav("/Users/ariannagirotto/Desktop/PRIMO TRIMESTRE/PROGETTO ADVANCED/Parenting Styles and Character Strengths.sav")

### Checking for missing values --> no missing values found
sum(is.na(data))

### Checking for duplicates --> no duplicates found
duplicate_counts <- table(apply(data, 1, paste, collapse = ","))
duplicates_with_counts <- subset(duplicate_counts, duplicate_counts > 1)
duplicates_with_counts

### Renaming columns to English 
new_names <- c("Courage", "Creativity", "Curiosity", "Justice", "Forgiveness", "Gratitude", "Honesty",
               "Hope", "Humor", "Judgment", "Kindness", "Leadership", "Love_of_Learning", "Love", 
               "Humility", "Perseverance", "Perspective", "Prudence", "Self_Regulation", "Social_Intelligence",
               "Spirituality", "Teamwork", "Zest")

colnames(data)[(ncol(data) - 22):ncol(data)] <- new_names

new_initial_names <- c("Consent_to_Participate", "Gender", "Age",
                       "Faculty", "Year_of_Study")

colnames(data)[2:6] <- new_initial_names

### Dropping unnecessary columns (1,2 and from 7 to 182)
data <- data[, -c(1, 2, 7:182)]

### Standardization of the variables (from 0 to 1)
data_normalized_robust <- data %>%
  mutate(across(-Gender, ~ (.-min(., na.rm = TRUE)) / 
                  (max(., na.rm = TRUE) - min(., na.rm = TRUE))))


########## -------------------------------------------------------- ###########

############################ DESCRIPTIVE STATISTIC #############################

###### a) Data summary
data_summary <- descr(data, 
                      stats = c("mean", "sd", "min", "max"), 
                      transpose = TRUE, 
                      order = "preserve")

print(data_summary)

###### b) Heatmap

cor_matrix <- cor(data_normalized_robust[, sapply(data_normalized_robust, is.numeric)], 
                  use = "pairwise.complete.obs")

cor_matrix_long <- melt(cor_matrix)

ggplot(data = cor_matrix_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() + 
  scale_fill_viridis_c(direction = -1) +  
  geom_text(aes(label = round(value, 2)), vjust = 1, color = "white", size= 4) +  
  labs(title = "Correlation Heatmap", x = 'Variables', y = 'Variables') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", hjust = 0.5))


########## -------------------------------------------------------- ###########

############################## TRAIN AND TEST SET ##############################

data_1 <- data_normalized_robust
data_normalized_robust

### Splitting with training test at 70 %
set.seed(123)  
train_index <- createDataPartition(data_1$EVMT, p = 0.7, list = FALSE)

# Creating train e test set
train <- data_1[train_index, ]
test <- data_1[-train_index, ]

# Checking the dimensions
cat("Training set dimension:", nrow(train), "\n")
cat("Test set dimension:", nrow(test), "\n")

### Checking the distribution of the response variable
table_train <- table(train$EVMT)
table_test <- table(test$EVMT)

# Proportions
cat("TRAIN distribution:\n")
prop.table(table_train)

cat("\nTEST distribution:\n")
prop.table(table_test)


########## -------------------------------------------------------- ############

################## FORWARD SELECTION FOR REGRESSION ANALYSIS ###################

########## 1.Calculate Average Proportion of Outliers Using Bootstrap ##########

### 1.1 Define Dependent and Independent Variables 

dep_vars <- c("Courage", "Creativity", "Curiosity", "Justice", "Forgiveness", 
              "Gratitude", "Honesty", "Hope", "Humor", "Judgment", "Kindness", 
              "Leadership", "Love_of_Learning", "Love", "Humility", "Perseverance", 
              "Perspective", "Prudence", "Self_Regulation", "Social_Intelligence", 
              "Spirituality", "Teamwork", "Zest")

ind_vars <- c("DAT", "KAT", "OAT", "DBT", "KBT", "OBT")


### 1.2 Define the Function to Calculate Average Proportion of Outliers

calculate_outlier_proportion <- function(train_data, n_bootstrap = 100, alpha = 0.975) {
  if (!all(sapply(train_data, is.numeric))) {
    stop("All columns in the dataset must be numeric.")
  }
  
  # Bootstrap to calculate the proportion of outliers
  bootstrap_proportions <- replicate(n_bootstrap, {
    train_sample <- train_data[sample(1:nrow(train_data), replace = TRUE), ]
    
    fit_MCD <- covMcd(train_sample)
    robust_distances <- mahalanobis(
      x = train_sample,
      center = fit_MCD$raw.center,
      cov = fit_MCD$raw.cov
    )
    
    threshold <- qchisq(p = alpha, df = ncol(train_sample))
    
    outliers <- robust_distances > threshold
    
    mean(outliers)
  })
  
  mean_proportion_outliers <- mean(bootstrap_proportions)
  
  cat(sprintf("Average proportion of outliers: %.2f%%\n", mean_proportion_outliers * 100))
  return(mean_proportion_outliers)
}

### 1.3 Initialize a Results Data Frame

outlier_proportions <- data.frame(
  Dependent_Variable = character(),
  Average_Outlier_Proportion = numeric(),
  stringsAsFactors = FALSE
)

### 1.4 Iterate Over Each Dependent Variable

for(dep_var in dep_vars){
  cat("Processing Dependent Variable:", dep_var, "\n")
  
  subset_data <- train %>%
    dplyr::select(all_of(c(dep_var, ind_vars)))
  
  if(!all(sapply(subset_data, is.numeric))){
    warning(paste("Skipping", dep_var, "- Non-numeric columns detected."))
    next
  }
  
  avg_outlier_prop <- calculate_outlier_proportion(subset_data, n_bootstrap = 100, alpha = 0.975)
  
  outlier_proportions <- outlier_proportions %>%
    dplyr::bind_rows(data.frame(
      Dependent_Variable = dep_var,
      Average_Outlier_Proportion = avg_outlier_prop,
      stringsAsFactors = FALSE
    ))
  
  cat("\n")  
}

### 1.5 Present the Results

outlier_proportions_table <- outlier_proportions %>%
  dplyr::mutate(Average_Outlier_Proportion = sprintf("%.2f%%", Average_Outlier_Proportion * 100))

print(
  outlier_proportions_table %>%
    knitr::kable(format = "markdown", caption = "Average Proportion of Outliers per Dependent Variable")
)


########## 2. Calculate the Robust and Not Robust Linear Regressions ###########

### 2.1 Define Dependent and Independent Variables

dep_vars <- c("Courage", "Creativity", "Curiosity", "Justice", "Forgiveness", 
              "Gratitude", "Honesty", "Hope", "Humor", "Judgment", "Kindness", 
              "Leadership", "Love_of_Learning", "Love", "Humility", "Perseverance", 
              "Perspective", "Prudence", "Self_Regulation", "Social_Intelligence", 
              "Spirituality", "Teamwork", "Zest")

ind_vars <- c("DAT", "KAT", "OAT", "DBT", "KBT", "OBT")


### 2.3 Subset Training and Test Data

train_sub <- train[, c(dep_vars, ind_vars)]
test_sub <- test[, c(dep_vars, ind_vars)]

### 2.4 Function for Significance Codes
significance_code <- function(pval) {
  if (is.na(pval)) return("")
  if (pval < 0.001) return("***")
  else if (pval < 0.01) return("**")
  else if (pval < 0.05) return("*")
  else if (pval < 0.1) return(".")
  else return("")
}

### 2.5 Variable Selection with Forward Selection and AIC

select_variables <- function(train_data, response_var, predictor_vars, robust = FALSE) {
  null_model <- lm(as.formula(paste(response_var, "~ 1")), data = train_data)
  full_model <- lm(as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + "))), data = train_data)
  
  forward_model <- stepAIC(
    null_model,
    scope = list(lower = null_model, upper = full_model),
    direction = "forward",
    trace = FALSE
  )
  
  if (robust) {
    selected_predictors <- names(coef(forward_model))[-1]  
    if (length(selected_predictors) == 0) {
      return(list(forward_model = forward_model, lts_model = NULL))
    }
    robust_formula <- as.formula(paste(response_var, "~", paste(selected_predictors, collapse = " + ")))
    lts_model <- ltsReg(robust_formula, data = train_data)
    return(list(forward_model = forward_model, lts_model = lts_model))
  }
  
  return(forward_model)
}

### 2.6 Bootstrap Validation for Variable Selection
bootstrap_validation <- function(train_data, response_var, predictor_vars, n_bootstrap = 100, robust = FALSE) {
  set.seed(123)
  
  bootstrap_results <- replicate(n_bootstrap, {
    bootstrap_sample <- train_data[sample(1:nrow(train_data), replace = TRUE), ]
    
    null_model <- lm(as.formula(paste(response_var, "~ 1")), data = bootstrap_sample)
    full_model <- lm(as.formula(paste(response_var, "~", paste(predictor_vars, collapse = " + "))), data = bootstrap_sample)
    
    forward_model <- stepAIC(
      null_model,
      scope = list(lower = null_model, upper = full_model),
      direction = "forward",
      trace = FALSE
    )
    
    if (robust) {
      selected_variables <- names(coef(forward_model))[-1]
      if (length(selected_variables) == 0) {
        return(NULL)
      }
      robust_formula <- as.formula(paste(response_var, "~", paste(selected_variables, collapse = " + ")))
      
      out <- tryCatch({
        lts_model <- ltsReg(robust_formula, data = bootstrap_sample)
        names(coef(lts_model))[-1]
      }, error = function(e) {
        NULL
      })
      return(out)
    } else {
      selected_variables <- names(coef(forward_model))[-1]
      return(selected_variables)
    }
  }, simplify = FALSE)
  
  if (robust) {
    selected_variables_count <- table(unlist(bootstrap_results))
  } else {
    selected_variables_count <- table(unlist(bootstrap_results))
  }
  return(selected_variables_count)
}

### 2.7 Build Final Model
build_final_model <- function(train_data, response_var, selected_predictors, robust = FALSE) {
  if (length(selected_predictors) == 0) {
    final_model <- lm(as.formula(paste(response_var, "~ 1")), data = train_data)
    return(final_model)
  }
  
  if (robust) {
    robust_formula <- as.formula(paste(response_var, "~", paste(selected_predictors, collapse = " + ")))
    final_model <- ltsReg(robust_formula, data = train_data)
  } else {
    final_model <- lm(as.formula(paste(response_var, "~", paste(selected_predictors, collapse = " + "))), data = train_data)
  }
  
  return(final_model)
}

### 2.8 Final Train Model Evaluation 
final_train_model_evaluation <- function(final_model, train_data, response_var, robust = FALSE) {
  if (robust) {
    predictors <- names(coef(final_model))[-1]
    if (length(predictors) > 0) {
      pred_matrix <- as.matrix(train_data[, predictors, drop = FALSE])
    } else {
      pred_matrix <- matrix(0, nrow = nrow(train_data), ncol = 0)
    }
    intercept <- coef(final_model)[1]
    coefficients <- coef(final_model)[-1]
    if (length(coefficients) == 0) { 
      train_predictions <- rep(intercept, nrow(train_data))
    } else {
      train_predictions <- intercept + pred_matrix %*% coefficients
    }
  } else {
    train_predictions <- predict(final_model, newdata = train_data)
  }
  
  actuals <- train_data[[response_var]]
  
  mse_train <- mean((actuals - train_predictions)^2)
  
  rmse_train <- sqrt(mse_train)
  
  summary_final_model <- summary(final_model)
  r_squared_adj_train <- summary_final_model$adj.r.squared
  
  return(list(R_squared_adj = r_squared_adj_train, RMSE = rmse_train))
}


### 2.9 Test Set Evaluation
test_set_evaluation <- function(final_model, test_data, response_var, robust = FALSE) {
  if (robust) {
    predictors <- names(coef(final_model))[-1]
    intercept <- coef(final_model)[1]
    coefficients <- coef(final_model)[-1]
    if (length(predictors) > 0) {
      pred_matrix <- as.matrix(test_data[, predictors, drop = FALSE])
      test_predictions <- intercept + pred_matrix %*% coefficients
    } else {
      test_predictions <- rep(intercept, nrow(test_data))
    }
  } else {
    test_predictions <- predict(final_model, newdata = test_data)
  }
  
  actuals <- test_data[[response_var]]
  
  mse_val <- mean((actuals - test_predictions)^2)
  
  rmse_val <- sqrt(mse_val)
  
  r_squared_test <- 1 - sum((actuals - test_predictions)^2) / sum((actuals - mean(actuals))^2)
  
  n <- length(actuals)
  p <- length(coef(final_model)) - 1  
  
  r_squared_adj_test <- 1 - (1 - r_squared_test) * (n - 1) / (n - p - 1)
  
  return(list(R_squared_adj = r_squared_adj_test, RMSE = rmse_val))
}



##### 3. Modeling Process for All Dependent Variables - NOT ROBUST VERSION #####

final_models <- list()
final_results_list <- list()

for (response_var in dep_vars) {
  cat("--------------------------------------------------\n")
  cat("Processing dependent variable:", response_var, "\n")
  
  # Bootstrap validation for variable selection
  bootstrap_res <- bootstrap_validation(
    train_data = train_sub, 
    response_var = response_var, 
    predictor_vars = ind_vars, 
    n_bootstrap = 100, 
    robust = FALSE
  )
  
  selected_vars <- names(bootstrap_res)[bootstrap_res > 50]
  cat("Variables selected (>50 times):", 
      if(length(selected_vars) > 0) paste(selected_vars, collapse = ", ") else "None", "\n")
  
  # Build final model
  final_model <- build_final_model(
    train_data = train_sub, 
    response_var = response_var, 
    selected_predictors = selected_vars, 
    robust = FALSE
  )
  
  final_models[[response_var]] <- final_model
  
  cat("Summary of final model for", response_var, ":\n")
  print(summary(final_model))
  
  # Train evaluation
  train_eval <- final_train_model_evaluation(
    final_model, 
    train_data = train_sub, 
    response_var = response_var, 
    robust = FALSE
  )
  
  # Test evaluation
  test_eval <- test_set_evaluation(
    final_model, 
    test_data = test_sub, 
    response_var = response_var, 
    robust = FALSE
  )
  
  result <- tibble(
    Response = response_var,
    Selected_Vars = if(length(selected_vars) > 0) paste(selected_vars, collapse = ", ") else "None",
    Train_R_squared_adj = round(train_eval$R_squared_adj, 4),
    Train_RMSE = round(train_eval$RMSE, 4),
    Test_R_squared_adj = round(test_eval$R_squared_adj, 4),
    Test_RMSE = round(test_eval$RMSE, 4)
  )
  
  final_results_list[[response_var]] <- result
}

final_results <- bind_rows(final_results_list)

### 3.1 Display the Final Results Table

final_results_table <- final_results %>%
  dplyr::select(Response, Selected_Vars, 
         Train_R_squared_adj, Train_RMSE,
         Test_R_squared_adj, Test_RMSE) 

print(
  final_results_table %>%
    kable(format = "markdown", caption = "Model Evaluation Results for Dependent Variables")
)

### 3.2 Create a Table of Coefficients with Significance Codes for Each Final Model

coeff_results_list <- lapply(names(final_models), function(response_var) {
  model_obj <- final_models[[response_var]]
  summ <- summary(model_obj)
  coefs <- summ$coefficients
  p_values <- coefs[,4]
  signif_codes <- sapply(p_values, significance_code)
  
  tibble(
    Response = response_var,
    Term = rownames(coefs),
    Estimate = round(coefs[,1], 4),
    Std_Error = round(coefs[,2], 4),
    t_value = round(coefs[,3], 4),
    p_value = round(coefs[,4], 4),
    Signif = signif_codes
  )
})

coeff_results <- bind_rows(coeff_results_list)

### 3.3 Pivot the Coefficients Table to Have Independent Variables as Columns


coeff_pivot <- coeff_results %>%
  mutate(Term = ifelse(Term == "(Intercept)", "Intercept", Term)) %>%
  mutate(Estimate_Signif = ifelse(Signif != "", 
                                  paste0(Estimate, Signif), 
                                  as.character(Estimate))) %>%
  dplyr::select(Response, Term, Estimate_Signif) %>%
  pivot_wider(names_from = Term, values_from = Estimate_Signif)

required_columns <- c("Intercept", ind_vars)
missing_columns <- setdiff(required_columns, colnames(coeff_pivot))
if(length(missing_columns) > 0){
  coeff_pivot[missing_columns] <- NA
}

coeff_pivot <- coeff_pivot %>%
  dplyr::select(Response, Intercept, all_of(ind_vars))


### 3.4 Display the Coefficients Table

print(
  coeff_pivot %>%
    kable(format = "markdown", caption = "Coefficients and Significance Codes for Each Final Model")
)


####### 4. Modeling Process for All Dependent Variables -  ROBUST VERSION ######

final_models <- list()  
final_results_list <- list()  

for (response_var in dep_vars) {
  cat("--------------------------------------------------\n")
  cat("Processing dependent variable:", response_var, "\n")
  
  # Bootstrap validation for variable selection
  bootstrap_res <- bootstrap_validation(
    train_data = train_sub, 
    response_var = response_var, 
    predictor_vars = ind_vars, 
    n_bootstrap = 100, 
    robust = TRUE
  )
  
  selected_vars <- names(bootstrap_res)[bootstrap_res > 50]
  cat("Variables selected (>50 times):", 
      if(length(selected_vars) > 0) paste(selected_vars, collapse = ", ") else "None", "\n")
  
  # Build final model
  final_model <- build_final_model(
    train_data = train_sub, 
    response_var = response_var, 
    selected_predictors = selected_vars, 
    robust = TRUE
  )
  
  final_models[[response_var]] <- final_model
  
  cat("Summary of final model for", response_var, ":\n")
  print(summary(final_model))
  
  # Train evaluation
  train_eval <- final_train_model_evaluation(
    final_model, 
    train_data = train_sub, 
    response_var = response_var, 
    robust = TRUE
  )
  
  # Test evaluation
  test_eval <- test_set_evaluation(
    final_model, 
    test_data = test_sub, 
    response_var = response_var, 
    robust = TRUE
  )
  
  result <- tibble(
    Response = response_var,
    Selected_Vars = if(length(selected_vars) > 0) paste(selected_vars, collapse = ", ") else "None",
    Train_R_squared_adj = round(train_eval$R_squared_adj, 4),
    Train_RMSE = round(train_eval$RMSE, 4),
    Test_R_squared_adj = round(test_eval$R_squared_adj, 4),
    Test_RMSE = round(test_eval$RMSE, 4)
  )
  
  final_results_list[[response_var]] <- result
}

final_results <- bind_rows(final_results_list)

### 4.1 Display the Final Results Table

final_results_table <- final_results %>%
  dplyr::select(Response, Selected_Vars, 
         Train_R_squared_adj, Train_RMSE,
         Test_R_squared_adj, Test_RMSE) 

# Print the final results table
print(
  final_results_table %>%
    kable(format = "markdown", caption = "Model Evaluation Results for Dependent Variables")
)

### 4.2 Create a Table of Coefficients with Significance Codes for Each Final Model

coeff_results_list <- lapply(names(final_models), function(response_var) {
  model_obj <- final_models[[response_var]]
  summ <- summary(model_obj)
  coefs <- summ$coefficients
  p_values <- coefs[,4]
  signif_codes <- sapply(p_values, significance_code)
  
  tibble(
    Response = response_var,
    Term = rownames(coefs),
    Estimate = round(coefs[,1], 4),
    Std_Error = round(coefs[,2], 4),
    t_value = round(coefs[,3], 4),
    p_value = round(coefs[,4], 4),
    Signif = signif_codes
  )
})

coeff_results <- bind_rows(coeff_results_list)

### 4.3 Pivot the Coefficients Table to Have Independent Variables as Columns

coeff_pivot <- coeff_results %>%
  mutate(Term = ifelse(Term == "(Intercept)", "Intercept", Term)) %>%
  mutate(Estimate_Signif = ifelse(Signif != "", 
                                  paste0(Estimate, Signif), 
                                  as.character(Estimate))) %>%
  dplyr::select(Response, Term, Estimate_Signif) %>%
  pivot_wider(names_from = Term, values_from = Estimate_Signif)

required_columns <- c("Intercept", ind_vars)
missing_columns <- setdiff(required_columns, colnames(coeff_pivot))
if(length(missing_columns) > 0){
  coeff_pivot[missing_columns] <- NA
}

coeff_pivot <- coeff_pivot %>%
  dplyr::select(Response, Intercept, all_of(ind_vars))

### 4.4 Display the Coefficients Table

print(
  coeff_pivot %>%
    kable(format = "markdown", caption = "Coefficients and Significance Codes for Each Final Model")
)

############# 5. Outlier Detection and Summary for All Variables ###############

variables <- c("Courage", "Creativity", "Curiosity", "Gratitude", "Honesty", "Hope", "Humor", 
               "Judgment", "Kindness", "Leadership", "Love_of_Learning", "Love", "Humility", 
               "Perseverance", "Perspective", "Prudence", "Self_Regulation", "Social_Intelligence", 
               "Spirituality", "Teamwork", "Zest")

outlier_summary <- data.frame(
  Variable = character(),
  RobustDistance = numeric(),
  Residuals = numeric(),
  Outlier = factor(levels = c("Leverage & Residual", "Leverage", "Residual", "Regular"))
)

for (var in variables) {
  
  model <- final_models[[var]]
  
  if (is.null(model) | !inherits(model, "ltsReg")) {
    selected_vars <- names(coef(model))[-1]  
    formula <- as.formula(paste(var, "~", paste(selected_vars, collapse = " + ")))
    
    robust_model <- ltsReg(formula, data = train_sub)
    
    final_models[[paste0(var, "_Robust")]] <- robust_model
    cat(paste("Robust LTS regression model for", var, "has been created.\n"))
  } else {
    robust_model <- model
  }
  
  selected_vars <- names(coef(robust_model))[-1]  
  
  missing_vars <- setdiff(selected_vars, colnames(train_sub))
  if (length(missing_vars) > 0) {
    stop(paste("The following predictors for", var, "are missing in 'train_sub':", 
               paste(missing_vars, collapse = ", ")))
  }
  
  train_data <- train_sub[, c(var, selected_vars)]
  
  mcd_result <- covMcd(train_data[, selected_vars])  
  robust_center <- mcd_result$center
  robust_cov <- mcd_result$cov
  
  robust_distances <- sqrt(mahalanobis(
    x = train_data[, selected_vars],
    center = robust_center,
    cov = robust_cov
  ))
  
  residuals <- residuals(robust_model)
  std_residuals <- residuals / sd(residuals)
  
  df <- length(selected_vars)  # Gradi di libertà
  threshold_robust_distance <- sqrt(qchisq(0.975, df = df))  
  threshold_residuals <- 2.5  
  
  outlier_distance <- robust_distances > threshold_robust_distance
  outlier_residuals <- abs(std_residuals) > threshold_residuals
  
  outlier_data <- data.frame(
    Variable = var,
    RobustDistance = robust_distances,
    Residuals = std_residuals,
    Outlier = factor(
      ifelse(outlier_distance & outlier_residuals, "Leverage & Residual",
             ifelse(outlier_distance, "Leverage",
                    ifelse(outlier_residuals, "Residual", "Regular"))),
      levels = c("Leverage & Residual", "Leverage", "Residual", "Regular")
    )
  )
  
  outlier_summary <- rbind(outlier_summary, outlier_data)
}


cat("Outlier Summary for All Variables:\n")
table(outlier_summary$Variable, outlier_summary$Outlier)

head(outlier_summary)


######## 5. Outlier Detection and Outlier Map for 'Courage' Regression #########

### 5.1 Prepare the Robust Regression Model for 'Courage'

courage_model <- final_models[["Courage"]]

if (is.null(courage_model) | !inherits(courage_model, "ltsReg")) {
  selected_vars_courage <- names(coef(courage_model))[-1]  
  formula_courage <- as.formula(paste("Courage ~", paste(selected_vars_courage, collapse = " + ")))
  
  final_model_r_courage <- ltsReg(formula_courage, data = train_sub)
  
  final_models[["Courage_Robust"]] <- final_model_r_courage
  
  cat("Robust LTS regression model for 'Courage' has been created and stored as 'Courage_Robust' in 'final_models'.\n")
} else {
  final_model_r_courage <- courage_model
  cat("Robust LTS regression model for 'Courage' is already available in 'final_models'.\n")
}


### 5.2 Define Selected Variables for 'Courage'

selected_vars_courage <- names(coef(final_model_r_courage))[-1]

missing_vars <- setdiff(selected_vars_courage, colnames(train_sub))
if (length(missing_vars) > 0) {
  stop(paste("The following selected predictors for 'Courage' are missing in 'train_sub':",
             paste(missing_vars, collapse = ", ")))
}


### 5.3 Outlier Detection

train_courage <- train_sub[, c("Courage", selected_vars_courage)]

mcd_result_courage <- covMcd(train_courage[, selected_vars_courage])  # Robust MCD estimation
robust_center_courage <- mcd_result_courage$center
robust_cov_courage <- mcd_result_courage$cov

robust_distances_courage <- sqrt(mahalanobis(
  x = train_courage[, selected_vars_courage],
  center = robust_center_courage,
  cov = robust_cov_courage
))

residuals_courage <- residuals(final_model_r_courage)
std_residuals_courage <- residuals_courage / sd(residuals_courage)

df_courage <- length(selected_vars_courage)  
threshold_robust_distance_courage <- sqrt(qchisq(0.975, df = df_courage))  
threshold_residuals_courage <- 2.5  

outlier_distance_courage <- robust_distances_courage > threshold_robust_distance_courage
outlier_residuals_courage <- abs(std_residuals_courage) > threshold_residuals_courage

outlier_data_courage <- data.frame(
  RobustDistance = robust_distances_courage,
  Residuals = std_residuals_courage,
  Outlier = factor(
    ifelse(outlier_distance_courage & outlier_residuals_courage, "Leverage & Residual",
           ifelse(outlier_distance_courage, "Leverage",
                  ifelse(outlier_residuals_courage, "Residual", "Regular")))
  )
)


### 5.4 Outlier Summary

cat("Outlier Summary for 'Courage' Regression Model:\n")
print(table(outlier_data_courage$Outlier))


### 5.5 Outlier Map Visualization

train_courage$ID <- rownames(train_courage)
outlier_data_courage <- outlier_data_courage %>%
  dplyr::mutate(ID = train_courage$ID)

outlier_colors <- c(
  "Regular" = "deeppink",
  "Leverage" = "orange",
  "Residual" = "green",
  "Leverage & Residual" = "purple"
)

outlier_map_courage <- ggplot(outlier_data_courage, aes(x = RobustDistance, y = Residuals, color = Outlier)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_hline(yintercept = c(-threshold_residuals_courage, threshold_residuals_courage), 
             color = "red", linetype = "dashed") +
  geom_vline(xintercept = threshold_robust_distance_courage, 
             color = "blue", linetype = "dashed") +
  labs(
    title = "Outlier Map (MCD-Based) for 'Courage' Regression Model",
    x = "Robust Distance (MCD)",
    y = "Standardized Residuals",
    color = "Outlier Type"
  ) +
  scale_color_manual(values = outlier_colors) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  ) +
  geom_text(
    data = subset(outlier_data_courage, Outlier != "Regular"),
    aes(label = ID),
    vjust = -1,
    size = 3,
    color = "deeppink"
  )

outlier_map_courage


### 5.6 Display a Table of Outliers

if(any(outlier_data_courage$Outlier != "Regular")){
  outliers_courage <- outlier_data_courage %>%
    dplyr::filter(Outlier != "Regular") %>%
    dplyr::select(RobustDistance, Residuals, Outlier, ID)
  
  outliers_courage <- train_courage %>%
    dplyr::filter(ID %in% outliers_courage$ID) %>%
    dplyr::left_join(outliers_courage, by = "ID") %>%
    dplyr::select(ID, Courage, all_of(selected_vars_courage), RobustDistance, Residuals, Outlier)
  
  print(
    outliers_courage %>%
      knitr::kable(format = "markdown", caption = "Outliers in 'Courage' Regression Model")
  )
} else {
  cat("No outliers detected for 'Courage' based on the defined thresholds.\n")
}

########## -------------------------------------------------------- ############

############################ MODEL BASED CLUSTERING ############################

# 1) PCA
# 2) Mclust and Tclust
# 3) UMAP with K-means Clustering and Robust with Trimmed K-Means


data_c <- data

### Dropping unnecessary columns (1,2 and from 7 to 182)
data_c <- data_c[, -c(1, 2, 3, 4)]

data_normalized_robust_c <- data_c %>%
  mutate(across(where(is.numeric), ~ (.-min(., na.rm = TRUE)) / 
                  (max(., na.rm = TRUE) - min(., na.rm = TRUE))))

set.seed(122)
########## -------------------------------------------------------- ############
#data_c
###################### 1) PCA: Reducing Dimensionality #########################

# Assuming 'data_normalized_robust' is your preprocessed data frame
pca_result <- prcomp(data_normalized_robust_c, center = TRUE, scale. = TRUE)

# Summary of PCA to check explained variance
summary(pca_result)

# Scree plot to visualize explained variance
plot(pca_result, type = "l", main = "Scree Plot")

# Retain the first 7 principal components (as per your initial analysis)
pca_data <- as.data.frame(pca_result$x[, 1:7])  

################################## 2) Clustering ###############################

################# a) Non-Robust Gaussian Mixture Model (Mclust) ################

# Fit Mclust model
fit_mclust <- Mclust(data = pca_data)
summary(fit_mclust)

# Extract cluster assignments
mclust_clusters <- fit_mclust$classification

################## b) Robust Gaussian Mixture Model (tclust) ###################

# Determine number of clusters from Mclust for consistency
num_clusters <- fit_mclust$G

# Fit tclust model --> alpha=0.05 trims 5% of observations as potential outliers
robust_fit <- tclust(pca_data, k = num_clusters, alpha = 0.05, restr.fact = 1, equal.weights = FALSE)

# Extract cluster assignments from tclust
robust_clusters <- robust_fit$cluster

# Assuming 'robust_clusters' contains the cluster assignments from the tclust model
cluster_counts <- table(robust_clusters)

# Display the number of observations per cluster
print(cluster_counts)


######################## 3) Compute Silhouette Scores ##########################

# Compute distance matrix using Euclidean distance
dist_matrix <- dist(pca_data)

# Silhouette for Mclust
sil_mclust <- silhouette(mclust_clusters, dist_matrix)
avg_sil_mclust <- mean(sil_mclust[, 3])

# Silhouette for tclust
sil_robust <- silhouette(robust_clusters, dist_matrix)
avg_sil_robust <- mean(sil_robust[, 3])

# Print average silhouette scores
cat("Average Silhouette Score (Mclust): ", round(avg_sil_mclust, 3), "\n")
cat("Average Silhouette Score (tclust):  ", round(avg_sil_robust, 3), "\n\n")

##################### 4) Plotting Clustering Results ############################

# Prepare data for plotting
plot_data <- pca_data
plot_data$Mclust_Cluster <- as.factor(mclust_clusters)

# Handle robust clusters: 0 indicates outliers
plot_data$Robust_Cluster <- ifelse(robust_clusters == 0, "Outliers", as.character(robust_clusters))
plot_data$Robust_Cluster <- factor(plot_data$Robust_Cluster, 
                                   levels = c("Outliers", as.character(1:num_clusters)))

# Select first two principal components for 2D plotting
plot_data$PC1 <- pca_result$x[, 1]
plot_data$PC2 <- pca_result$x[, 2]

# Plot for Mclust
p1 <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Mclust_Cluster)) + 
  geom_point(alpha = 0.8, size = 3.5, shape = 16) + 
  stat_ellipse(aes(group = Mclust_Cluster), level = 0.95, linetype = "solid") +  # Add ellipses
  labs(title = paste(toupper("Mclust Clustering"), "\nSilhouette Score:", round(avg_sil_mclust, 3)),
       color = "Cluster") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


colors_vector <- c("Outliers" = "#fede00", "1" = "#e57f84", "2" = "#4297a0")

# Plot for tclust
p2 <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Robust_Cluster)) + 
  geom_point(alpha = 0.8, size = 3.5, shape = 16) + 
  scale_color_manual(values = colors_vector) + 
  stat_ellipse(aes(group = Robust_Cluster), level = 0.95, linetype = "solid") +  # Add ellipses
  ggtitle(paste("TCLUST CLUSTERING\nSilhouette Score:", round(avg_sil_robust, 3))) + 
  theme_minimal() + 
  xlab("Component 1") + 
  ylab("Component 2") + 
  labs(color = "Cluster") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Display plots
print(p1)
print(p2)



########################### Full Code from UMAP Clustering to Radar Charts ###########################

data_normalized_robust_c <- na.omit(data_normalized_robust_c)

######## 4) Dimensionality Reduction with UMAP and K-means Clustering ##########

### Direct Application of UMAP to normalized data
umap_result <- umap(data_normalized_robust_c)
umap_data <- data.frame(umap_result$layout)
colnames(umap_data) <- c("X1", "X2")

# Clustering with K-means --> 3 clusters
kmeans_result <- kmeans(umap_data, centers = 3)
umap_data$cluster <- as.factor(kmeans_result$cluster)

### a) Compute average values for each variable in each k-means cluster
data_with_kmeans_clusters <- data_normalized_robust_c %>% mutate(cluster = umap_data$cluster)

# Summarize numeric columns by cluster
kmeans_cluster_means <- data_with_kmeans_clusters %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

cat("K-means Cluster-wise Averages:\n")
print(kmeans_cluster_means, n = Inf, width = Inf)

### b) Cluster dimension for K-means
kmeans_cluster_sizes <- umap_data %>%
  group_by(cluster) %>%
  summarise(size = n())

cat("\nK-means Cluster Sizes:\n")
print(kmeans_cluster_sizes)

### c) Visualization of UMAP with different colors for each cluster
ggplot(umap_data, aes(x = X1, y = X2, color = cluster)) + 
  geom_point(alpha = 0.8, size = 3.5, shape = 16) + 
  stat_ellipse(aes(color = cluster), level = 0.95) +  # Aggiungi ellissi di confidenza
  ggtitle("UMAP with Clustering K-means") + 
  theme_minimal() + 
  xlab("Component 1") + 
  ylab("Component 2") + 
  labs(color = "Cluster") + 
  scale_color_manual(values = c("1" = "#e57f84", "2" = "#4297a0", "3" = "#b99095")) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

### d) Silhouette for K-means
silhouette_kmeans <- silhouette(as.integer(umap_data$cluster), dist(umap_data[, c("X1", "X2")]))
cat("K-means Silhouette Average Width:", mean(silhouette_kmeans[, 3]), "\n")

### e) Radar Chart for K-means Clusters
original_data <- data_normalized_robust_c

# Identify variables: exclude any clustering columns
variables <- setdiff(names(original_data), c("cluster", "trimmed_cluster"))

# Prepare data for radar chart (K-means)
radar_data_kmeans <- kmeans_cluster_means %>%
  column_to_rownames("cluster") %>%
  as.data.frame()

# Reorder columns to match 'variables' exactly
radar_data_kmeans <- radar_data_kmeans[, variables, drop = FALSE]

max_vals_kmeans <- apply(original_data[, variables], 2, max, na.rm = TRUE)
min_vals_kmeans <- apply(original_data[, variables], 2, min, na.rm = TRUE)

radar_limits_kmeans <- rbind(max_vals_kmeans, min_vals_kmeans)
radar_data_combined_kmeans <- rbind(radar_limits_kmeans, radar_data_kmeans)

# Ensure numeric
radar_data_combined_kmeans <- radar_data_combined_kmeans %>%
  mutate(across(everything(), as.numeric))

# Define colors for K-means clusters
cluster_colors_kmeans <- c("1" = "#e57f84", "2" = "#4297a0", "3" = "#b99095")

# Plot the radar chart for K-means
radarchart(radar_data_combined_kmeans, 
           axistype = 1, 
           pcol = cluster_colors_kmeans,
           plwd = 5, 
           plty = 1,
           title = "K-means Cluster-wise Variable Profiles",
           cglcol = "grey", 
           cglty = 1, 
           axislabcol = "grey", 
           caxislabels = rep("", ncol(radar_data_combined_kmeans)), 
           cglwd = 0.8,
           vlcex = 0.8)

legend("right", 
       legend = paste("Cluster", rownames(radar_data_kmeans)), 
       col = cluster_colors_kmeans, 
       lwd = 4, 
       cex = 0.8, 
       bty = "n",
       title = "Legend",
       title.cex = 0.9)


################# ------- Robust with Trimmed K-Means ------- ##################

# Trimming Parameters --> Trimming Percentage of 5%
alpha <- 0.05 

# Application of Trimmed K-means to UMAP results
fit_trimmed_kmeans <- tclust(
  x = umap_data[, c("X1", "X2")], 
  k = 3,                          
  alpha = alpha,                 
  restr.fact = 1
)

# Assign trimmed k-means cluster labels to the UMAP data
umap_data$trimmed_cluster <- as.factor(fit_trimmed_kmeans$cluster)

### Rename Cluster 0 to "outliers"
umap_data <- umap_data %>%
  mutate(trimmed_cluster = dplyr::recode(trimmed_cluster, `0` = "outliers"))

### a) Compute Average Values for Each Variable in Each Trimmed Cluster

# Append trimmed_cluster to the original data with k-means clusters
data_with_trimmed_clusters <- data_with_kmeans_clusters %>%
  mutate(trimmed_cluster = umap_data$trimmed_cluster)

# Summarize numeric columns by cluster
trimmed_cluster_means <- data_with_trimmed_clusters %>%
  group_by(trimmed_cluster) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

cat("Trimmed K-means Cluster-wise Averages:\n")
print(trimmed_cluster_means, n = Inf, width = Inf)

### b) Cluster dimension for trimmed clusters
trimmed_cluster_sizes <- umap_data %>%
  group_by(trimmed_cluster) %>%
  summarise(size = n())

print("Trimmed K-means Cluster Sizes:")
print(trimmed_cluster_sizes)

### c) Visualization of UMAP with different colors for each Trimmed Cluster
ggplot(umap_data, aes(x = X1, y = X2, color = trimmed_cluster)) + 
  geom_point(alpha = 0.8, size = 3.5, shape = 16) + 
  ggtitle("UMAP with Trimmed K-means Clustering") + 
  theme_bw() +  # Usa uno sfondo bianco
  xlab("Component 1") + 
  ylab("Component 2") + 
  labs(color = "Trimmed Cluster") + 
  scale_color_manual(values = c("1" = "#e57f84", "2" = "#4297a0", "3" = "#b99095", "outliers" = "#fede00")) + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_blank()
  ) +
  stat_ellipse(aes(group = trimmed_cluster), level = 0.95, geom = "polygon", alpha = 0.2, linetype = "solid", fill = NA)  # Rimuove lo sfondo grigio dell'ellisse

### d) Silhouette for trimmed clusters
silhouette_trimmed <- silhouette(as.integer(as.factor(umap_data$trimmed_cluster)), dist(umap_data[, c("X1", "X2")]))
cat("Trimmed K-means Silhouette Average Width:", mean(silhouette_trimmed[, 3]), "\n")

### e) Radar Chart for Trimmed K-means Clusters
trimmed_radar_data <- trimmed_cluster_means %>%
  column_to_rownames("trimmed_cluster") %>%
  as.data.frame()

# Ensure that trimmed_radar_data has the same variables in the same order
trimmed_radar_data <- trimmed_radar_data[, variables, drop = FALSE]

max_vals_trimmed <- apply(original_data[, variables], 2, max, na.rm = TRUE)
min_vals_trimmed <- apply(original_data[, variables], 2, min, na.rm = TRUE)

radar_limits_trimmed <- rbind(max_vals_trimmed, min_vals_trimmed)
radar_data_combined_trimmed <- rbind(radar_limits_trimmed, trimmed_radar_data)

# Ensure numeric
radar_data_combined_trimmed <- radar_data_combined_trimmed %>%
  mutate(across(everything(), as.numeric))

# Identify trimmed clusters
trimmed_cluster_levels <- rownames(trimmed_radar_data)

# Define colors for Trimmed K-means clusters
trimmed_colors <- c("1" = "#e57f84", "2" = "#4297a0", "3" = "#b99095", "outliers" = "#fede00")

# Plot the radar chart for Trimmed K-means
radarchart(radar_data_combined_trimmed, 
           axistype = 1, 
           pcol = trimmed_colors[trimmed_cluster_levels],
           plwd = 5, 
           plty = 1,
           title = "Trimmed K-means Cluster-wise Variable Profiles",
           cglcol = "grey", 
           cglty = 1, 
           axislabcol = "grey", 
           caxislabels = rep("", ncol(radar_data_combined_trimmed)), 
           cglwd = 0.8,
           vlcex = 0.8)

legend("right", 
       legend = paste("Cluster", trimmed_cluster_levels), 
       col = trimmed_colors[trimmed_cluster_levels], 
       lwd = 4, 
       cex = 0.8, 
       bty = "n",
       title = "Legend",
       title.cex = 0.9)
