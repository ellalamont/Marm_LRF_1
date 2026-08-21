# Validate function similar to rms::validate (https://rdrr.io/cran/rms/src/R/validate.lrm.s) but that works with firth logistic regression
# Also Calibrate functions
# E. Lamont
# 8/4/26

library(logistf)
library(pROC)

validate_firth <- function(formula, data, outcome, B = 500, seed = 42){
  
  ## Bootstrapping optimism correction for a firth logistic regression model
  ## Uses stratified bootstrapping
  ## Based on the rms::validate function 
  
  set.seed(seed)
  y <- data[[outcome]]
  
  #-----------------------------
  # Helper function to calculate metrics
  ## Takes a fitted logistf model (fit) and new data to predict on
  #-----------------------------
  get_metrics <- function(fit, newdata){
    p <- predict(fit, newdata = newdata, type = "response") # predicted probabilities
    lp <- predict(fit, newdata = newdata, type = "link") # linear predictor for calibration slope. Log-odds prediction before converting it into a probability
    y <- newdata[[outcome]] # true outcomes
    auc <- as.numeric(pROC::auc(y, p))
    dxy <- 2*auc - 1
    brier <- mean((y - p)^2)
    # Calibration checks whether observed outcome is similar to predicted probability
    cal <- try( # Calibration intercept and slope. This is weak calibration.
      logistf(y ~ lp, pl = FALSE, # Fits a new firth logistic regression to check if model agrees with reality
              control = logistf.control(maxit = 1000, maxstep = 0.5)), 
      silent = TRUE)
    if(inherits(cal,"try-error") || # || is logical OR of length 1
      is.na(coef(cal)[2]) ||
      abs(coef(cal)[2]) > 10){
      return(c(AUC = NA, Dxy = NA, Brier = NA, Intercept = NA, Slope = NA))
    }
    intercept <- coef(cal)[1] # Perfect would be 0. -1 means predictions too high
    slope <- coef(cal)[2] # Perfect would be 1. <1 is overfit
    return(c(AUC = auc, 
             Dxy = dxy,
             Brier = brier,
             Intercept = intercept,
             Slope = slope))
  }
  
  #-----------------------------
  # Apparent performance
  #-----------------------------
  fit <- logistf(formula, data = data) # Fits model using the full dataset
  apparent <- get_metrics(fit, data) # Evaluates model of full dataset. This is the optimistic performance
  
  #-----------------------------
  # Bootstrap
  #-----------------------------
  optimism <- matrix(NA, nrow = B, ncol = length(apparent)) # Creats an empty matrix
  colnames(optimism) <- names(apparent)
  relapse <- which(y == 1)
  cure <- which(y == 0)
  successful <- 0 # Counts successful bootstrap models

  for(i in 1:B){ # stratified bootstrap loop
    boot_ids <- c( # Resamples relapse and cure patients separately to maintain a class balance
      sample(relapse,
             length(relapse), # Resamples with replacement the same number of relapses as original dataset
             replace = TRUE),
      sample(cure,
             length(cure), # Resamples with replacement the same number of cures as original datset
             replace = TRUE))
    boot_data <- data[boot_ids, ] # Create the specific bootstrap dataset for this loop
    fit_boot <- try( # New firth model with the bootstrap dataset
      logistf(formula, data = boot_data, pl = FALSE, control = logistf.control(maxit = 1000, maxstep = 0.5)), silent = TRUE)
    if(inherits(fit_boot,"try-error"))
      next
    successful <- successful + 1
    train_perf <- get_metrics(fit_boot, boot_data) # Evaluates bootstrap model on boostrap dataset
    test_perf <- get_metrics(fit_boot, data) # Evaluates bootstrap model on full dataset. This estimates generalization
    optimism[successful,] <- train_perf - test_perf # Measures the amount of optimism for this bootstrap model
  } # End of loop
  
  optimism <- optimism[1:successful,,drop=FALSE]
  mean_optimism <- colMeans(optimism, na.rm = TRUE) # average the optimism across all bootstrap models.
  corrected <- apparent - mean_optimism # Gets the corrected optimism
  
  results <- data.frame(
    Metric = names(apparent),
    Apparent = as.numeric(apparent),
    Optimism = as.numeric(mean_optimism),
    Corrected = as.numeric(corrected))
  
  # Store additional information
  attr(results,"successful_bootstraps") <- successful
  attr(results,"requested_bootstraps") <- B

  return(results)
}

# ROC is discrimination - Can the model rank patients correctly?
# Calibration is: are the predicted probabilites numerically correct?
# Ex. model predicts 100 patients have a 70% chance of relapse. Perfect calibration would be when 70 patients actually relapse

calibrate_firth <- function(formula, data, B = 200,
                            grid = seq(0.01,0.99,length=200),
                            span = 0.75, seed = 42){
  
  ## Using LOESS calibration curve. This is weak calibration
  
  
  if(!is.null(seed))
    set.seed(seed)
  
  ## Fit original model
  fit <- logistf(formula, data=data, pl=FALSE) # Fit the original model on all the data
  p_orig <- predict(fit, type="response") # Get the predicted probabilities
  yname <- all.vars(formula)[1] # Dependent variable
  y <- data[[yname]]
  
  #-----------------------------
  # Apparent calibration
  ## Calibration on same data used to build the model
  #-----------------------------
  app_fit <- loess(y ~ p_orig, span=span, degree=1) # Fit LOESS calibration curve
  apparent <- pmin(pmax(predict(app_fit, newdata=data.frame(p_orig=grid)), 0),1) # Determines the observed probability at each predicted risk. Keeps its all between 0 and 1.
  
  #-----------------------------
  # Bootstrap
  #-----------------------------
  optimism <- matrix(NA,B,length(grid)) # New matrix with each row being one bootstrap calibration correction
  successful <- 0 # Count successful models
  for(b in seq_len(B)){
    idx <- sample(seq_len(nrow(data)), replace=TRUE) # Sample patients with replacement (not stratified)
    boot_data <- data[idx,] # Create bootstrap dataset from sampled patients
    fit_boot <- try( # Fit the firth regression model using the bootstrap dataset
      logistf(formula, data=boot_data, pl=FALSE, control=logistf.control(maxit=1000)),
      silent=TRUE)
    if(inherits(fit_boot,"try-error"))
      next
    successful <- successful + 1

    p_boot <- predict(fit_boot, newdata = boot_data, type="response") # Predict the bootstrap model on the bootstrap dataset for this loop. This is the training performance
    p_test <- predict(fit_boot, newdata=data, type="response") # Predict the bootstrap model on the full dataset. This is the test performance
    
    # LOESS: Creates a smooth estimate of relationship between what model predicted and what actually happened. Gives an estimated calibration at every point along the range, not just individual data points.
    boot_curve <- try( # Fit LOESS calibration curve on the bootstrap data
      loess(boot_data[[yname]] ~ p_boot, span=span, degree=1), silent=TRUE)
    
    test_curve <- try( # Fit LOESS calibration curve on the full data
      loess(y ~ p_test, span=span, degree=1), silent=TRUE)
    
    if(inherits(boot_curve,"try-error") ||
       inherits(test_curve,"try-error"))
      next
    
    # At each predicted probability, determine what the LOESS curve estimates the observed probability is. 
    boot_cal <- pmin(pmax( # Predict on bootstrap data
      predict(boot_curve, newdata=data.frame(p_boot=grid)), 0),1)
    
    test_cal <- pmin(pmax( # Predict on full data
      predict(test_curve, newdata=data.frame(p_test=grid)), 0),1)
    
    optimism[successful,] <- boot_cal - test_cal # Calculate the optimsim
  } # End bootstrap loop
  
  optimism <- optimism[ complete.cases(optimism), , drop=FALSE]
  mean_optimism <- colMeans(optimism) # Average optimism for all bootstrap samples
  corrected <- pmin(pmax(apparent - mean_optimism, 0),1) # Get the corrected values
  
  results <- data.frame(Predicted = grid,
                    Apparent = apparent,
                    Corrected = corrected,
                    Ideal = grid)
  
  attr(results, "successful_bootstraps") <- nrow(optimism)
  class(results) <- c("calibrate_firth", class(results))
  results
}



plot.calibrate_firth <- function(x,..., main = "Calibration plot"){
  plot(x$Predicted,
       x$Corrected, # Estimate of how model would perform on new data
       type="l",
       lwd=3,
       col="blue",
       xlim=c(0,1),
       ylim=c(0,1),
       xlab="Predicted probability",
       ylab="Observed probability",
       main = main)
  lines(x$Predicted,
        x$Apparent,
        lty=2,
        lwd=2)
  lines(x$Predicted, x$Ideal, col="grey42", lty=1) # The ideal calibration line
  
  legend("topleft",
         c("Corrected",
           "Apparent",
           "Ideal"),
         lwd=c(3,2,1),
         lty=c(1,2,1),
         cex = 0.8,
         col=c("blue","black","grey42"),
         bty="n")
}


