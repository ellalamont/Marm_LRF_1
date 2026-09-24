# Generate Concensus Genes through 4 different methods
# E. Lamont
# 8/6/26

source("Import_data.R")
source("Function_Firth_Validate.R")
library(rmda)
library(caret)

###########################################################
##################### ORGANIZE DATA #######################

tmp_log2tpm <- GoodSamples60_log2tpmf %>% 
  dplyr::select(!contains("Rv")) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("SampleID2")

model_df <- GoodSamples60_pipeSummary %>% 
  mutate(Outcome = ifelse(Outcome2 == 0, "control", "case")) %>%
  dplyr::select(SampleID2, Outcome) %>%
  mutate(Outcome = factor(Outcome, levels = c("case", "control"))) %>%
  inner_join(tmp_log2tpm, by = "SampleID2") %>%
  column_to_rownames("SampleID2") %>%
  na.omit()

# 88 samples
# 4030 genes

model_X <- model_df %>% dplyr::select(-Outcome)
model_Y <- model_df$Outcome

###########################################################
############# LOGSTIC REGRESSION MODEL (CV) ###############

# Define training control
set.seed(23) 
LR_train.control <- trainControl(
  method = "repeatedcv", 
  classProbs = TRUE, 
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  number = 4,
  repeats = 30)
LR_lassoGrid <- expand.grid(
  alpha = 1,
  lambda = 10^seq(-3, 1, length = 10))

# Train the model
set.seed(23)
log2tpm_LR_model <- caret::train(
  x = model_X,
  y = model_Y,
  method = "glmnet",
  metric = "ROC",
  tuneGrid = LR_lassoGrid,
  preProcess = c("nzv", "center", "scale"),
  trControl = LR_train.control)
print(log2tpm_LR_model)

# Collect the genes used
log2tpm_LR_coef_df <- coef(log2tpm_LR_model$finalModel, s = log2tpm_LR_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
log2tpm_LR_coef_df <- log2tpm_LR_coef_df[log2tpm_LR_coef_df[, 1] != 0, , drop = FALSE]
colnames(log2tpm_LR_coef_df)[1] <- "V1"
log2tpm_LR_Genes <- rownames(log2tpm_LR_coef_df)
log2tpm_LR_Genes
# [1] "(Intercept)" "Rv0451c"     "Rv0459"      "Rv0793"      "Rv1062"      "Rv1440"     
# [7] "Rv3851"   

# Predict case probability using final fitted lasso model (on the same dataset)
log2tpm_LR_final_predictions <- predict(log2tpm_LR_model, newdata = model_X, type = "prob")

# Combine predictions with outcomes
log2tpm_LR_final_predictions_df <- data.frame(
  SampleID2 = rownames(model_df),
  Outcome = model_Y,
  case_probability = log2tpm_LR_final_predictions$case) %>%
  mutate(Predicted = ifelse(case_probability >= 0.5, "case", "control")) %>%
  mutate(Predicted = factor(Predicted, levels = c("case", "control"))) %>%
  mutate(Outcome2 = ifelse(Outcome == "case", 1, 0))

# Plot the model ROC
log2tpm_LR_ROC <- roc(response = log2tpm_LR_model$pred$obs,
              predictor = log2tpm_LR_model$pred$case,
              levels = c("control", "case"),  # (Control, case)
              direction = "<")
plot(log2tpm_LR_ROC, print.auc = TRUE)
auc(log2tpm_LR_ROC)
# Area under the curve: 0.4911

###########################################################
########## LOGSTIC REGRESSION MODEL (BOOTSTRAP) ###########

# Define training control
set.seed(23) 
LR.Boot_train.control <- trainControl(
  method = "boot632", # Optimistic corrected bootstrap
  number = 1000, # Increased
  classProbs = TRUE, 
  summaryFunction = twoClassSummary, 
  savePredictions = "final")
## Tuning Grid ##
LR.Boot_lassoGrid <- expand.grid(
  alpha = 1,
  lambda = 10^seq(-3, -0.5, length = 8))

# Train the model
set.seed(23)
log2tpm_LR.Boot_model <- caret::train(
  x = model_X,
  y = model_Y,
  method = "glmnet",
  metric = "ROC",
  tuneGrid = LR.Boot_lassoGrid,
  preProcess = c("nzv", "center", "scale"),
  trControl = LR.Boot_train.control)
print(log2tpm_LR.Boot_model)
# The final values used for the model were alpha = 1 and lambda = 0.005179475.

# Collect the genes used
log2tpm_LR.Boot_coef_df <- coef(log2tpm_LR.Boot_model$finalModel, s = log2tpm_LR.Boot_model$bestTune$lambda) %>%
  as.matrix() %>%
  as.data.frame()
log2tpm_LR.Boot_coef_df <- log2tpm_LR.Boot_coef_df[log2tpm_LR.Boot_coef_df[, 1] != 0, , drop = FALSE]
colnames(log2tpm_LR.Boot_coef_df)[1] <- "V1"
log2tpm_LR.Boot_Genes <- rownames(log2tpm_LR.Boot_coef_df)
log2tpm_LR.Boot_Genes
# [1] "(Intercept)" "Rv0023"      "Rv0105c"     "Rv0407"      "Rv0432"      "Rv0434"     
# [7] "Rv0446c"     "Rv0474"      "Rv0689c"     "Rv0720"      "Rv0837c"     "Rv0876c"    
# [13] "Rv0881"      "Rv0892"      "Rv0927c"     "Rv1028A"     "Rv1031"      "Rv1112"     
# [19] "Rv1202"      "Rv1328"      "Rv1440"      "Rv1482c"     "Rv1516c"     "Rv1538c"    
# [25] "Rv1647"      "Rv1652"      "Rv1835c"     "Rv1840c"     "Rv1896c"     "Rv1913"     
# [31] "Rv2107"      "Rv2108"      "Rv2232"      "Rv2240c"     "Rv2420c"     "Rv2699c"    
# [37] "Rv2802c"     "Rv2830c"     "Rv2873"      "Rv2943A"     "Rv2966c"     "Rv3082c"    
# [43] "Rv3123"      "Rv3208A"     "Rv3221c"     "Rv3332"      "Rv3334"      "Rv3351c"    
# [49] "Rv3389c"     "Rv3391"      "Rv3405c"     "Rv3495c"     "Rv3555c"     "Rv3650"     
# [55] "Rv3775" 

# Predict case probability using final fitted lasso model (on the same dataset)
log2tpm_LR.Boot_final_predictions <- predict(log2tpm_LR.Boot_model, newdata = model_X, type = "prob")

# Combine predictions with outcomes
log2tpm_LR.Boot_final_predictions_df <- data.frame(
  SampleID2 = rownames(model_df),
  Outcome = model_Y,
  case_probability = log2tpm_LR.Boot_final_predictions$case) %>%
  mutate(Predicted = ifelse(case_probability >= 0.5, "case", "control")) %>%
  mutate(Predicted = factor(Predicted, levels = c("case", "control"))) %>%
  mutate(Outcome2 = ifelse(Outcome == "case", 1, 0))

# Plot the model ROC
log2tpm_LR.Boot_ROC <- roc(response = log2tpm_LR.Boot_model$pred$obs,
                   predictor = log2tpm_LR.Boot_model$pred$case,
                   levels = c("control", "case"),  # (Control, case)
                   direction = "<")
plot(log2tpm_LR.Boot_ROC, print.auc = TRUE)
auc(log2tpm_LR.Boot_ROC)
# Area under the curve: 0.3738


###########################################################
################### GENESELECT R MODEL ####################
library(GeneSelectR2)

model_X_matrix <- model_X %>% as.matrix()

log2tpm_GeneSelect_model <- geneselectr2_fit(
  X = model_X_matrix,
  y = model_Y,
  bio_mode = "none",
  regularization_method = "elastic_net",
  score_weights = c(3, 3, 1), # (Stability, Prediction, Biology)
  alpha = 0.5,
  K = 4,
  R = 30,
  random_seed = 42,
  verbose = TRUE)

# Choose genes with high stability (pi > 0.5)
log2tpm_GeneSelect_GeneScores <- log2tpm_GeneSelect_model$gene_scores
log2tpm_GeneSelect_Genes <- log2tpm_GeneSelect_GeneScores$gene[log2tpm_GeneSelect_GeneScores$pi_exact > 0.5]
cat(sprintf("Stable signature (pi > 0.5): %d genes\n", length(log2tpm_GeneSelect_Genes)))
# Stable signature (pi > 0.5): 17 genes
log2tpm_GeneSelect_Genes
# [1] "Rv1440" "Rv1062"


# Build the model with just these genes (not done automatically)
log2tpm_X_SelectedGenes <- model_X_matrix[, log2tpm_GeneSelect_Genes]
set.seed(42)
GeneSelect_train.control <- trainControl(
  method = "repeatedcv",
  number = 4,
  repeats = 30,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE)
GeneSelect_tune_grid <- expand.grid(
  alpha = seq(0,1,by=0.1), # Elastic net
  lambda = 10^seq(-3,0, length=50))
log2tpm_GeneSelect_Finalmodel <- train(
  x = log2tpm_X_SelectedGenes,
  y = model_Y,
  method = "glmnet",
  metric = "ROC",
  tuneGrid = GeneSelect_tune_grid,
  trControl = GeneSelect_train.control,
  preProcess = c("center","scale"))
print(log2tpm_GeneSelect_Finalmodel)

# Predict case probability using final fitted model (on the same dataset)
log2tpm_GeneSelect_final_predictions <- predict(log2tpm_GeneSelect_Finalmodel, newdata = model_X, type = "prob")

# Combine predictions with outcomes
log2tpm_GeneSelect_final_predictions_df <- data.frame(
  SampleID2 = rownames(model_df),
  Outcome = model_Y,
  case_probability = log2tpm_GeneSelect_final_predictions$case) %>%
  mutate(Predicted = ifelse(case_probability >= 0.5, "case", "control")) %>%
  mutate(Predicted = factor(Predicted, levels = c("case", "control"))) %>%
  mutate(Outcome2 = ifelse(Outcome == "case", 1, 0))

# Plot the model ROC
log2tpm_GeneSelect_ROC <- roc(response = log2tpm_GeneSelect_Finalmodel$pred$obs,
                      predictor = log2tpm_GeneSelect_Finalmodel$pred$case,
                      levels = c("control", "case"),  # (Control, case)
                      direction = "<")
plot(log2tpm_GeneSelect_ROC, print.auc = TRUE)
auc(log2tpm_GeneSelect_ROC)
# Area under the curve: 0.6777

###########################################################
##################### SWITCHBOX MODEL #####################
# Top scoring pairs
library(switchBox)
library(logistf)

# Helper function
make_pair_features <- function(expr, pairs){
  ## Converts gene pairs into binary predictors ##
  out <- matrix(0, nrow = nrow(expr), ncol = nrow(pairs))
  colnames(out) <- apply(pairs,1,paste,collapse="_")
  rownames(out) <- rownames(expr)
  for(i in seq_len(nrow(pairs))){
    g1 <- pairs[i,1]
    g2 <- pairs[i,2]
    out[,i] <- as.integer(expr[,g1] > expr[,g2])}
  as.data.frame(out)}

# Feature selection on All samples
set.seed(42)
log2tpm_feature_folds <- createMultiFolds(model_Y, k=4, times=30)
log2tpm_pair_counts <- list()
for(i in seq_along(log2tpm_feature_folds)){
  cat("Feature fold:", i,"\n")
  train <- log2tpm_feature_folds[[i]]
  test <- setdiff(seq_len(nrow(model_X)), train)
  train_expr <- t(as.matrix(model_X[train,]))
  test_expr <- t(as.matrix(model_X[test,]))
  
  fit <- SWAP.GetKTSP.TrainTestResults(
    trainMat=train_expr,
    trainGroup=model_Y[train],
    testMat=test_expr,
    testGroup=model_Y[test],
    predictions=TRUE,
    decision_values=TRUE)
  
  selected_pairs <- rownames(fit$classifier$TSPs)
  for(p in selected_pairs){
    if(is.null(log2tpm_pair_counts[[p]])){
      log2tpm_pair_counts[[p]] <- 0}
    log2tpm_pair_counts[[p]] <- log2tpm_pair_counts[[p]] + 1}
} # End

# Select Final Stable Pairs
log2tpm_SwitchBox_pair_stability <- data.frame(Pair = names(log2tpm_pair_counts),
                                       Count = base::unlist(log2tpm_pair_counts)) %>%
  mutate(Frequency = Count / length(log2tpm_feature_folds)) %>%
  arrange(desc(Frequency))
print(head(log2tpm_SwitchBox_pair_stability,20))

# Choose final number of pairs
# n_final_pairs <- 5
# final_pairs <- pair_stability %>%
# slice_head(n=n_final_pairs)
log2tpm_SwitchBox_final_pairs <- log2tpm_SwitchBox_pair_stability %>% 
  filter(Count >= 10) 
log2tpm_SwitchBox_final_pairs <- log2tpm_SwitchBox_final_pairs %>%
  mutate(Pair_sorted = sapply(strsplit(Pair,","), 
                              function(x) paste(sort(x), collapse=","))) %>%
  distinct(Pair_sorted, .keep_all=TRUE)
print(log2tpm_SwitchBox_final_pairs)
# Pair Count  Frequency     Pair_sorted
# Rv3110,Rv2801A   Rv3110,Rv2801A    39 0.32500000  Rv2801A,Rv3110
# Rv3574,Rv2598     Rv3574,Rv2598    24 0.20000000   Rv2598,Rv3574
# Rv1516c,Rv0600c Rv1516c,Rv0600c    22 0.18333333 Rv0600c,Rv1516c
# Rv1476,Rv2341     Rv1476,Rv2341    16 0.13333333   Rv1476,Rv2341
# Rv3030,Rv2593c   Rv3030,Rv2593c    16 0.13333333  Rv2593c,Rv3030
# Rv1500,Rv2999     Rv1500,Rv2999    13 0.10833333   Rv1500,Rv2999
# Rv3792,Rv1330c   Rv3792,Rv1330c    12 0.10000000  Rv1330c,Rv3792
# Rv0913c,Rv3226c Rv0913c,Rv3226c    11 0.09166667 Rv0913c,Rv3226c
# Rv0364,Rv0963c   Rv0364,Rv0963c    10 0.08333333  Rv0364,Rv0963c

# Grab the unique genes in all these pairs
log2tpm_SwitchBox_Genes <- log2tpm_SwitchBox_final_pairs %>%
  separate_rows(Pair, sep = ",") %>%
  distinct(Pair) %>%
  pull(Pair) %>% sort()
log2tpm_SwitchBox_Genes
# [1] "Rv0364"  "Rv0600c" "Rv0913c" "Rv0963c" "Rv1330c" "Rv1476"  "Rv1500"  "Rv1516c" "Rv2341" 
# [10] "Rv2593c" "Rv2598"  "Rv2801A" "Rv2999"  "Rv3030"  "Rv3110"  "Rv3226c" "Rv3574"  "Rv3792" 


# Organize pairs for building the model
log2tpm_SwitchBox_final_pair_matrix <- do.call(rbind, strsplit(log2tpm_SwitchBox_final_pairs$Pair,","))

# Train a final classifier
log2tpm_SwitchBox_FinalTrain <- SWAP.Train.KTSP(
  inputMat = t(as.matrix(model_X)),
  phenoGroup = model_Y,
  RestrictedPairs = log2tpm_SwitchBox_final_pair_matrix)

# Predict based on the final classifier
log2tpm_SwitchBox_FinalModel <- SWAP.KTSP.Classify(
  inputMat = t(as.matrix(model_X)),
  classifier = log2tpm_SwitchBox_FinalTrain)

# Try to see what other things I can find with this package
log2tpm_SwitchBox_FinalStats <- SWAP.KTSP.Statistics(
  inputMat = t(as.matrix(model_X)),
  classifier = log2tpm_SwitchBox_FinalTrain)
heatmap(1*log2tpm_SwitchBox_FinalStats$comparisons, scale="none", labRow= rownames(model_X))

log2tpm_SwitchBox_Results <- SWAP.GetKTSP.Result(
  inputMat = t(as.matrix(model_X)),
  Groups = model_Y,
  classifier = log2tpm_SwitchBox_FinalTrain,
  predictions=TRUE, decision_values=TRUE)
log2tpm_SwitchBox_Results$roc
# Area under the curve: 0.9965

###########################################################
#################### COLLECT FINAL GENES ##################

# Get the genes that show up in all the models
log2tpm_ConsensusGenes <- Reduce(intersect, list(log2tpm_LR_Genes, log2tpm_LR.Boot_Genes, log2tpm_GeneSelect_Genes, log2tpm_SwitchBox_Genes))
log2tpm_ConsensusGenes
# character(0)

# Get the genes that show up in 3/4 models
log2tpm_geneLists <- list(log2tpm_LR = log2tpm_LR_Genes, log2tpm_LR.Boot = log2tpm_LR.Boot_Genes, log2tpm_GeneSelect = log2tpm_GeneSelect_Genes, log2tpm_SwitchBox = log2tpm_SwitchBox_Genes)
log2tpm_allGenes <- unlist(log2tpm_geneLists)
log2tpm_geneCounts <- table(log2tpm_allGenes)
ConsensusGenes_3plus <- names(log2tpm_geneCounts[log2tpm_geneCounts >= 3])
ConsensusGenes_3plus
# "Rv1440" 

# Get the genes that show up in 2/3 models (removign the LR.Boot because it is really similar to the LR CV)
log2tpm_geneLists2 <- list(log2tpm_LR = log2tpm_LR_Genes, log2tpm_GeneSelect = log2tpm_GeneSelect_Genes, log2tpm_SwitchBox = log2tpm_SwitchBox_Genes)
log2tpm_allGenes2 <- unlist(log2tpm_geneLists2)
log2tpm_geneCounts2 <- table(log2tpm_allGenes2)
ConsensusGenes_2plus <- names(log2tpm_geneCounts2[log2tpm_geneCounts2 >= 2])
ConsensusGenes_2plus
# [1] "Rv1062" "Rv1440"

# How correlated are my genes?
cor(model_df[,c("Rv1062", "Rv1440")], use="complete.obs")
# Rv1062    Rv1440
# Rv1062 1.0000000 0.5948229
# Rv1440 0.5948229 1.0000000

###########################################################
################# ORGANIZE DATA w/METADATA ################

model_df2 <- model_df %>% dplyr::select(Outcome, all_of(ConsensusGenes_2plus)) %>%
  mutate(Outcome = factor(Outcome, levels = c("case", "control")))

model_X2 <- model_df2 %>% dplyr::select(-Outcome)
model_Y2 <- model_df2$Outcome

# Change outcome to 0 and 1
model_df2 <- model_df2 %>%
  mutate(Outcome2 = ifelse(Outcome == "control", 0, 1))


###########################################################
######################### 2 GENES #########################
# 9/24/26

# Firth logistic regression
Genes_firth <- logistf(Outcome2 ~ Rv1062 + Rv1440, data = model_df2)
Genes_firth
summary(Genes_firth)
# logistf(formula = Outcome2 ~ Rv1062 + Rv1440, data = model_df2)
# 
# Model fitted by Penalized ML
# Coefficients:
#   coef     se(coef)    lower 0.95    upper 0.95    Chisq         p method
# (Intercept) -2.8270673397 4.950354e-01 -3.978293e+00 -1.9260229446      Inf 0.0000000      2
# Rv1062       0.0001488819 7.957071e-05 -2.304183e-05  0.0003139841 2.511278 0.1130342      2
# Rv1440       0.0003314159 1.395803e-04  5.062652e-05  0.0006529757 5.190231 0.0227142      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=17.18378 on 2 df, p=0.000185605, n=88
# Wald test = 34.20496 on 2 df, p = 3.7367e-08

Genes_firth_validate <- validate_firth(Outcome2 ~ Rv1062 + Rv1440, data = model_df2, outcome = "Outcome2", B = 2000)
Genes_firth_validate
# Metric    Apparent     Optimism   Corrected
# 1                   AUC  0.70484061  0.020294090  0.68454652
# 2                   Dxy  0.40968123  0.040588181  0.36909305
# 3                 Brier  0.07294944 -0.008184632  0.08113407
# 4 Intercept.(Intercept) -0.00693627  0.072161473 -0.07909774
# 5              Slope.lp  1.00550884  0.064780860  0.94072798
attr(Genes_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Genes_firth_df <- model_df2
Genes_firth_df$Predicted_prob <- predict(Genes_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Genes_firth_calibrate <- calibrate_firth(Outcome2 ~ Rv1062 + Rv1440, data = model_df2, B = 2000)
Genes_firth_calibrate
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Rv2385_v1.pdf", width = 5, height = 4)
plot(Genes_firth_calibrate, main = "Rv1062 + Rv1440")
# dev.off()

# ROC - Can the model rank patients correctly?
Genes_firth_ROC <- roc(response = Genes_firth_df$Outcome2, predictor = Genes_firth_df$Predicted_prob)
plot(Genes_firth_ROC, print.auc = TRUE, col = "blue")
auc(Genes_firth_ROC)
# Area under the curve: 0.7048
set.seed(42)
ci.auc(Genes_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 0.4652-0.9126 (2000 stratified bootstrap replicates)







###########################################################
################ DECISION CURVES TOGETHER #################

# 8/14/26: NOT YET DONE FOR THIS SET

# base_firth_DC, base.Arm_DC_firth, Base.Arm.Rv3747.Rv0161.Rv2075c_firth_DC, Rv1016_firth_DC, Rv2075_firth_DC, Rv2385_firth_DC, Rv3747_firth_DC. Rv3747.Rv2075_firth_DC, Rv3747.Rv2385_firth_DC, Rv3747.Rv0161.Rv2075c_firth_DC

plot_decision_curve(
  legend.position = "topright",
  xlim = c(0,1),
  list(base_firth_DC,
       base.Arm_DC_firth),
  curve.names = c("Base",
                  "Base + Arm"))

# pdf("Figures/PredictiveModeling/W2/DecisionCurves/DC_v1.pdf", width = 6, height = 5)
plot_decision_curve(
  legend.position = "none",
  xlim = c(0,0.6),
  confidence.intervals = FALSE,
  list(base_firth_DC,
       base.Arm_DC_firth,
       Base.Arm.Rv3747.Rv0161.Rv2075c_firth_DC,
       Rv3747.Rv0161.Rv2075c_firth_DC,
       Arm.Rv3747.Rv0161.Rv2075c_firth_DC),
  curve.names = c("Base",
                  "Base + Arm",
                  "Base + Arm + Rv3747+Rv2075c+Rv0161",
                  "Rv3747 + Rv2075c + Rv0161",
                  "Arm + Rv3747+Rv2075c+Rv0161"),
  col = c("#000000", "#8B008B", "#8C564B", "#0072B2", "orange"))
legend("bottomleft",
       legend = c("Base",
                  "Base + Arm",
                  "Base + Arm + Rv3747+Rv2075c+Rv0161",
                  "Rv3747 + Rv2075c + Rv0161",
                  "Treat all", 
                  "Treat none"),
       col = c("#000000", "#8B008B", "#8C564B", "#0072B2", "grey", "black"),
       lwd = 2, cex = 0.6, bg = "white", bty = "o")
# dev.off()

# pdf("Figures/PredictiveModeling/W2/DecisionCurves/DC_v2.pdf", width = 6, height = 5)
plot_decision_curve(
  legend.position = "none",
  xlim = c(0,0.6),
  confidence.intervals = FALSE,
  list(Rv3747.Rv0161.Rv2075c_firth_DC,
       Rv3747_firth_DC,
       Rv0161_firth_DC,
       Rv2075c_firth_DC,
       Rv2385_firth_DC),
  col = c("#0072B2", "#D55E00", "#009E73", "yellow", "#CC79A7"),
  curve.names = c("3 Genes",
                  "Rv3747",
                  "Rv0161",
                  "Rv2075c",
                  "Rv2385c"))
legend("bottomleft",
       legend = c("3 Genes",
                  "Rv3747",
                  "Rv0161",
                  "Rv2075c",
                  "Rv2385c"),
       col = c("#0072B2", "#D55E00", "#009E73", "#CC79A7","yellow", "grey", "black"),
       lwd = 2, cex = 0.6, bg = "white", bty = "o")
# dev.off()


###########################################################
################### ROC CURVES TOGETHER ###################

my_fav_colors <- c(`Base` = "#000000", 
                   `Base.Arm` = "#8B008B", 
                   `Base.Arm.Rv3747.Rv0161.Rv2075c`= "#8C564B", 
                   `Rv3747.Rv0161.Rv2075c` = "#0072B2", 
                   `Rv3747.Rv1138c.Rv2075c` = "#87CEEB", 
                   `Rv3747.Rv0161.Rv1138c.Rv2075c` = "darkblue",
                   `Rv3747` = "#D55E00", 
                   `Rv0161` = "#009E73",
                   `Rv2075c` = "#CC79A7",
                   `Rv1138c` = "#F0E442")

# pdf("Figures/PredictiveModeling/W2/ROC/ROC_v1.pdf", width = 6, height = 5)
plot(base_firth_ROC, col = my_fav_colors[1], lwd = 2, lty = 1, 
     print.auc = FALSE)
lines(base.Arm_firth_ROC, col = my_fav_colors[2], lty = 1)
lines(Base.Arm.Rv3747.Rv0161.Rv2075c_firth_ROC, col = my_fav_colors[3], lty = 2)
lines(Rv3747.Rv0161.Rv2075c_firth_ROC, col = my_fav_colors[4], lty = 3)
lines(Rv3747.Rv1138c.Rv2075c_firth_ROC, col = my_fav_colors[5], lty = 4)
# lines(Rv3747_firth_ROC, col = my_fav_colors[6], lty = 5)
# lines(Rv0161_firth_ROC, col = my_fav_colors[7], lty = 5)
# lines(Rv2075c_firth_ROC, col = my_fav_colors[8], lty = 5)
# lines(Rv1138c_firth_ROC, col = my_fav_colors[9], lty = 5)
legend("bottomright", 
       legend = c(sprintf("Base (AUC = %.3f)", auc(base_firth_ROC)),
                  sprintf("Base+Arm (AUC = %.3f)", auc(base.Arm_firth_ROC)),
                  sprintf("Base+Arm+Rv3747+Rv0161+Rv2075c (AUC = %.3f)", auc(Base.Arm.Rv3747.Rv0161.Rv2075c_firth_ROC)),
                  sprintf("Rv3747+Rv0161.Rv2075c (AUC = %.3f)", auc(Rv3747.Rv0161.Rv2075c_firth_ROC)),
                  sprintf("Rv3747+Rv1138c+Rv2075c (AUC = %.3f)", auc(Rv3747.Rv1138c.Rv2075c_firth_ROC))),
       col = my_fav_colors,
       lty = c(1,1,2,3,4),
       lwd = 2, cex = 0.7, bty = "n")
# dev.off()



###########################################################
################### ANOVA COMPARE MODELS ##################

# base_firth, base.Arm_firth, Base.Arm.Rv3747.Rv0161.Rv2075c_firth, Rv1016_firth, Rv2075_firth, Rv2385_firth, Rv3747_firth, Rv3747.Rv2075_firth, Rv3747.Rv2385_firth, Rv3747.Rv0161.Rv2075c_firth

# Using a Penalized likelihood ratio test
# Does the larger model significantly improve fit compared with the smaller model?
# Hierachical nested models should use "nested", others should use "PLR"

anova(base.Arm_firth, Base.Arm.Rv3747.Rv0161.Rv2075c_firth, method = "nested")
# Comparison of logistf models:
#   Formula ChiSquared 
# 1 Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age + Arm2 + Arm3 + Rv3747 + Rv0161 + Rv2075c  23.037016 
# 2                             Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age + Arm2 + Arm3   7.622569 
# 
# Method:  nested 
# Chi-Squared:  15.41445   df= 3   P= 0.00149464 

# "Addition of the three genes significantly improved model fit compared with the clinical model alone"


anova(Rv3747.Rv0161.Rv2075c_firth, Base.Arm.Rv3747.Rv0161.Rv2075c_firth, method = "nested")
# Comparison of logistf models:
#   Formula ChiSquared 
# 1 Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age + Arm2 + Arm3 + Rv3747 + Rv0161 + Rv2075c   23.03702 
# 2                                               Outcome2 ~ Rv3747 + Rv0161 + Rv2075c   20.95215 
# 
# Method:  nested 
# Chi-Squared:  2.08487   df= 6   P= 0.9117298 

# Adding the clinical variables to the 3-gene model did not significantly improve model fit.


anova(Rv3747_firth, Rv3747.Rv0161.Rv2075c_firth, method = "nested")
# Comparison of logistf models:
#   Formula ChiSquared 
# 1 Outcome2 ~ Rv3747 + Rv0161 + Rv2075c   27.20335 
# 2                    Outcome2 ~ Rv3747   19.96357 
# 
# Method:  nested 
# Chi-Squared:  7.239785   df= 2   P= 0.02678556 

# "The addition of 2 genes to the Rv3747 model improved model fit"

anova(Rv2075c_firth, Rv3747.Rv0161.Rv2075c_firth, method = "nested")
# Comparison of logistf models:
#   Formula ChiSquared 
# 1 Outcome2 ~ Rv3747 + Rv0161 + Rv2075c   27.20335 
# 2                   Outcome2 ~ Rv2075c   15.34056 
# 
# Method:  nested 
# Chi-Squared:  11.8628   df= 2   P= 0.002654765 


anova(Rv0161_firth, Rv3747.Rv0161.Rv2075c_firth, method = "nested")
# Comparison of logistf models:
#   Formula ChiSquared 
# 1 Outcome2 ~ Rv3747 + Rv0161 + Rv2075c   27.20335 
# 2                    Outcome2 ~ Rv0161   18.08080 
# 
# Method:  nested 
# Chi-Squared:  9.122554   df= 2   P= 0.01044871 

anova(Rv3747.Rv0161.Rv2075c_firth, Rv3747.Rv0161.Rv1138c.Rv2075c_firth, method = "nested")
# Comparison of logistf models:
#   Formula ChiSquared 
# 1 Outcome2 ~ Rv3747 + Rv0161 + Rv1138c + Rv2075c   26.62719 
# 2           Outcome2 ~ Rv3747 + Rv0161 + Rv2075c   25.62074 
# 
# Method:  nested 
# Chi-Squared:  1.006444   df= 1   P= 0.3157563 

# Adding Rv1138c doesn't actually help...

###########################################################
################### TPR and FPR TOGETHER ##################

# base_firth_DC, base.Arm_DC_firth, Base.Arm.Rv3747.Rv0161.Rv2075c_firth_DC, Rv1016_firth_DC, Rv2075_firth_DC, Rv2385_firth_DC, Rv3747_firth_DC. Rv3747.Rv2075_firth_DC, Rv3747.Rv2385_firth_DC, Rv3747.Rv0161.Rv2075c_firth_DC, Rv1138c_firth_DC

my_fav_colors <- c(`Base` = "#000000", 
                   `Base.Arm` = "#8B008B", 
                   `Base.Arm.Rv3747.Rv0161.Rv2075c`= "#8C564B", 
                   `Rv3747.Rv0161.Rv2075c` = "#0072B2", 
                   `Rv3747.Rv1138c.Rv2075c` = "#87CEEB", 
                   `Rv3747.Rv0161.Rv1138c.Rv2075c` = "darkblue",
                   `Rv3747` = "#D55E00", 
                   `Rv0161` = "#009E73",
                   `Rv2075c` = "#CC79A7",
                   `Rv1138c` = "#F0E442")

# Plot TPR and FPR (Alternative to ROC)
# Input is the decision curve output (which is the adjusted probabilities)
# Want the FPR to fall faster than the TPR at different risk thresholds
plot_roc_components(base.Arm_DC_firth, col = my_fav_colors[2], confidence.intervals = FALSE)
plot_roc_components(Base.Arm.Rv3747.Rv0161.Rv2075c_firth_DC, col = my_fav_colors[3], 
                    confidence.intervals = FALSE, add = TRUE)
plot_roc_components(Rv3747.Rv0161.Rv2075c_firth_DC, col = "#0072B2", confidence.intervals = FALSE)
plot_roc_components(Rv3747_firth_DC, col = "#D55E00", confidence.intervals = FALSE)
plot_roc_components(Rv0161_firth_DC, col = "#009E73", confidence.intervals = FALSE)
plot_roc_components(Rv2075c_firth_DC, col = "#CC79A7", confidence.intervals = FALSE)
plot_roc_components(Rv2385_firth_DC, col = "yellow", confidence.intervals = FALSE)
plot_roc_components(Rv1138c_firth_DC, col = "green", confidence.intervals = FALSE)

# Generate plot with all together
# First have to make a function to remove the ends of the data or there will be wierd lines on the graph
plot_component <- function(dc_object, col) {
  dat <- dc_object$derived.data
  # remove reference strategies
  dat <- dat[!dat$model %in% c("All", "None"), ]
  lines(dat$thresholds, dat$TPR, col = col, lwd = 2)
  lines(dat$thresholds, dat$FPR, col = col, lwd = 2, lty = 2)
}

# first model
dat <- base_firth_DC$derived.data
dat <- dat[!dat$model %in% c("All", "None"), ]

# pdf("Figures/PredictiveModeling/W2/TPR_FPR/TPR.FPR_v1.pdf", width = 7, height = 5)
plot(dat$thresholds, dat$TPR,
     type = "l", col = my_fav_colors[1], lwd = 2,
     ylim = c(0,1), xlim = c(0,1),
     xlab = "Risk Threshold", ylab = "Probability", main = "TPR and FPR")
lines(dat$thresholds, dat$FPR, col = my_fav_colors[1], lwd = 2, lty = 2)
plot_component(base.Arm_DC_firth, my_fav_colors[2])
plot_component(Base.Arm.Rv3747.Rv0161.Rv2075c_firth_DC, my_fav_colors[3])
plot_component(Rv3747.Rv0161.Rv2075c_firth_DC, my_fav_colors[4])
plot_component(Rv3747.Rv1138c.Rv2075c_firth_DC, my_fav_colors[5])
# plot_component(Rv3747_firth_DC, "#D55E00")
# plot_component(Rv0161_firth_DC, "#009E73")
# plot_component(Rv2075c_firth_DC, "#CC79A7")
# plot_component(Rv2385_firth_DC, "yellow")
# plot_component(Rv1138c_firth_DC, "green")
legend("topright",
       legend=c("Base",
                "Base+Arm",
                "Base+Arm+Rv3747+Rv0161+Rv2075c",
                "Rv3747+Rv0161+Rv2075c",
                "Rv3747+Rv1138c+Rv2075c"),
       col=c(my_fav_colors[1], my_fav_colors[2], my_fav_colors[3], my_fav_colors[4], my_fav_colors[5]), 
       lwd=2, bty = "n", cex = 0.6)
legend("top", legend=c("TPR","FPR"),
       lty=c(1,2), lwd=2, bty = "n", cex = 0.6)
# dev.off()

###########################################################
############### PROBABILITY BOXPLOTS TOGETHER #############

# base_firth_df, base.Arm_firth_df, Base.Arm.Rv3747.Rv0161.Rv2075c_firth_df, Rv1016_firth_df, Rv2075_firth_df, Rv2385_firth_df, Rv3747_firth_df, Rv3747.Rv2075_firth_df, Rv3747.Rv2385_firth_df, Rv3747.Rv0161.Rv2075c_firth_df, Rv1138c_firth_df

# Plot basics
my_plot_themes <- theme_bw() +
  # theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none", legend.text=element_text(size=10),
        legend.title = element_blank(),
        plot.title = element_text(size=10), 
        axis.title.x = element_blank(), 
        # axis.text.x = element_text(angle = 0, size=10, vjust=0, hjust=0.5),
        axis.text.x = element_text(angle = 45, size=10, vjust=1, hjust=1),
        axis.title.y = element_text(size=10),
        axis.text.y = element_text(size=10), 
        plot.subtitle = element_text(size=9), 
        plot.margin = margin(2, 2, 2, 2))
facet_themes <- theme(strip.background=element_rect(fill="#ededed", linewidth = 0.9),
                      strip.text = element_text(size = 10))

my_fav_colors <- c(`Base` = "#000000", 
                   `Base.Arm` = "#8B008B", 
                   `Base.Arm.Rv3747.Rv0161.Rv2075c`= "#8C564B", 
                   `Rv3747.Rv0161.Rv2075c` = "#0072B2", 
                   `Rv3747.Rv1138c.Rv2075c` = "#87CEEB", 
                   `Rv3747.Rv0161.Rv1138c.Rv2075c` = "darkblue",
                   `Rv3747` = "#D55E00", 
                   `Rv0161` = "#009E73",
                   `Rv2075c` = "#CC79A7",
                   `Rv1138c` = "#F0E442")

All_ApparentPredictedProbabilities_df <- base_firth_df %>% 
  rownames_to_column("SampleID2") %>%
  select(SampleID2, Outcome, Adjusted_prob) %>%
  dplyr::rename(Base = Adjusted_prob) %>%
  mutate(Base.Arm = base.Arm_firth_df$base.Arm_firth_popprob) %>%
  mutate(Base.Arm.Rv3747.Rv0161.Rv2075c = Base.Arm.Rv3747.Rv0161.Rv2075c_firth_df$Adjusted_prob) %>%
  mutate(Rv3747.Rv0161.Rv2075c = Rv3747.Rv0161.Rv2075c_firth_df$Adjusted_prob) %>%
  mutate(Rv3747.Rv1138c.Rv2075c = Rv3747.Rv1138c.Rv2075c_firth_df$Adjusted_prob) %>%
  # mutate(Rv3747.Rv0161.Rv1138c.Rv2075c = Rv3747.Rv0161.Rv1138c.Rv2075c_firth_df$Adjusted_prob) %>%
  mutate(Rv3747 = Rv3747_firth_df$Adjusted_prob) %>%
  mutate(Rv0161 = Rv0161_firth_df$Adjusted_prob) %>%
  mutate(Rv2075c = Rv2075c_firth_df$Adjusted_prob) %>%
  # mutate(Rv2385 = Rv2385_firth_df$Adjusted_prob) %>%
  mutate(Rv1138c = Rv1138c_firth_df$Adjusted_prob) %>%
  # pivot_longer(cols = Base:Rv1138c, names_to = "Model", values_to = "PredictedProbs") %>%
  pivot_longer(cols = Base:Rv1138c, names_to = "Model", values_to = "PredictedProbs") %>%
  mutate(Outcome = factor(Outcome, levels = c("control", "case"))) %>%
  mutate(Model = factor(Model, levels = c("Base", "Base.Arm", "Base.Arm.Rv3747.Rv0161.Rv2075c", "Rv3747.Rv0161.Rv2075c", "Rv3747.Rv1138c.Rv2075c", "Rv3747.Rv0161.Rv1138c.Rv2075c", "Rv3747", "Rv0161", "Rv2075c", "Rv1138c")))


PredictedProbabilities_Boxplot <- All_ApparentPredictedProbabilities_df %>%
  filter(!Model %in% c("Rv3747", "Rv0161", "Rv2075c", "Rv1138c")) %>%
  ggplot(aes(x = Model, y = PredictedProbs, fill = Model)) +
  scale_fill_manual(values = my_fav_colors) +
  geom_boxplot(alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) + 
  # geom_text_repel(aes(label = SampleID2), size = 2) + 
  facet_wrap(~Outcome) + 
  # scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.1)) + 
  labs(y = "Predicted probability of case") +
  my_plot_themes + facet_themes
PredictedProbabilities_Boxplot
# ggsave(PredictedProbabilities_Boxplot,
#        file = "Figures/PredictiveModeling/W2/ModelProbabilities_v1.pdf",
#        width = 7, height = 5, units = "in")

PredictedProbabilities_Boxplot_SingleGenes <- All_ApparentPredictedProbabilities_df %>%
  filter(Model %in% c("Rv3747.Rv0161.Rv2075c", "Rv3747.Rv1138c.Rv2075c", "Rv3747", "Rv0161", "Rv2075c", "Rv1138c")) %>%
  ggplot(aes(x = Model, y = PredictedProbs, fill = Model)) +
  scale_fill_manual(values = my_fav_colors) +
  geom_boxplot(alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) + 
  # geom_text_repel(aes(label = SampleID2), size = 2) + 
  facet_wrap(~Outcome) + 
  # scale_y_continuous(limits = c(0, 0.9), breaks = seq(0, 0.9, 0.1)) + 
  labs(y = "Predicted probability of case") +
  my_plot_themes + facet_themes
PredictedProbabilities_Boxplot_SingleGenes
# ggsave(PredictedProbabilities_Boxplot_SingleGenes,
#        file = "Figures/PredictiveModeling/W2/ModelProbabilities_SingleGenes_v1.pdf",
#        width = 7, height = 5, units = "in")




