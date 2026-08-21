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

###########################################################
#################### COLLECT FINAL GENES ##################

# Get the genes that show up in all the models
log2tpm_ConsensusGenes <- Reduce(intersect, list(log2tpm_LR_Genes, log2tpm_LR.Boot_Genes, log2tpm_GeneSelect_Genes, log2tpm_SwitchBox_Genes))
log2tpm_ConsensusGenes

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
# "Rv0161"  "Rv0238"  "Rv1138c" "Rv1267c" "Rv1934c" "Rv1952"  "Rv2075c" "Rv2092c" "Rv2415c" "Rv3747" 

# Adding Rv0269c becuase that goes up in case and Rv0089 because I like it
ConsensusGenes2 <- c(ConsensusGenes, "Rv0269c", "Rv0089", "Rv1955", "Rv2385")

# How correlated are my genes?
cor(model_df2[,c("Rv2385", "Rv2075c", "Rv3747", "Rv0161", "Rv1138c")], use="complete.obs")
# Rv2385    Rv2075c     Rv3747     Rv0161    Rv1138c
# Rv2385   1.0000000  0.6334566  0.5814414  0.5550559 -0.2861174
# Rv2075c  0.6334566  1.0000000  0.5738714  0.3719474 -0.4217958
# Rv3747   0.5814414  0.5738714  1.0000000  0.3790249 -0.1770300
# Rv0161   0.5550559  0.3719474  0.3790249  1.0000000 -0.4965064
# Rv1138c -0.2861174 -0.4217958 -0.1770300 -0.4965064  1.0000000

###########################################################
################# ORGANIZE DATA w/METADATA ################

model_df2 <- GoodSamples60_pipeSummary %>% 
  filter(Type == "Week 2 sputum") %>%
  dplyr::select(SampleID2, Outcome, TTD, XpertCT_wk0, BMI, Age, Arm, main_lineage) %>% # Keep current metrics/biomarkers
  mutate(XpertCT_wk0 = as.numeric(XpertCT_wk0)) %>%
  mutate(Outcome = factor(Outcome, levels = c("case", "control"))) %>%
  inner_join(tmp_log2tpm, by = "SampleID2") %>%
  column_to_rownames("SampleID2") %>%
  na.omit()

# Collapsing arm to make binary, A have worse cavitation than B or C
model_df2 <- model_df2 %>%
  mutate(Arm2 = ifelse(Arm == "A", "A", "BorC")) %>%
  mutate(Arm3 = ifelse(Arm == "C", "C", "AorB"))

# Collapsing lineage to make binary, There is only 1 lineage3, combining with lineage2 
model_df2 <- model_df2 %>%
  mutate(Lineage2 = ifelse(main_lineage == "lineage4", "lineage4", "lineage2or3"))

model_df2 <- model_df2 %>% dplyr::select(Outcome, TTD, XpertCT_wk0, BMI, Age, Arm2, Arm3, main_lineage,  all_of(ConsensusGenes2), )

model_X2 <- model_df2 %>% dplyr::select(-Outcome)
model_Y2 <- model_df2$Outcome

# Change outcome to 0 and 1
model_df2 <- model_df2 %>%
  mutate(Outcome2 = ifelse(Outcome == "control", 0, 1))


###########################################################
#################### CHECK FOR LINEARITY ##################
# 8/6/26: Read Tim Sterling's paper and he was talking about checking continuous variables for linearity so here I am
# *!*!*!* Need to think about this some more !*!*!*!*
library(rms)

# Chunk test from rms package?
dd <- datadist(model_df2)
options(datadist = "dd")

fit <- lrm(
  Outcome2 ~ rcs(Age, 3) + rcs(BMI, 3) +
    rcs(TTD, 3) + rcs(XpertCT_wk0, 3), # + 
  # rcs(Rv3747, 3) + rcs(Rv0161, 3) + rcs(Rv2075c, 3),
  data = model_df2
)

anova(fit)
# Wald Statistics          Response: Outcome2 
# 
# Factor          Chi-Square d.f. P     
# Age             0.08       2    0.9586
#   Nonlinear     0.02       1    0.8852
# BMI             3.30       2    0.1919
#   Nonlinear     2.00       1    0.1574
# TTD             4.12       2    0.1277
#   Nonlinear     4.12       1    0.0425
# XpertCT_wk0     2.29       2    0.3188
#   Nonlinear     0.39       1    0.5330
# TOTAL NONLINEAR 4.35       4    0.3605
# TOTAL           6.02       8    0.6452

plot(Predict(fit, TTD))
plot(Predict(fit, Age))


# "Linearity of continuous predictors was assessed using restricted cubic splines with three knots and Wald chunk tests. No evidence of overall nonlinearity was observed (overall nonlinear Wald test, P = 0.36). Although TTD showed evidence of a nonlinear association (P = 0.043), spline terms were not retained because of the limited sample size and absence of overall improvement in model fit."

ggplot(model_df2, aes(TTD, Outcome2)) +
  geom_jitter(height = 0.05) +
  geom_smooth(
    method = "glm",
    method.args = list(family = binomial),
    formula = y ~ splines::ns(x, df = 3)
  )



###########################################################
########################## BASE ###########################
# Updated 8/14/26

base_firth <- logistf(Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age, data = model_df2)
base_firth
summary(base_firth)
# logistf(formula = Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age, data = model_df2)
# 
# Model fitted by Penalized ML
# Coefficients:
#                     coef   se(coef)  lower 0.95 upper 0.95      Chisq         p method
# (Intercept)  1.44339408 3.22854333 -4.14922553 13.8261773 0.19657373 0.6575008      2
# TTD          0.02028666 0.09576925 -0.18959962  0.2159827 0.04378891 0.8342468      2
# XpertCT_wk0 -0.06423823 0.11427512 -0.33167420  0.1556415 0.31998080 0.5716192      2
# BMI         -0.10587740 0.15170159 -0.69919046  0.1580985 0.52300036 0.4695646      2
# Age          0.03093571 0.03380139 -0.03754642  0.1079358 0.83476810 0.3608974      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=2.333807 on 4 df, p=0.674621, n=29
# Wald test = 3.272337 on 4 df, p = 0.5133271

# Use Bootstrapping to validate
base_firth_validate <- validate_firth(Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age, data = model_df2, outcome = "Outcome2", B = 2000)
base_firth_validate
#                   Metric    Apparent    Optimism  Corrected
# 1                   AUC 0.67676768  0.1338258  0.54294192
# 2                   Dxy 0.35353535  0.2676515  0.08588384
# 3                 Brier 0.21065548 -0.0598897  0.27054517
# 4 Intercept.(Intercept) 0.02805212  0.2925095 -0.26445735
# 5              Slope.lp 1.22805330  0.7744758  0.45357752
attr(base_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
base_firth_df <- model_df2
base_firth_df$Predicted_prob <- predict(base_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
base_firth_calibrate <- calibrate_firth(Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/base_v1.pdf", width = 5, height = 4)
plot(base_firth_calibrate, main = "TTD + Xpert_W0 + BMI + Age")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(base_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
base_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
base_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                data=base_firth_df,
                                fitted.risk=TRUE, 
                                study.design="cohort", 
                                bootstraps=2000)
plot_decision_curve(base_firth_DC,  curve.names = "base_firth", lty = 1)

# ROC - Can the model rank patients correctly?
base_firth_ROC <- roc(response = base_firth_df$Outcome2, predictor = base_firth_df$Predicted_prob)
plot(base_firth_ROC, print.auc = TRUE, col = "blue")
auc(base_firth_ROC)
# Area under the curve: 0.6768
set.seed(42)
ci.auc(base_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 0.4444-0.8737 (2000 stratified bootstrap replicates)



###########################################################
######################### BASE+ARM ########################
# Updated 8/12/26

# Firth logistic regression
base.Arm_firth <- logistf(Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age + Arm2 + Arm3, data = model_df2)
base.Arm_firth
summary(base.Arm_firth)
#                     coef   se(coef)  lower 0.95  upper 0.95     Chisq       p method
# (Intercept)  2.36909136 3.17948830 -3.33424977 15.35448569 0.5853912 0.4442069      2
# TTD          0.02864642 0.11427775 -0.26773699  0.25962406 0.0504163 0.8223405      2
# XpertCT_wk0 -0.19631308 0.15217580 -0.64234949  0.08663131 1.7802404 0.1821198      2
# BMI         -0.12263327 0.14995698 -0.74815775  0.16051963 0.6988194 0.4031807      2
# Age          0.04261795 0.03536373 -0.02891248  0.13160031 1.4437394 0.2295352      2
# Arm2BorC     1.14523026 1.32084981 -1.71309620  4.15156547 0.7041422 0.4013953      2
# Arm3C        1.11927567 1.37073190 -1.71076632  4.20225715 0.6480737 0.4208022      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=8.27142 on 6 df, p=0.2188852, n=29
# Wald test = 7.166245 on 6 df, p = 0.3057463

# Get the odds ratios
exp(cbind(OR = coef(base.Arm_firth), confint(base.Arm_firth)))
#                     OR  Lower 95%    Upper 95%
# (Intercept) 10.6876766 0.03564132 4.659812e+06
# TTD          1.0290607 0.76510898 1.296443e+00
# XpertCT_wk0  0.8217549 0.52605501 1.090495e+00
# BMI          0.8845880 0.47323757 1.174121e+00
# Age          1.0435391 0.97150149 1.140652e+00
# Arm2BorC     3.1431650 0.18030666 6.353338e+01
# Arm3C        3.0626350 0.18072724 6.683702e+01

base.Arm_firth_validate <- validate_firth(Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age + Arm2 + Arm3, data = model_df2, outcome = "Outcome2", B = 2000)
base.Arm_firth_validate
#                 Metric    Apparent    Optimism  Corrected
# 1                   AUC  0.84343434  0.13012992  0.7133044
# 2                   Dxy  0.68686869  0.26025984  0.4266088
# 3                 Brier  0.15982048 -0.07917636  0.2389968
# 4 Intercept.(Intercept) -0.05205247  0.17230161 -0.2243541
# 5              Slope.lp  1.23025227  0.96290780  0.2673445
attr(base.Arm_firth_validate, "successful_bootstraps")
# 1733

# Make predictions based on the model
base.Arm_firth_df <- model_df2
base.Arm_firth_df$base.Arm_firth_prob <- predict(base.Arm_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
base.Arm_firth_calibrate <- calibrate_firth(Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age + Arm2 + Arm3, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/base.Arm_v1.pdf", width = 5, height = 4)
plot(base.Arm_firth_calibrate, main = "Base+Arm")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(base.Arm_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
base.Arm_firth_df$base.Arm_firth_popprob <- p_adj

# "Predicted probabilities were recalibrated to a 5% population case prevalence prior to decision curve analysis to account for intentional enrichment of case cases in the study cohort. Threshold performance curves therefore reflect sensitivity and false-positive rates at population-calibrated case risk thresholds."

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
base.Arm_DC_firth <- decision_curve(Outcome2 ~ base.Arm_firth_popprob, 
                                    data=base.Arm_firth_df,
                                    fitted.risk=TRUE, 
                                    study.design="cohort", 
                                    bootstraps=2000)
plot_decision_curve(base.Arm_DC_firth,  curve.names = "base.Arm firth", lty = 1)

# Plot risk distribution
base.Arm_firth_df %>% 
  ggplot(aes(x=base.Arm_firth_popprob, fill=factor(Outcome))) + # Can also plot base.Arm_firth_prob
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c(case = "#bc5300", control = "#0072B2")) + 
  labs(x = "Predicted case risk (5% population prevalence-adjusted)",
       y = "Density",
       title = "Distribution of predicted case risk") +
  theme_classic()

# Discrimination measures: Quantify how well does the risk model discriminate/separate cases and controls? MRD, AARD, AUC

# Above Average Risk Difference: mean difference between observed and predicted risk among patients whose predicted risk is above the average predicted risk
above_avg <- base.Arm_firth_df$base.Arm_firth_prob > mean(base.Arm_firth_df$base.Arm_firth_prob) # Define above-average risk patients
above_avg_risk_difference <- mean(model_df2$Outcome2[above_avg]) - mean(base.Arm_firth_df$base.Arm_firth_prob[above_avg])
above_avg_risk_difference # 0.01768879

# ROC - Can the model rank patients correctly?
base.Arm_firth_ROC <- roc(response = base.Arm_firth_df$Outcome2, predictor = base.Arm_firth_df$base.Arm_firth_prob)
plot(base.Arm_firth_ROC, print.auc = TRUE, col = "blue")
auc(base.Arm_firth_ROC)
# Area under the curve: 0.8434
set.seed(42)
ci.auc(base.Arm_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 0.6768-0.9596 (2000 stratified bootstrap replicates)
## *** IS THIS THE RIGHT WAY TO GET CONFIDENCE INTERVALS??


# Plot TPR and FPR (Alternative to ROC)
# Input is the decision curve output (which is the adjusted probabilities)
# Want the FPR to fall faster than the TPR at different risk thresholds
plot_roc_components(base.Arm_DC_firth, col = c("black", "red"))

# Clinical Impact Plot - not sure what this means
plot_clinical_impact(base.Arm_DC_firth, xlim = c(0, .20), col = c("black", "blue"))

# Summary graphs
summary(base.Arm_firth_df$base.Arm_firth_prob)
summary(base.Arm_firth_df$base.Arm_firth_popprob)
hist(base.Arm_firth_df$base.Arm_firth_prob,
     breaks=10,
     main="Original Firth probabilities",
     xlab="Predicted case probability")
hist(base.Arm_firth_df$base.Arm_firth_popprob,
     breaks=10,
     main="Population-adjusted probabilities",
     xlab="Predicted case probability")
summary(base.Arm_firth_df$base.Arm_firth_prob)
hist(base.Arm_firth_df$base.Arm_firth_prob)

###########################################################
########################## Rv2385 #########################
# Added 8/12/26

# Firth logistic regression
Rv2385_firth <- logistf(Outcome2 ~ Rv2385, data = model_df2)
Rv2385_firth
summary(Rv2385_firth)
# coef se(coef) lower 0.95 upper 0.95    Chisq            p method
# (Intercept)  0.9400042 0.5800077 -0.1605997  2.5479172  2.768364 0.0961443605      2
# Rv2385      -0.6915674 0.3151072 -2.3527547 -0.2336899 12.802717 0.0003461163      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=12.80272 on 1 df, p=0.0003461163, n=29
# Wald test = 4.93871 on 1 df, p = 0.02626164

Rv2385_firth_validate <- validate_firth(Outcome2 ~ Rv2385, data = model_df2, outcome = "Outcome2", B = 2000)
Rv2385_firth_validate
#                 Metric      Apparent      Optimism   Corrected
# 1                   AUC 9.242424e-01  0.0008383838 0.92340404
# 2                   Dxy 8.484848e-01  0.0016767677 0.84680808
# 3                 Brier 1.360103e-01 -0.0169040401 0.15291436
# 4 Intercept.(Intercept) 5.492091e-07 -0.0606887251 0.06068927
# 5              Slope.lp 9.999985e-01  0.1624187230 0.83757977
attr(Rv2385_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Rv2385_firth_df <- model_df2
Rv2385_firth_df$Predicted_prob <- predict(Rv2385_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Rv2385_firth_calibrate <- calibrate_firth(Outcome2 ~ Rv2385, data = model_df2, B = 2000)
Rv2385_firth_calibrate
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Rv2385_v1.pdf", width = 5, height = 4)
plot(Rv2385_firth_calibrate, main = "Rv2385")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Rv2385_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Rv2385_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Rv2385_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                  data=Rv2385_firth_df,
                                  fitted.risk=TRUE, 
                                  study.design="cohort", 
                                  bootstraps=2000)
plot_decision_curve(Rv2385_firth_DC,  curve.names = "Rv2385_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Rv2385_firth_ROC <- roc(response = Rv2385_firth_df$Outcome2, predictor = Rv2385_firth_df$Predicted_prob)
plot(Rv2385_firth_ROC, print.auc = TRUE, col = "blue")
auc(Rv2385_firth_ROC)
# Area under the curve: 0.9242
set.seed(42)
ci.auc(Rv2385_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 0.798-1 (2000 stratified bootstrap replicates)

###########################################################
######################### Rv1138c #########################
# Added 8/13/26

# Firth logistic regression
Rv1138c_firth <- logistf(Outcome2 ~ Rv1138c, data = model_df2)
Rv1138c_firth
summary(Rv1138c_firth)
#                   coef  se(coef) lower 0.95 upper 0.95    Chisq            p method
# (Intercept) -2.0469281 0.7247103 -3.7397902 -0.7802491 11.29721 0.0007762371      2
# Rv1138c      0.4348161 0.1475323  0.1718225  0.7727013 11.65520 0.0006402328      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=11.6552 on 1 df, p=0.0006402328, n=29
# Wald test = 9.426335 on 1 df, p = 0.002138914

Rv1138c_firth_validate <- validate_firth(Outcome2 ~ Rv1138c, data = model_df2, outcome = "Outcome2", B = 2000)
Rv1138c_firth_validate
#                 Metric      Apparent      Optimism   Corrected
# 1                   AUC  8.989899e-01  0.0003890834  0.89860082
# 2                   Dxy  7.979798e-01  0.0007781669  0.79720163
# 3                 Brier  1.482603e-01 -0.0146987920  0.16295913
# 4 Intercept.(Intercept) -4.228946e-08  0.0134037402 -0.01340378
# 5              Slope.lp  9.999999e-01  0.0201420037  0.97985792
attr(Rv1138c_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Rv1138c_firth_df <- model_df2
Rv1138c_firth_df$Predicted_prob <- predict(Rv1138c_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Rv1138c_firth_calibrate <- calibrate_firth(Outcome2 ~ Rv1138c, data = model_df2, B = 2000)
Rv1138c_firth_calibrate
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Rv1138c_v1.pdf", width = 5, height = 4)
plot(Rv1138c_firth_calibrate, main = "Rv1138c")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Rv1138c_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Rv1138c_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Rv1138c_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                  data=Rv1138c_firth_df,
                                  fitted.risk=TRUE, 
                                  study.design="cohort", 
                                  bootstraps=2000)
plot_decision_curve(Rv1138c_firth_DC,  curve.names = "Rv1138c_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Rv1138c_firth_ROC <- roc(response = Rv1138c_firth_df$Outcome2, predictor = Rv1138c_firth_df$Predicted_prob)
plot(Rv1138c_firth_ROC, print.auc = TRUE, col = "blue")
auc(Rv1138c_firth_ROC)
# Area under the curve: 0.9899
set.seed(42)
ci.auc(Rv1138c_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 0.7677-0.9848 (2000 stratified bootstrap replicates)

###########################################################
########################## Rv3747 #########################
# Updated 8/13/26

# Firth logistic regression
Rv3747_firth <- logistf(Outcome2 ~ Rv3747, data = model_df2)
Rv3747_firth
summary(Rv3747_firth)
#                   coef se(coef) lower 0.95 upper 0.95    Chisq            p method
# (Intercept)  3.790842  1.50847   1.292399   8.145966 11.05993 8.821359e-04      2
# Rv3747      -3.266807  1.25387  -6.797440  -1.135726 22.11901 2.562589e-06      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=22.11901 on 1 df, p=2.562589e-06, n=29
# Wald test = 6.944557 on 1 df, p = 0.008407465

Rv3747_firth_validate <- validate_firth(Outcome2 ~ Rv3747, data = model_df2, outcome = "Outcome2", B = 2000)
Rv3747_firth_validate
#                 Metric      Apparent      Optimism   Corrected
# 1                   AUC  9.797980e-01  0.0004040404  0.979393939
# 2                   Dxy  9.595960e-01  0.0008080808  0.958787879
# 3                 Brier  5.399058e-02 -0.0137716776  0.067762253
# 4 Intercept.(Intercept) -2.395257e-16  0.0089168572 -0.008916857
# 5              Slope.lp  1.000000e+00  0.0499612098  0.950038790
attr(Rv3747_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Rv3747_firth_df <- model_df2
Rv3747_firth_df$Predicted_prob <- predict(Rv3747_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Rv3747_firth_calibrate <- calibrate_firth(Outcome2 ~ Rv3747, data = model_df2, B = 2000)
Rv3747_firth_calibrate
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Rv3747_v1.pdf", width = 5, height = 4)
plot(Rv3747_firth_calibrate, main = "Rv3747")
# dev.off()


# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Rv3747_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Rv3747_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Rv3747_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                  data=Rv3747_firth_df,
                                  fitted.risk=TRUE, 
                                  study.design="cohort", 
                                  bootstraps=2000)
plot_decision_curve(Rv3747_firth_DC,  curve.names = "Rv3747_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Rv3747_firth_ROC <- roc(response = Rv3747_firth_df$Outcome2, predictor = Rv3747_firth_df$Predicted_prob)
plot(Rv3747_firth_ROC, print.auc = TRUE, col = "blue")
auc(Rv3747_firth_ROC)
# Area under the curve: 0.9798
set.seed(42)
ci.auc(Rv3747_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 0.9242-1 (2000 stratified bootstrap replicates)



###########################################################
######################### Rv2075c #########################
# Updated 8/13/26

# Firth logistic regression
Rv2075c_firth <- logistf(Outcome2 ~ Rv2075c, data = model_df2)
Rv2075c_firth
summary(Rv2075c_firth)
#                     coef  se(coef) lower 0.95 upper 0.95    Chisq            p method
# (Intercept)  1.7790945 0.7784898  0.3875992  3.5807510  6.536003 1.057126e-02      2
# Rv2075c     -0.5099922 0.1581548 -0.8999835 -0.2383791 16.816585 4.117189e-05      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=16.81659 on 1 df, p=4.117189e-05, n=29
# Wald test = 10.43495 on 1 df, p = 0.001236529

Rv2075c_firth_validate <- validate_firth(Outcome2 ~ Rv2075c, data = model_df2, outcome = "Outcome2", B = 2000)
Rv2075c_firth_validate
#                 Metric      Apparent      Optimism  Corrected
# 1                   AUC 9.393939e-01 -8.838384e-05 0.93948232
# 2                   Dxy 8.787879e-01 -1.767677e-04 0.87896465
# 3                 Brier 1.059672e-01 -1.273950e-02 0.11870671
# 4 Intercept.(Intercept) 3.807584e-08 -3.165348e-02 0.03165351
# 5              Slope.lp 9.999999e-01  6.172818e-02 0.93827175
attr(Rv2075c_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Rv2075c_firth_df <- model_df2
Rv2075c_firth_df$Predicted_prob <- predict(Rv2075c_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Rv2075c_firth_calibrate <- calibrate_firth(Outcome2 ~ Rv2075c, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Rv2075c_v1.pdf", width = 5, height = 4)
plot(Rv2075c_firth_calibrate, main = "Rv2075c")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Rv2075c_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Rv2075c_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Rv2075c_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                   data=Rv2075c_firth_df,
                                   fitted.risk=TRUE, 
                                   study.design="cohort", 
                                   bootstraps=2000)
plot_decision_curve(Rv2075c_firth_DC,  curve.names = "Rv2075c_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Rv2075c_firth_ROC <- roc(response = Rv2075c_firth_df$Outcome2, predictor = Rv2075c_firth_df$Predicted_prob)
plot(Rv2075c_firth_ROC, print.auc = TRUE, col = "blue")
auc(Rv2075c_firth_ROC)
# Area under the curve: 0.9394
set.seed(42)
ci.auc(Rv2075c_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 0.8332-1 (2000 stratified bootstrap replicates)

###########################################################
########################### Rv0161 ########################

# Firth logistic regression
Rv0161_firth <- logistf(Outcome2 ~ Rv0161, data = model_df2)
Rv0161_firth
summary(Rv0161_firth)
#                     coef  se(coef) lower 0.95 upper 0.95    Chisq           p method
# (Intercept)  1.9320601 0.8082331  0.4952078   4.022031  7.444338 6.363632e-03      2
# Rv0161      -0.7569889 0.2751607 -1.7064614  -0.324150 19.836320 8.436432e-06      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=19.83632 on 1 df, p=8.436432e-06, n=29
# Wald test = 7.899802 on 1 df, p = 0.004944021

Rv0161_firth_validate <- validate_firth(Outcome2 ~ Rv0161, data = model_df2, outcome = "Outcome2", B = 2000)
Rv0161_firth_validate
#                 Metric      Apparent      Optimism  Corrected
# 1                   AUC 9.595960e-01 -0.0004772727 0.96007323
# 2                   Dxy 9.191919e-01 -0.0009545455 0.92014646
# 3                 Brier 8.729436e-02 -0.0136964767 0.10099083
# 4 Intercept.(Intercept) 8.368190e-07 -0.0271150201 0.02711586
# 5              Slope.lp 9.999990e-01  0.0907939008 0.90920509
attr(Rv0161_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Rv0161_firth_df <- model_df2
Rv0161_firth_df$Predicted_prob <- predict(Rv0161_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Rv0161_firth_calibrate <- calibrate_firth(Outcome2 ~ Rv0161, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Rv0161_v1.pdf", width = 5, height = 4)
plot(Rv0161_firth_calibrate, main = "Rv0161")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Rv0161_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Rv0161_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Rv0161_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                  data=Rv0161_firth_df,
                                  fitted.risk=TRUE, 
                                  study.design="cohort", 
                                  bootstraps=2000)
plot_decision_curve(Rv0161_firth_DC,  curve.names = "Rv0161_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Rv0161_firth_ROC <- roc(response = Rv0161_firth_df$Outcome2, predictor = Rv0161_firth_df$Predicted_prob)
plot(Rv0161_firth_ROC, print.auc = TRUE, col = "blue")
auc(Rv0161_firth_ROC)
# Area under the curve: 0.9596
set.seed(42)
ci.auc(Rv0161_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 0.8687-1 (2000 stratified bootstrap replicates)

###########################################################
##################### Rv3747 + Rv2075c ####################
# Added 8/12/26

# Firth logistic regression
Rv3747.Rv2075c_firth <- logistf(Outcome2 ~ Rv3747 + Rv2075c, data = model_df2)
Rv3747.Rv2075c_firth
summary(Rv3747.Rv2075c_firth)
#                   coef  se(coef) lower 0.95  upper 0.95       Chisq            p method
# (Intercept)  5.3011602 2.4031790   1.672052 21.58955619 12.542998 0.0003976938      2
# Rv3747      -3.3251568 1.5283132 -12.520086 -0.61143309  7.977606 0.0047359530      2
# Rv2075c     -0.4002594 0.2589533  -2.025039  0.06472683  2.707275 0.0998915009      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=23.78765 on 2 df, p=6.832468e-06, n=29
# Wald test = 5.110335 on 2 df, p = 0.07767924

Rv3747.Rv2075c_firth_validate <- validate_firth(Outcome2 ~ Rv3747 + Rv2075c, data = model_df2, outcome = "Outcome2", B = 2000)
Rv3747.Rv2075c_firth_validate
#                   Metric    Apparent     Optimism  Corrected
# 1                   AUC  0.99494949  0.009431148  0.98551835
# 2                   Dxy  0.98989899  0.018862297  0.97103669
# 3                 Brier  0.03182053 -0.025191011  0.05701154
# 4 Intercept.(Intercept) -0.03631741 -0.004808007 -0.03150941
# 5              Slope.lp  1.12128068  0.300321877  0.82095880
attr(Rv3747.Rv2075c_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Rv3747.Rv2075c_firth_df <- model_df2
Rv3747.Rv2075c_firth_df$Predicted_prob <- predict(Rv3747.Rv2075c_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Rv3747.Rv2075c_firth_calibrate <- calibrate_firth(Outcome2 ~ Rv3747 + Rv2075c, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Rv3747.Rv2075c_v1.pdf", width = 5, height = 4)
plot(Rv3747.Rv2075c_firth_calibrate, main = "Rv3747 + Rv2075c")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Rv3747.Rv2075c_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Rv3747.Rv2075c_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Rv3747.Rv2075c_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                          data=Rv3747.Rv2075c_firth_df,
                                          fitted.risk=TRUE, 
                                          study.design="cohort", 
                                          bootstraps=2000)
plot_decision_curve(Rv3747.Rv2075c_firth_DC,  curve.names = "Rv3747.Rv2075c_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Rv3747.Rv2075c_firth_ROC <- roc(response = Rv3747.Rv2075c_firth_df$Outcome2, predictor = Rv3747.Rv2075c_firth_df$Predicted_prob)
plot(Rv3747.Rv2075c_firth_ROC, print.auc = TRUE, col = "blue")
auc(Rv3747.Rv2075c_firth_ROC)
# Area under the curve: 0.9949
set.seed(42)
ci.auc(Rv3747.Rv2075c_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 0.9697-1 (2000 stratified bootstrap replicates)

# Clinical Impact Plot - not sure what this means
plot_clinical_impact(Rv3747.Rv2075c_firth_DC, xlim = c(0, .20), col = c("black", "blue"))


# Check that these genes are not just a proxy for treatment arm
kruskal.test(Rv0161 ~ Arm2, data=model_df2)
kruskal.test(Rv2075c ~ Arm2, data=model_df2)
kruskal.test(Rv3747 ~ Arm2, data=model_df2)
kruskal.test(Rv0161 ~ Arm3, data=model_df2)
kruskal.test(Rv2075c ~ Arm3, data=model_df2)
kruskal.test(Rv3747 ~ Arm3, data=model_df2)
kruskal.test(Rv0161 ~ Outcome2, data=subset(model_df2, Arm2=="A"))
kruskal.test(Rv0161 ~ Outcome2, data=subset(model_df2, Arm2=="BorC"))
kruskal.test(Rv3747 ~ Outcome2, data=subset(model_df2, Arm2=="A"))
kruskal.test(Rv3747 ~ Outcome2, data=subset(model_df2, Arm2=="BorC"))


###########################################################
##################### Rv3747 + Rv2385 #####################
# Added 8/13/26

# Firth logistic regression
Rv3747.Rv2385_firth <- logistf(Outcome2 ~ Rv3747 + Rv2385, data = model_df2)
Rv3747.Rv2385_firth
summary(Rv3747.Rv2385_firth)
#                   coef  se(coef) lower 0.95  upper 0.95       Chisq            p method
# (Intercept)  5.2722351 2.3781393   1.681672 17.4181781 12.293710 0.0004544874      2
# Rv3747      -3.7934932 1.6607459 -11.455333 -1.0800028 11.092671 0.0008666955      2
# Rv2385      -0.4471503 0.2998847  -3.045967  0.1093481  2.464373 0.1164541883      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=23.2327 on 2 df, p=9.017436e-06, n=29
# Wald test = 5.692311 on 2 df, p = 0.05806714

Rv3747.Rv2385_firth_validate <- validate_firth(Outcome2 ~ Rv3747 + Rv2385, data = model_df2, outcome = "Outcome2", B = 2000)
Rv3747.Rv2385_firth_validate
#                   Metric    Apparent     Optimism  Corrected
# 1                   AUC  0.98989899  0.01785371  0.97204528
# 2                   Dxy  0.97979798  0.03570742  0.94409056
# 3                 Brier  0.03155968 -0.02595452  0.05751421
# 4 Intercept.(Intercept) -0.13075792  0.11229020 -0.24304812
# 5              Slope.lp  1.09224381  0.17508297  0.91716084
attr(Rv3747.Rv2385_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Rv3747.Rv2385_firth_df <- model_df2
Rv3747.Rv2385_firth_df$Predicted_prob <- predict(Rv3747.Rv2385_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Rv3747.Rv2385_firth_calibrate <- calibrate_firth(Outcome2 ~ Rv3747 + Rv2385, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Rv3747.Rv2385_v1.pdf", width = 5, height = 4)
plot(Rv3747.Rv2385_firth_calibrate, main = "Rv3747 + Rv2075c")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Rv3747.Rv2385_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Rv3747.Rv2385_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Rv3747.Rv2385_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                         data=Rv3747.Rv2385_firth_df,
                                         fitted.risk=TRUE, 
                                         study.design="cohort", 
                                         bootstraps=2000)
plot_decision_curve(Rv3747.Rv2385_firth_DC,  curve.names = "Rv3747.Rv2385_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Rv3747.Rv2385_firth_ROC <- roc(response = Rv3747.Rv2385_firth_df$Outcome2, predictor = Rv3747.Rv2385_firth_df$Predicted_prob)
plot(Rv3747.Rv2385_firth_ROC, print.auc = TRUE, col = "blue")
auc(Rv3747.Rv2385_firth_ROC)
# Area under the curve: 0.9899
set.seed(42)
ci.auc(Rv3747.Rv2385_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 0.9545-1 (2000 stratified bootstrap replicates)

###########################################################
################ Rv3747 + Rv0161 + Rv2075c ################
# Added 8/13/26

# Firth logistic regression
Rv3747.Rv0161.Rv2075c_firth <- logistf(Outcome2 ~ Rv3747 + Rv0161 + Rv2075c, data = model_df2)
Rv3747.Rv0161.Rv2075c_firth
summary(Rv3747.Rv0161.Rv2075c_firth)
#                   coef  se(coef) lower 0.95  upper 0.95       Chisq            p method
# (Intercept)  4.7578797 2.2226891   1.634979 34.43585833 14.7782341 0.0001209232      2
# Rv3747      -0.1151323 0.2485085  -4.328966  0.54098690  0.1882756 0.6643557336      2
# Rv0161      -0.5112694 0.2164944  -2.754451 -0.05923748  4.8245323 0.0280574767      2
# Rv2075c     -0.5250751 0.2863655  -4.094819  0.04912236  3.2629493 0.0708615692      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=27.20335 on 3 df, p=5.336827e-06, n=29
# Wald test = 9.510712 on 3 df, p = 0.02321767

Rv3747.Rv0161.Rv2075c_firth_validate <- validate_firth(Outcome2 ~ Rv3747 + Rv0161 + Rv2075c, data = model_df2, outcome = "Outcome2", B = 2000)
Rv3747.Rv0161.Rv2075c_firth_validate
#                   Metric    Apparent     Optimism  Corrected
# 1                   AUC  1.00000000  0.001118687 0.99888131
# 2                   Dxy  1.00000000  0.002237374 0.99776263
# 3                 Brier  0.01804828 -0.010870646 0.02891892
# 4 Intercept.(Intercept) -0.01282859 -0.128122670 0.11529408
# 5              Slope.lp  1.41275391 -0.086770256 1.49952416
attr(Rv3747.Rv0161.Rv2075c_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Rv3747.Rv0161.Rv2075c_firth_df <- model_df2
Rv3747.Rv0161.Rv2075c_firth_df$Predicted_prob <- predict(Rv3747.Rv0161.Rv2075c_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Rv3747.Rv0161.Rv2075c_firth_calibrate <- calibrate_firth(Outcome2 ~ Rv3747 + Rv0161 + Rv2075c, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Rv3747.Rv0161.Rv2075c_v1.pdf", width = 5, height = 4)
plot(Rv3747.Rv0161.Rv2075c_firth_calibrate, main = "Rv3747 + Rv0161 + Rv2075c")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Rv3747.Rv0161.Rv2075c_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Rv3747.Rv0161.Rv2075c_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Rv3747.Rv0161.Rv2075c_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                      data=Rv3747.Rv0161.Rv2075c_firth_df,
                                      fitted.risk=TRUE, 
                                      study.design="cohort", 
                                      bootstraps=2000)
plot_decision_curve(Rv3747.Rv0161.Rv2075c_firth_DC,  curve.names = "Rv3747.Rv0161.Rv2075c_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Rv3747.Rv0161.Rv2075c_firth_ROC <- roc(response = Rv3747.Rv0161.Rv2075c_firth_df$Outcome2, predictor = Rv3747.Rv0161.Rv2075c_firth_df$Predicted_prob)
plot(Rv3747.Rv0161.Rv2075c_firth_ROC, print.auc = TRUE, col = "blue")
auc(Rv3747.Rv0161.Rv2075c_firth_ROC)
# Area under the curve: 1
set.seed(42)
ci.auc(Rv3747.Rv0161.Rv2075c_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 1-1 (2000 stratified bootstrap replicates)

###########################################################
################ Rv3747 + Rv1138c + Rv2075c ################
# Updated 8/14/26

# Firth logistic regression
Rv3747.Rv1138c.Rv2075c_firth <- logistf(Outcome2 ~ Rv3747 + Rv1138c + Rv2075c, data = model_df2)
Rv3747.Rv1138c.Rv2075c_firth
summary(Rv3747.Rv1138c.Rv2075c_firth)
#                   coef  se(coef) lower 0.95  upper 0.95       Chisq            p method
# (Intercept)  0.1896737 1.2006416 -3.4844960 3.62966904 0.02035335 0.88655462      2
# Rv3747      -0.9408102 0.5576070 -4.0422900 0.05167249 3.19919114 0.07367470      2
# Rv1138c      1.0970457 0.5473741  0.1637378 3.98741942 4.75252072 0.02925541      2
# Rv2075c     -0.2101254 0.2619000 -2.2385677 0.58192692 0.44312206 0.50561941      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=26.23121 on 3 df, p=8.531265e-06, n=29
# Wald test = 7.825992 on 3 df, p = 0.0497482

Rv3747.Rv1138c.Rv2075c_firth_validate <- validate_firth(Outcome2 ~ Rv3747 + Rv1138c + Rv2075c, data = model_df2, outcome = "Outcome2", B = 2000)
Rv3747.Rv1138c.Rv2075c_firth_validate
#                   Metric    Apparent     Optimism  Corrected
# 1                   AUC  1.0000000  0.002666667  0.99733333
# 2                   Dxy  1.0000000  0.005333333  0.99466667
# 3                 Brier  0.0112750 -0.015047029  0.02632202
# 4 Intercept.(Intercept) -0.2502603  0.015694102 -0.26595440
# 5              Slope.lp  1.4254421 -0.247753154  1.67319522
attr(Rv3747.Rv1138c.Rv2075c_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Rv3747.Rv1138c.Rv2075c_firth_df <- model_df2
Rv3747.Rv1138c.Rv2075c_firth_df$Predicted_prob <- predict(Rv3747.Rv1138c.Rv2075c_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Rv3747.Rv1138c.Rv2075c_firth_calibrate <- calibrate_firth(Outcome2 ~ Rv3747 + Rv1138c + Rv2075c, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Rv3747.Rv1138c.Rv2075c_v1.pdf", width = 5, height = 4)
plot(Rv3747.Rv1138c.Rv2075c_firth_calibrate, main = "Rv3747 + Rv2385 + Rv2075c")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Rv3747.Rv1138c.Rv2075c_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Rv3747.Rv1138c.Rv2075c_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Rv3747.Rv1138c.Rv2075c_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                      data=Rv3747.Rv1138c.Rv2075c_firth_df,
                                      fitted.risk=TRUE, 
                                      study.design="cohort", 
                                      bootstraps=2000)
plot_decision_curve(Rv3747.Rv1138c.Rv2075c_firth_DC,  curve.names = "Rv3747.Rv1138c.Rv2075c_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Rv3747.Rv1138c.Rv2075c_firth_ROC <- roc(response = Rv3747.Rv1138c.Rv2075c_firth_df$Outcome2, predictor = Rv3747.Rv1138c.Rv2075c_firth_df$Predicted_prob)
plot(Rv3747.Rv1138c.Rv2075c_firth_ROC, print.auc = TRUE, col = "blue")
auc(Rv3747.Rv1138c.Rv2075c_firth_ROC)
# Area under the curve: 1
set.seed(42)
ci.auc(Rv3747.Rv1138c.Rv2075c_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 1-1 (2000 stratified bootstrap replicates)


###########################################################
############ Rv3747 + Rv0161+ Rv1138c + Rv2075c ###########
# Updated 8/14/26

# Firth logistic regression
Rv3747.Rv0161.Rv1138c.Rv2075c_firth <- logistf(Outcome2 ~ Rv3747 + Rv0161 + Rv1138c + Rv2075c, data = model_df2)
Rv3747.Rv0161.Rv1138c.Rv2075c_firth
summary(Rv3747.Rv0161.Rv1138c.Rv2075c_firth)
#                   coef  se(coef) lower 0.95  upper 0.95       Chisq            p method
# (Intercept)  2.69717448 1.9829403 -1.7671912 16.0549864 1.51913364 0.2177506      2
# Rv3747      -0.08489781 0.2472991 -2.3879497  0.9116303 0.07504919 0.7841222      2
# Rv0161      -0.35433798 0.2012110 -1.3611223  0.2283150 1.74025451 0.1871066      2
# Rv1138c      0.28462016 0.2378478 -0.3841507  2.3724294 1.00644375 0.3157563      2
# Rv2075c     -0.40991125 0.2550987 -2.7426601  0.2182260 2.00703802 0.1565707      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=26.62719 on 4 df, p=2.364443e-05, n=29
# Wald test = 11.58488 on 4 df, p = 0.02072056

Rv3747.Rv0161.Rv1138c.Rv2075c_firth_validate <- validate_firth(Outcome2 ~ Rv3747 + Rv0161 + Rv1138c + Rv2075c, data = model_df2, outcome = "Outcome2", B = 2000)
Rv3747.Rv0161.Rv1138c.Rv2075c_firth_validate
#                   Metric    Apparent     Optimism  Corrected
# 1                   AUC  1.00000000  0.0009065657 0.99909343
# 2                   Dxy  1.00000000  0.0018131313 0.99818687
# 3                 Brier  0.01347070 -0.0092931785 0.02276388
# 4 Intercept.(Intercept) -0.04859951 -0.1206967384 0.07209722
# 5              Slope.lp  1.47297601 -0.1617848487 1.63476086
attr(Rv3747.Rv0161.Rv1138c.Rv2075c_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Rv3747.Rv0161.Rv1138c.Rv2075c_firth_df <- model_df2
Rv3747.Rv0161.Rv1138c.Rv2075c_firth_df$Predicted_prob <- predict(Rv3747.Rv0161.Rv1138c.Rv2075c_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Rv3747.Rv0161.Rv1138c.Rv2075c_firth_calibrate <- calibrate_firth(Outcome2 ~ Rv3747 + Rv0161 + Rv1138c + Rv2075c, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Rv3747.Rv0161.Rv1138c.Rv2075c_v1.pdf", width = 5, height = 4)
plot(Rv3747.Rv0161.Rv1138c.Rv2075c_firth_calibrate, main = "Rv3747 + Rv0161 + Rv2385 + Rv2075c")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Rv3747.Rv0161.Rv1138c.Rv2075c_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Rv3747.Rv0161.Rv1138c.Rv2075c_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Rv3747.Rv0161.Rv1138c.Rv2075c_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                                  data=Rv3747.Rv0161.Rv1138c.Rv2075c_firth_df,
                                                  fitted.risk=TRUE, 
                                                  study.design="cohort", 
                                                  bootstraps=2000)
plot_decision_curve(Rv3747.Rv0161.Rv1138c.Rv2075c_firth_DC,  curve.names = "Rv3747.Rv0161.Rv1138c.Rv2075c_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Rv3747.Rv0161.Rv1138c.Rv2075c_firth_ROC <- roc(response = Rv3747.Rv0161.Rv1138c.Rv2075c_firth_df$Outcome2, predictor = Rv3747.Rv0161.Rv1138c.Rv2075c_firth_df$Predicted_prob)
plot(Rv3747.Rv0161.Rv1138c.Rv2075c_firth_ROC, print.auc = TRUE, col = "blue")
auc(Rv3747.Rv0161.Rv1138c.Rv2075c_firth_ROC)
# Area under the curve: 1
set.seed(42)
ci.auc(Rv3747.Rv0161.Rv1138c.Rv2075c_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 1-1 (2000 stratified bootstrap replicates)

###########################################################
################# BASE + ARM + THREE GENES ################
# Updated 8/14/26

# Firth logistic regression
Base.Arm.Rv3747.Rv0161.Rv2075c_firth <- logistf(Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age + Arm2 + Arm3 + Rv3747 + Rv0161 + Rv2075c, data = model_df2)
Base.Arm.Rv3747.Rv0161.Rv2075c_firth
summary(Base.Arm.Rv3747.Rv0161.Rv2075c_firth)
#                     coef   se(coef)  lower 0.95 upper 0.95       Chisq         p method
# (Intercept)  4.004467089 4.45112253 -12.14521759 19.71098956 0.5661604537 0.45178870      2
# TTD          0.065340818 0.11947360  -0.36670584  0.53477205 0.2198978377 0.63911777      2
# XpertCT_wk0 -0.003553512 0.18367118  -0.91727785  0.55208565 0.0002558914 0.98723710      2
# BMI         -0.168780599 0.20388921  -1.05031935  0.92885760 0.4567502760 0.49914667      2
# Age          0.039562920 0.04204398  -0.04224964  0.24210611 0.8169677953 0.36606821      2
# Arm2BorC    -0.770761477 1.64256053  -7.26870959  4.78998775 0.1750726786 0.67564235      2
# Arm3C        1.578744535 1.81362663  -3.82657313  8.42336922 0.5670240887 0.45144389      2
# Rv3747       0.085105592 0.25640588  -1.27508855  0.83780546 0.0736167472 0.78614180      2
# Rv0161      -0.386503067 0.16989865  -1.13436917 -0.07432395 6.3570066834 0.01169186      2
# Rv2075c     -0.356209207 0.19682019  -1.23856208  0.02375656 3.3409527135 0.06757547      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=23.03702 on 9 df, p=0.006113719, n=29
# Wald test = 13.52617 on 9 df, p = 0.1402069

# Odds Ratios
exp(cbind(OR = coef(Base.Arm.Rv3747.Rv0161.Rv2075c_firth), confint(Base.Arm.Rv3747.Rv0161.Rv2075c_firth)))
#                   OR    Lower 95%    Upper 95%
# (Intercept) 54.8425904 5.313724e-06 3.633909e+08
# TTD          1.0675228 6.930135e-01 1.707059e+00
# XpertCT_wk0  0.9964528 3.996053e-01 1.736872e+00
# BMI          0.8446942 3.498260e-01 2.531615e+00
# Age          1.0403560 9.586304e-01 1.273929e+00
# Arm2BorC     0.4626606 6.970108e-04 1.202999e+02
# Arm3C        4.8488644 2.178414e-02 4.552215e+03
# Rv3747       1.0888320 2.794062e-01 2.311289e+00
# Rv0161       0.6794286 3.216249e-01 9.283709e-01
# Rv2075c      0.7003261 2.898006e-01 1.024041e+00

Base.Arm.Rv3747.Rv0161.Rv2075c_firth_validate <- validate_firth(Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age + Arm2 + Arm3 + Rv3747 + Rv0161 + Rv2075c, data = model_df2, outcome = "Outcome2", B = 2000)
Base.Arm.Rv3747.Rv0161.Rv2075c_firth_validate
#                   Metric    Apparent    Optimism  Corrected
# 1                   AUC  1.0000000  0.01621815  0.9837818
# 2                   Dxy  1.0000000  0.03243631  0.9675637
# 3                 Brier  0.0202700 -0.02849431  0.0487643
# 4 Intercept.(Intercept) -0.4361693  0.26030743 -0.6964768
# 5              Slope.lp  1.8243063  0.01721456  1.8070918
attr(Base.Arm.Rv3747.Rv0161.Rv2075c_firth_validate, "successful_bootstraps")
# 1733

# Make predictions based on the model
Base.Arm.Rv3747.Rv0161.Rv2075c_firth_df <- model_df2
Base.Arm.Rv3747.Rv0161.Rv2075c_firth_df$Predicted_prob <- predict(Base.Arm.Rv3747.Rv0161.Rv2075c_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Base.Arm.Rv3747.Rv0161.Rv2075c_firth_calibrate <- calibrate_firth(Outcome2 ~ TTD + XpertCT_wk0 + BMI + Age + Arm2 + Arm3 + Rv3747 + Rv0161 + Rv2075c, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Base.Arm.Rv3747.Rv0161.Rv2075c_v1.pdf", width = 5, height = 4)
plot(Base.Arm.Rv3747.Rv0161.Rv2075c_firth_calibrate, main = "Base+Arm+Rv3747.Rv0161.Rv2075c")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Base.Arm.Rv3747.Rv0161.Rv2075c_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Base.Arm.Rv3747.Rv0161.Rv2075c_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Base.Arm.Rv3747.Rv0161.Rv2075c_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                               data=Base.Arm.Rv3747.Rv0161.Rv2075c_firth_df,
                                               fitted.risk=TRUE, 
                                               study.design="cohort", 
                                               bootstraps=2000)
plot_decision_curve(Base.Arm.Rv3747.Rv0161.Rv2075c_firth_DC,  curve.names = "Base.Arm.Rv3747.Rv0161.Rv2075c_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Base.Arm.Rv3747.Rv0161.Rv2075c_firth_ROC <- roc(response = Base.Arm.Rv3747.Rv0161.Rv2075c_firth_df$Outcome2, predictor = Base.Arm.Rv3747.Rv0161.Rv2075c_firth_df$Predicted_prob)
plot(Base.Arm.Rv3747.Rv0161.Rv2075c_firth_ROC, print.auc = TRUE, col = "blue")
auc(Base.Arm.Rv3747.Rv0161.Rv2075c_firth_ROC)
# Area under the curve: 1
set.seed(42)
ci.auc(Base.Arm.Rv3747.Rv0161.Rv2075c_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 1-1 (2000 stratified bootstrap replicates)
# Warning message:
#   In ci.auc.roc(Base.Arm.Rv3747.Rv0161.Rv2075c_firth_ROC, method = "bootstrap",  :
#                   ci.auc() of a ROC curve with AUC == 1 is always 1-1 and can be misleading.

# Summary Graphs
boxplot(Rv0161 ~ Outcome2, data=model_df2)
boxplot(Rv2075c ~ Outcome2, data=model_df2)
boxplot(Rv3747 ~ Outcome2, data=model_df2)
table(model_df2$Arm2, model_df2$Outcome2)
table(model_df2$Arm3, model_df2$Outcome2)
boxplot(Rv3747 ~ Outcome2 + Arm3, data=model_df2)

###########################################################
#################### ARM + THREE GENES ####################
# Updated 8/14/26

# Firth logistic regression
Arm.Rv3747.Rv0161.Rv2075c_firth <- logistf(Outcome2 ~ Arm2 + Arm3 + Rv3747 + Rv0161 + Rv2075c, data = model_df2)
Arm.Rv3747.Rv0161.Rv2075c_firth
summary(Arm.Rv3747.Rv0161.Rv2075c_firth)
# logistf(formula = Outcome2 ~ Arm2 + Arm3 + Rv3747 + Rv0161 + 
#           Rv2075c, data = model_df2)
# 
# Model fitted by Penalized ML
# Coefficients:
#   coef  se(coef)   lower 0.95 upper 0.95       Chisq          p method
# (Intercept)  2.5510774 1.4660156  -0.01350323 18.5120391 3.794715476 0.05141461      2
# Arm2BorC    -0.2041663 2.0731477 -17.06499717  5.5965252 0.001201298 0.97235106      2
# Arm3C        1.6439965 2.1638194  -7.33249149 14.5185072 0.145564427 0.70281071      2
# Rv3747      -0.2240251 0.2240012  -2.11805010  0.2240693 0.918684364 0.33782065      2
# Rv0161      -0.4322717 0.1732527  -1.57134727 -0.1107233 6.273214490 0.01225766      2
# Rv2075c     -0.2565253 0.1881667  -2.18712370  0.1236040 1.760857124 0.18451764      2
# 
# Method: 1-Wald, 2-Profile penalized log-likelihood, 3-None
# 
# Likelihood ratio test=25.28633 on 5 df, p=0.0001226689, n=29
# Wald test = 11.21297 on 5 df, p = 0.04731707

Arm.Rv3747.Rv0161.Rv2075c_firth_validate <- validate_firth(Outcome2 ~ Arm2 + Arm3 + Rv3747 + Rv0161 + Rv2075c, data = model_df2, outcome = "Outcome2", B = 2000)
Arm.Rv3747.Rv0161.Rv2075c_firth_validate
#                   Metric    Apparent    Optimism  Corrected
# 1                   AUC  1.00000000  0.004709530  0.99529047
# 2                   Dxy  1.00000000  0.009419061  0.99058094
# 3                 Brier  0.02086237 -0.017416000  0.03827837
# 4 Intercept.(Intercept) -0.18582170  0.179028972 -0.36485067
# 5              Slope.lp  1.64374572  0.032549439  1.61119628
attr(Arm.Rv3747.Rv0161.Rv2075c_firth_validate, "successful_bootstraps")
# 1733

# Make predictions based on the model
Arm.Rv3747.Rv0161.Rv2075c_firth_df <- model_df2
Arm.Rv3747.Rv0161.Rv2075c_firth_df$Predicted_prob <- predict(Arm.Rv3747.Rv0161.Rv2075c_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
Arm.Rv3747.Rv0161.Rv2075c_firth_calibrate <- calibrate_firth(Outcome2 ~ Arm2 + Arm3 + Rv3747 + Rv0161 + Rv2075c, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/Arm.Rv3747.Rv0161.Rv2075c_v1.pdf", width = 5, height = 4)
plot(Arm.Rv3747.Rv0161.Rv2075c_firth_calibrate, main = "Arm+Rv3747.Rv0161.Rv2075c")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Arm.Rv3747.Rv0161.Rv2075c_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Arm.Rv3747.Rv0161.Rv2075c_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Arm.Rv3747.Rv0161.Rv2075c_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                               data=Arm.Rv3747.Rv0161.Rv2075c_firth_df,
                                               fitted.risk=TRUE, 
                                               study.design="cohort", 
                                               bootstraps=2000)
plot_decision_curve(Arm.Rv3747.Rv0161.Rv2075c_firth_DC,  curve.names = "Arm.Rv3747.Rv0161.Rv2075c_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Arm.Rv3747.Rv0161.Rv2075c_firth_ROC <- roc(response = Arm.Rv3747.Rv0161.Rv2075c_firth_df$Outcome2, predictor = Arm.Rv3747.Rv0161.Rv2075c_firth_df$Predicted_prob)
plot(Arm.Rv3747.Rv0161.Rv2075c_firth_ROC, print.auc = TRUE, col = "blue")
auc(Arm.Rv3747.Rv0161.Rv2075c_firth_ROC)
# Area under the curve: 1
set.seed(42)
ci.auc(Arm.Rv3747.Rv0161.Rv2075c_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 1-1 (2000 stratified bootstrap replicates)
# Warning message:
#   In ci.auc.roc(Arm.Rv3747.Rv0161.Rv2075c_firth_ROC, method = "bootstrap",  :
#                   ci.auc() of a ROC curve with AUC == 1 is always 1-1 and can be misleading.

# Summary Graphs
boxplot(Rv0161 ~ Outcome2, data=model_df2)
boxplot(Rv2075c ~ Outcome2, data=model_df2)
boxplot(Rv3747 ~ Outcome2, data=model_df2)
table(model_df2$Arm2, model_df2$Outcome2)
table(model_df2$Arm3, model_df2$Outcome2)
boxplot(Rv3747 ~ Outcome2 + Arm3, data=model_df2)


###########################################################
##################### TTD+Xpert FIRTH #####################

TTD.Xpert_firth <- logistf(Outcome2 ~ TTD + XpertCT_wk0, data = model_df2)
TTD.Xpert_firth
summary(TTD.Xpert_firth)
# Likelihood ratio test=0.967831 on 2 df, p=0.6163653, n=29

# Use Bootstrapping to validate
TTD.Xpert_firth_validate <- validate_firth(Outcome2 ~ TTD + XpertCT_wk0, data = model_df2, outcome = "Outcome2", B = 2000)
TTD.Xpert_firth_validate
#                   Metric    Apparent    Optimism  Corrected
# 1                   AUC 0.593434343  0.06601334  0.52742100
# 2                   Dxy 0.186868687  0.13202668  0.05484201
# 3                 Brier 0.226970566 -0.03188238  0.25885295
# 4 Intercept.(Intercept) 0.004821339  0.26550889 -0.26068755
# 5              Slope.lp 1.029388577  0.61491901  0.41446957
attr(TTD.Xpert_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
TTD.Xpert_firth_df <- model_df2
TTD.Xpert_firth_df$Predicted_prob <- predict(TTD.Xpert_firth, type="response")

# Determine model calibration
## Maybe don't have to do this, can use the intercept and slope from the validate function??
# source("Function_Firth_Validate.R")
TTD.Xpert_firth_calibrate <- calibrate_firth(Outcome2 ~ TTD + XpertCT_wk0, data = model_df2, B = 2000)
# pdf("Figures/PredictiveModeling/W2/CalibrationCurves/TTD.Xpert_v1.pdf", width = 5, height = 4)
plot(TTD.Xpert_firth_calibrate, main = "TTD + Xpert_W0")
# dev.off()

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(TTD.Xpert_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
TTD.Xpert_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
TTD.Xpert_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                     data=TTD.Xpert_firth_df,
                                     fitted.risk=TRUE, 
                                     study.design="cohort", 
                                     bootstraps=2000)
plot_decision_curve(TTD.Xpert_firth_DC,  curve.names = "TTD.Xpert_firth", lty = 1)

# ROC - Can the model rank patients correctly?
TTD.Xpert_firth_ROC <- roc(response = TTD.Xpert_firth_df$Outcome2, predictor = TTD.Xpert_firth_df$Predicted_prob)
plot(TTD.Xpert_firth_ROC, print.auc = TRUE, col = "blue")
auc(TTD.Xpert_firth_ROC)
# Area under the curve: 0.5934
set.seed(42)
ci.auc(TTD.Xpert_firth_ROC, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
# 95% CI: 0.3788-0.798 (2000 stratified bootstrap replicates)




###########################################################
##################### ARM ONLY FIRTH ######################

Arms_firth <- logistf(Outcome2 ~ Arm2 + Arm3, data = model_df2)
Arms_firth
summary(basic2_firth)
# Likelihood ratio test=2.035934 on 3 df, p=0.5649821, n=29

# Use Bootstrapping to validate
Arms_firth_validate <- validate_firth(Outcome2 ~Arm2 + Arm3, data = model_df2, outcome = "Outcome2", B = 2000)
Arms_firth_validate
# Metric    Apparent    Optimism  Corrected
# 1                   AUC 0.671717172  0.03731941  0.6343978
# 2                   Dxy 0.343434343  0.07463881  0.2687955
# 3                 Brier 0.210106537 -0.02675251  0.2368590
# 4 Intercept.(Intercept) 0.000879326  0.11115962 -0.1102803
# 5              Slope.lp 1.000925716  0.28424142  0.7166843
attr(Arms_firth_validate, "successful_bootstraps")
# 1733

# Make predictions based on the model
Arms_firth_df <- model_df2
Arms_firth_df$Predicted_prob <- predict(Arms_firth, type="response")

# Determine model calibration
# Can't do a LOESS smoother on this because it's just categorical variables

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Arms_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Arms_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilities clinically useful?
set.seed(42)
Arms_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                data=Arms_firth_df,
                                fitted.risk=TRUE, 
                                study.design="cohort", 
                                bootstraps=2000)
plot_decision_curve(Arms_firth_DC,  curve.names = "Arms_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Arms_firth_ROC <- roc(response = Arms_firth_df$Outcome2, predictor = Arms_firth_df$Predicted_prob)
plot(Arms_firth_ROC, print.auc = TRUE, col = "blue")
auc(Arms_firth_ROC)
# Area under the curve: 0.6717

###########################################################
################### LINEAGE ONLY FIRTH ####################

# Firth logistic regression
Lineage_firth <- logistf(Outcome2 ~ Lineage2, data = model_df2)
Lineage_firth
summary(Lineage_firth)
# (Intercept) Lineage2lineage4 
# -0.1251631       -0.7221347 
# 
# Likelihood ratio test=0.9366678 on 1 df, p=0.3331363, n=29

Lineage_firth_validate <- validate_firth(Outcome2 ~ Lineage2, data = model_df2, outcome = "Outcome2", B = 2000)
Lineage_firth_validate
# Metric     Apparent    Optimism  Corrected
# 1                   AUC 5.959596e-01  0.02156481  0.5743948
# 2                   Dxy 1.919192e-01  0.04312963  0.1487896
# 3                 Brier 2.273586e-01 -0.01731987  0.2446784
# 4 Intercept.(Intercept) 4.519022e-16  0.14392667 -0.1439267
# 5              Slope.lp 1.000000e+00  0.31593269  0.6840673
attr(Lineage_firth_validate, "successful_bootstraps")
# 2000

# Make predictions based on the model
Lineage_firth_df <- model_df2
Lineage_firth_df$Predicted_prob <- predict(Lineage_firth, type="response")

# Determine model calibration
# This doesn't work with just a binary predictor

# Re-calibrate the the prediction risks back to a 5% population prevalence.
# Because case cases were intentionally enriched, predicted probabilities were recalibrated to a 5% population prevalence before decision curve analysis.
sample_prev <- mean(model_df2$Outcome2==1)
pop_prev <- 0.05
p <- predict(Lineage_firth, type="response")
p_adj <- (p*(pop_prev/sample_prev)) / (p*(pop_prev/sample_prev) +
                                         (1-p)*((1-pop_prev)/(1-sample_prev))) # Bayes theorem Posterior odds=Likelihood ratio × Population odds
Lineage_firth_df$Adjusted_prob <- p_adj

# Decision curve using the re-calibrated probabilities - Are the probabilites clinically useful?
set.seed(42)
Lineage_firth_DC <- decision_curve(Outcome2 ~ Adjusted_prob, 
                                   data=Lineage_firth_df,
                                   fitted.risk=TRUE, 
                                   study.design="cohort", 
                                   bootstraps=2000)
plot_decision_curve(Lineage_firth_DC,  curve.names = "Lineage_firth", lty = 1)

# ROC - Can the model rank patients correctly?
Lineage_firth_ROC <- roc(response = Lineage_firth_df$Outcome2, predictor = Lineage_firth_df$Predicted_prob)
plot(Lineage_firth_ROC, print.auc = TRUE, col = "blue")
auc(Lineage_firth_ROC)
# Area under the curve: 0.596

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




