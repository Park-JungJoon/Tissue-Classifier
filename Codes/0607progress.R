library(dplyr)
library(ggplot2)
library(ISLR)
library(MASS)
library(glmnet)
library(randomForest)
library(gbm)
library(rpart)
library(boot)
library(data.table)
library(nnet)
library(splitstackshape)

#
binomial_deviance <- function(y_obs, yhat){
epsilon = 0.0001
yhat = ifelse(yhat < epsilon, epsilon, yhat)
yhat = ifelse(yhat > 1-epsilon, 1-epsilon, yhat)
a = ifelse(y_obs==0, 0, y_obs * log(y_obs/yhat))
b = ifelse(y_obs==1, 0, (1-y_obs) * log((1-y_obs)/(1-yhat)))
return(2*sum(a + b))
}

#gene cutting
gtex_taus <- read.table('/panpyro/alfa/jjpark/0526_r_glm_model/GTEX_TAU.tsv', sep = '\t', header = FALSE)
sorted_gtex_taus <- gtex_taus[order(gtex_taus$V2, decreasing = TRUE), ]
gtex_index_80 <- ceiling(0.8 * nrow(sorted_gtex_taus))
gtex_genes <- sorted_gtex_taus[1:gtex_index_80,]$V1
gtex_genes <- unlist(gtex_genes) 
gtex_genes <- c(gtex_genes, 'V1')

# gtex file read
gtex <- fread("/panpyro/alfa/jjpark/0526_r_glm_model/TIS_T_DROP_INTER_GTEX_PTN_GCT.tsv", sep2 = "\t", header = TRUE, data.table = FALSE)
gtex_trimmed <- subset(gtex, select = gtex_genes)

# gtex index preparing
set.seed(0601)
n <- nrow(gtex_trimmed)
idx <- 1:n
training_idx <- sample(idx, n * .60)
idx <- setdiff(idx, training_idx)
validate_idx <- sample(idx, length(idx) * 0.2)
test_idx <- setdiff(idx, validate_idx)

training_gtex<- gtex_trimmed[training_idx,]
validation_gtex<- gtex_trimmed[validate_idx,]
test_gtex<- gtex_trimmed[test_idx,]
# model matrix 
xx_gtex <- model.matrix(V1 ~. -1, gtex_trimmed)

#Linear Regression
x = xx_gtex[training_idx ,]
y = training_gtex$V1

linear_regression <- multinom(y ~ ., data = as.data.frame(x))


# lasso model
x = xx_gtex[training_idx ,]
y = training_gtex$V1
gtex_glmnet <- cv.glmnet(x,y, family = 'multinomial', alpha = 1)

# lasso model predict
x_valid <- xx_gtex[validate_idx, ]
y_valid <- validation_gtex$V1
yhat_valid <- predict(gtex_glmnet, newx = x_valid, type = 'response')

# Evaluate model performance
predicted_classes <- colnames(yhat_valid)[apply(yhat_valid, 1, which.max)]
accuracy <- sum(predicted_classes == y_valid) / length(y_valid)

#
pred_glmnet <- prediction(gtex_glmnet, y_valid)
perf_glmnet <- performance(pred_glment, measure ="tpr", x.measure = "fpr")
performance()
