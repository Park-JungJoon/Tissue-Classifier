# Library 
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
library(ROCR)
library(gridExtra)
library(modelr)
library(tictoc)
library(doParallel)
library(foreach)
# Funtions
binomial_deviance <- function(y_obs, yhat){
  epsilon = 0.0001
  yhat = ifelse(yhat < epsilon, epsilon, yhat)
  yhat = ifelse(yhat > 1-epsilon, 1-epsilon, yhat)
  a = ifelse(y_obs==0, 0, y_obs * log(y_obs/yhat))
  b = ifelse(y_obs==1, 0, (1-y_obs) * log((1-y_obs)/(1-yhat)))
  return(2*sum(a + b))
}

predict.cv.glmnet=function(object,newx,s=c("lambda.1se","lambda.min"),...){
  if(is.numeric(s))lambda=s
  else
    if(is.character(s)){
      s=match.arg(s)
      lambda=object[[s]]
      names(lambda)=s
    }
  else stop("Invalid form for s")
  predict(object$glmnet.fit,newx,s=lambda,...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
} 


# Loading data, sample genes
gene_matrix <- fread("/eevee/val/jjpark/PAPER_RNA_SEQ_ATLAS/d_ts/TIS_T_DROP_INTER_GTEX_PTN_GCT.tsv", sep2 = "\t", header = TRUE, data.table = FALSE)
allgene <- names(gene_matrix)
using_gene <- sample(allgene, 5000)
using_gene <- c("V1", using_gene)

df <- subset(gene_matrix, select = using_gene)
df$V1 <- ifelse(df$V1 == "Skin", 1, 0)

#index
set.seed(0621)
n <- nrow(df)
idx <- 1:n
training_idx <- sample(idx, n * .60)
idx <- setdiff(idx, training_idx)
validate_idx <- sample(idx, n * .20)
test_idx <- setdiff(idx, validate_idx)
training <- df[training_idx,]
validation <- df[validate_idx,]
test <- df[test_idx,]

# Logicstic regression
tic("logistic regression fitting")
lm_skin <- glm(V1 ~ ., data = training, family = binomial, na.action = na.omit)
toc()
save.image(file = "/home/jjpark/lm_lasso_tree.RData")
# Logistic regression model evaluation
y_obs <- validation$V1
yhat_lm <- predict(lm_liver, newdata = validation, type='response')
pred_lm <- prediction(as.numeric(yhat_lm), as.numeric(y_obs))
performance(pred_lm, "auc")@y.values[[1]]
binomial_deviance(y_obs, yhat_lm) 
plot(performance(pred_lm, measure="tpr", x.measure='fpr'))
save.image(file = "/home/jjpark/lm_lasso_tree.RData")

# LASSO Model

xx <- model.matrix(V1 ~ . -1, df)
x <- xx[training_idx,]
y <- as.numeric(as.character(training$V1))
glimpse(x)
df_cvfit <- cv.glmnet(x,y, family = 'binomial')

#save.image(file = "/home/jjpark/lm_lasso_tree.RData")

save.image(file = "/home/jjpark/lm_lasso_tree.RData")

coef(df_cvfit, s = c("lambda.1se"))
coef(df_cvfit, s = c("lambda.min"))

# LASSO evaluation
predict.cv.glmnet(df_cvfit, s = "lambda.min", newx = x[1:10,], type = 'response')
yhat_glmnet <- predict(df_cvfit, s="lambda.min", newx = xx[validate_idx,], type = 'response')
yhat_glmnet<- yhat_glmnet[,1]
pred_glmnet <- prediction(yhat_glmnet, y_obs)
performance(pred_glmnet, "auc")@y.values[[1]]
binomial_deviance(y_obs, yhat_glmnet)



# Tree Model
tic("tree")
df_tr <- rpart(V1 ~ ., data = training)
toc()
printcp(df_tr)
summary(df_tr)
opar <- par(mfrow = c(1,1), xpd = NA)
plot(df_tr)
text(df_tr, use.n = T)
par(opar)

# tree evalutation
yhat_tr <- predict(df_tr, validation)
#yhat_tr <- yhat_tr[,"1"]
pred_tr <- prediction(yhat_tr, y_obs)
performance(pred_tr, "auc")@y.values[[1]]
binomial_deviance(y_obs, yhat_tr)


# Random Forest
set.seed(0625)
tic("Random Forest")
df_rf <- randomForest(V1 ~., training)
toc()
save.image(file = "/home/jjpark/lm_lasso_tree.RData")
rftoc <- toc()
opar <- par(mfrow=c(1,2))
plot(df_rf)
varImpPlot(df_rf)
par(opar)

save.image(file = "/home/jjpark/all_model_done.RData")
###
library(randomForest)
install.packages('doParallel')
library(doParallel)
install.packages('foreach')
library(foreach)

# Set the number of threads you want to use
num_threads <- 28

# Register the parallel backend with the desired number of threads
cl <- makeCluster(num_threads)
registerDoParallel(cl)

tic('multithread rf')
# Fit the random forest model in parallel
df_rf <- foreach(ntree = rep(4, num_threads), .combine = combine, .packages = 'randomForest') %dopar% {
  randomForest(V1 ~ ., data = training, ntree = ntree)
}
toc()
save.image(file = "/home/jjpark/lm_lasso_tree_rf.RData")
# Stop the parallel backend
stopCluster(cl)
registerDoSEQ()


load(file = "/home/jjpark/lm_lasso_tree_rf.RData")
# Random Forest Evaluation
yhat_rf <- predict(df_rf, newdata=validation)   #, type='prob')[,'1']
pred_rf <- prediction(yhat_rf, y_obs)
performance(pred_rf, "auc")@y.values[[1]]
binomial_deviance(y_obs, yhat_rf)

# Boosting
set.seed(0625)
df_toy_gbm <-
  training %>%
  mutate(V1 = as.numeric(as.character(V1)))

df_gbm <- gbm(V1 ~ ., data = df_toy_gbm, distribution = 'bernoulli', n.trees = 348, n.minobsinnode =1, cv.folds = 3, verbose = T)

(best_iter = gbm.perf(df_gbm, method = "cv"))
# best_iter = 348
save.image(file = "/home/jjpark/lm_lasso_tree_rf.RData")
# Boosting Evaluation
yhat_gbm <- predict(df_gbm, n.trees = best_iter, newdata = validation, type = 'response')
pred_gbm <- prediction(yhat_gbm, y_obs)
performance(pred_gbm, "auc")@y.values[[1]]
binomial_deviance(y_obs, yhat_gbm)


# Final model test

data.frame(method=c('lm', 'glmnet', 'rf', 'gbm'),
           auc = c(performance(pred_lm, "auc")@y.values[[1]],
                   performance(pred_glmnet, "auc")@y.values[[1]],
                   performance(pred_rf, "auc")@y.values[[1]],
                   performance(pred_gbm, "auc")@y.values[[1]]),
           bin_dev = c(binomial_deviance(y_obs, yhat_lm),
                       binomial_deviance(y_obs, yhat_glmnet),
                       binomial_deviance(y_obs, yhat_rf),
                       binomial_deviance(y_obs, yhat_gbm)))


perf_lm <- performance(pred_lm, measure = "tpr", x.measure = "fpr")
perf_glmnet <- performance(pred_glmnet, measure = "tpr", x.measure = "fpr")
perf_rf <- performance(pred_tr, measure = "tpr", x.measure = "fpr")
perf_gbm <- performance(pred_gbm, measure = "tpr", x.measure = "fpr")

plot(perf_lm, col='black', main="ROC Curve")
plot(perf_glmnet, add=T, col='blue')
plot(perf_rf, add=T, col='red')
plot(perf_gbm, add=T, col='cyan')
abline(0,1)
legend('bottomright', inset=.1,
       legend=c("GLM", "glmnet", "RF", "GBM"),
       col=c('black', 'blue', 'red', 'cyan'), lty=1, lwd=2)

save.image(file = "/home/jjpark/lm_lasso_tree_rf_gbm.RData")

load(file = "/home/jjpark/all_model_done.RData")
