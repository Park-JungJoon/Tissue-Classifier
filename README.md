# Tissue-Classifier

# Report 06-07
+ Dataset : GTEx protein coding gene expression matrix, sample name replaced in tissue name of the sample.

  <img width="700" alt="image" src="https://github.com/Park-JungJoon/Tissue-Classifier/assets/97942772/80b81ced-1caa-4859-9fe7-2e10bbd97fa8">

  + 17,211 rows (samples), 19,152 columns (genes), V1 column is tissue's name.

## Problem 1 : Memory consumption (Stackoverflow error)

```
 xx_gtex <- model.matrix(V1 ~. -1, gtex_trimmed)
```

  - First stackoverflow error; Making model matrix. 
  - Founded out memory consumption append to columns not rows.
  - Gene pre-filtering required.

### Solve : Gene pre-filtering by tau values.

```
gtex_taus <- read.table('/panpyro/alfa/jjpark/0526_r_glm_model/GTEX_TAU.tsv', sep = '\t', header = FALSE)
sorted_gtex_taus <- gtex_taus[order(gtex_taus$V2, decreasing = TRUE), ]
gtex_index_80 <- ceiling(0.8 * nrow(sorted_gtex_taus))
gtex_genes <- sorted_gtex_taus[1:gtex_index_80,]$V1
gtex_genes <- unlist(gtex_genes) 
gtex_genes <- c(gtex_genes, 'V1')
```
  + Using 80% of total genes. 
  + Tau is a value in 0-1 scale, 1 means tissue-specific and 0 means non-tissue-specific. 

## Progress 1 : Built LASSO model

![image](https://github.com/Park-JungJoon/Tissue-Classifier/assets/97942772/106810ad-9d01-4ce2-92b5-13b016df0fb9)

+ Built LASSO model with filtered dataset.

## Problem 2 : Making accuracy funtion work with multinomial classifiers.

```
predicted_classes <- colnames(yhat_valid)[apply(yhat_valid, 1, which.max)]
accuracy <- sum(predicted_classes == y_valid) / length(y_valid)
```
+ Using functions to calculate accuracy.

## Problem 3 : Another filtering required for linear regression model.

+ Multinomial Linear Regression model is available with "nnet" library in R packages.

```
linear_regression_gtex <- multinom( y ~., data = as.data.frame(x))
Error in nnet.default(X,Y,w,mask = mask, size = 0, skip = TRUE, softmax = TRUE)
  too many (41073) weights
```
