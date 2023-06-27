# Tissue-Classifier

## 1. Data Collecting
+ Data downloaded from [GTEx](https://gtexportal.org/home/)
+ 17,382 Samples from 948 donors.
+ Samples were collected from 54 non-diseased tissue sites, from dead body.
+ Download gene count table, Used program list below.
  - Gencode version: v26
  - RNA-seq alignment: STAR v2.5.3a
  - Transcript quantification: RSEM v1.3.0
+ Additional meta data : age, sex.

</br>

## 2. Data Handling 
+ Sub-tissue integrated into a 25 major tissue. [Integrating Details](https://github.com/Park-JungJoon/Tissue-Classifier/blob/main/Codes/Tissue%20integrated.txt)
+ Gene randomly sampled for reduce memory consumption (19,151 -> 5,000)
+ Tissue random sampling was not performed. Skin, one of the major tissue with 1,809 samples, used as criteria tissue.
+ When tissues with a small number of samples were used, true positives were not sufficiently distributed in the validation set, which made it difficult to evaluate model performance, so tissue random sampling were not used.

<img width="700" alt="image" src="https://github.com/Park-JungJoon/Tissue-Classifier/assets/97942772/7f789ee2-6a2c-4460-b824-7365ab32c390">

  + 17,211 rows (samples), 5001 columns(gene), V1 column stands for tissue. Skin tissues are written as 1, others are 0.

</br>

## 3. Used Model - all process worked in R environment. 
+ [Logistic Regression](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/glm)
+ [LASSO](https://www.rdocumentation.org/packages/glmnet/versions/4.1-7/topics/glmnet)
+ [Tree](https://www.rdocumentation.org/packages/rpart/versions/4.1.19/topics/rpart)
+ [Random Forest](https://www.rdocumentation.org/packages/randomForest/versions/4.7-1.1/topics/randomForest)
+ [Boosting](https://www.rdocumentation.org/packages/gbm/versions/2.1.8.1/topics/gbm)

</br>

## 4. Data Description
+ All models stands for bionomial classification based on randomly sampled 5000 genes.
+ Classify standard tissue is skin, one of the most abundant tissue in GTEx dataset, with significant biologically difference.
+ The first column is 1 or 0, and if the sample is skin values are 1, all other values ​​are 0. Other cells in dataframe is a gene expression (normalized in GeTMM).
+ Random forest model take too long time, use 28 threads.

</br>

## 5. Result

Model | AUC | Binomial Deviance | Time
-|-|-|-
Linear Regression	|0.9989467	|127.34694|	6,681sec
LASSO	|0.9999938	|24.85979	|74sec	
Tree|	0.9971076|79.72645|	68sec	
Random Forest|	0.9999982|	34.31887	|31,766sec(28thr)	
Boosting	|0.9999991|	23.99763	|21,245sec

+ All model predicts well with 99% accuracy.
+ LASSO seems proper model with fast calculating.
+ Overfitting suspected.

Few genes in tree model associated with Skin.
[ENSG00000181458](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TMEM45A&keywords=ENSG00000181458)
[ENSG00000163064](https://www.genecards.org/cgi-bin/carddisp.pl?gene=EN1&keywords=ENSG00000163064#expression)
[ENSG00000156535](https://www.genecards.org/cgi-bin/carddisp.pl?gene=CD109&keywords=ENSG00000156535)

<img width="570" alt="image" src="https://github.com/Park-JungJoon/Tissue-Classifier/assets/97942772/0f4b8703-7711-4ce1-b58a-df1712065c63">


All Rdata saved as jjpark@eevee:/home/jjpark/all_model_done.RData
All codes saved as [Rmd file](https://github.com/Park-JungJoon/Tissue-Classifier/blob/main/Codes/classify_skin.Rmd).
