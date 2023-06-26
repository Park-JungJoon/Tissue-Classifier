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
+ Gene randomly sampled for reduce memory consumption (19,151 > 5,000)
+ Tissue random sampling was not performed. Skin, one of the major tissue with 1,809 samples, used as criteria tissue.
+ When tissues with a small number of samples were used, true positives were not sufficiently distributed in the validation set, which made it difficult to evaluate model performance, so tissue random sampling were not used.

<img width="700" alt="image" src="https://github.com/Park-JungJoon/Tissue-Classifier/assets/97942772/7f789ee2-6a2c-4460-b824-7365ab32c390">

  + 17,211 rows (samples), 5001 columns(gene), V1 column stands for tissue. Skin tissues are written as 1, others are 0.

## Used Model - all process worked in R environment. 
+ [Logistic Regression](https://en.wikipedia.org/wiki/Logistic_regression)
+ []
