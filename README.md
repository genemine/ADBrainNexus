# AD gene prediction
This repository provides source codes to build the model for predicting AD-associated genes and to make predictions. The model is built using ExtraTree. The feature vectors of training and test genes are extracted from the BrainNexus network, which is available at https://zenodo.org/record/5163897.

## 1. Usage
(1) After this repository is downloaded and unzipped, go into the folder. Then, make sure to unzip the mat2pred.txt.zip folder, which contains feature matrix of unlabeled genes to predict.

(2) We have created a python script, named 'AD_gene_pred.py', which includes all source codes (i) to build the model for predicting AD-associated genes and (ii) to make predictions.
Assuming that you are currently in the downloaded folder, just run the following command and you will be able to built a model and make predictions:

```bash
 python AD_gene_pred.py data/BrainNexus_AD.txt FGN 100
 
 python AD_gene_pred.py data/integrated_feature.txt integrated 100
 
 ```
In the above commands, the third input, i.e. 100, represents the times of random sampling to build sub-models: in each random sampling, a subset of negative genes (the same number as the positives) are randomly selected from the whole set of negatives. 5-fold CV is performed on the combined set of the positives and the sampled negatives. The sub-models built during CV are used to rank all human genes.

Note: The three modules, which are pandas, numpy and sklearn, need to be installed before running the script. 

(3) After running the script, the prediction results will be saved to the file named 'prediction.txt', in which each row contains a gene and the predicted probabilistic score ranging from 0 and 1 that measures the likelihood for the gene to be associated with AD.

## 2. Contact
If you have any questions, please contact us:<br>
Cuixiang Lin, lincxcsu@csu.edu.cn <br>
Jianxin Wang, jxwang@mail.csu.edu.cn<br>
