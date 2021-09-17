# Table of Contents

[Project Files](#project-files)

1. [Data Preparation](#data-preparation)</br>
    1A. [RNASeq Data](#rnaseq-data)
      a. [Quality Control](#quality-control)
      b. [Aggregating SRX to SRR](#aggregating)
      c. [Normalisation](#normalisation)
      d. [Clustering](#clustering)
      
    1B. [GO Annotation Data](#go-annotation-data)
      a. [Binary clasiffication](#binary-classification)
      b. [Training and Testing Data Split](#training-and-testing-data-split)
2. [Imputation](#imputation)</br>
    2A. [Correlation per cluster](#correlation-per-cluster)</br>
    2B. [Assigning weight to GO terms](#assigning-weight-to-GO-terms)</br>
    2C. [Assigning a threshold](#assigning-a-threshold)</br>
    2D. [Imputation process](#imputation-process)
3. [Optimisation](#optimisation)</br>
    3A. [Parameters to be optimised](#parameters-to-be-optimised)</br>
    3B. [Kfold validation](#kfold-validation)</br>
    3C. [Prediction Scores](#prediction-scores)</br>
      a. [Confusion Matrix](#confusion-matrix)</br>
4. [Functions](#functions)</br>
    4A. [Data preparation functions](#data-preparation-functions)</br>
    4B. [Imputation functions](#imputation-functions)</br>
    4C. [Optimisation functions](#optimisation-functions)</br>

</br></br>

## Project Files

How to run the algorithm:

  1. Create a folder named "Data" (case sensitive)
  2. Open and run YeastDataPrep.Rmd 
  3. Open and run YeastImputation-Blinded.Rmd 

Note that the `optimise_impute` function which performs the coarse sweep and kfold validation on Chunk 1 of YeastImputation-Blinded.Rmd have run times that last up to 6 hours because of Cluster 2 having ~2000 gene members. 

For a faster runtime, search "runtime" on functions.r and un-comment the line right below and run the `source("functions.R")` line.

## 1. Data Preparation
Data preparation is divided into two main processes: RNASeq counts data 
processing and Gene Ontology Annotation data processing. The processed data 
are saved as .rds files to be used in the downstream processes.</br>

### 1A. RNASeq Data
Saccharomyces cerevisiae count and quality metrics data were downloaded from DEE2

  - counts: http://dee2.io/mx/scerevisiae_se.tsv.bz2
  - qc: http://dee2.io/mx/scerevisiae_qc.tsv.bz2)
  
Downloaded count data in long format is converted to wide format for easier matrix manipulation.

#### *1Aa. Quality Control*
The quality data allows for filtering high quality datasets. 
Databases were filtered and quantified according to QC_SUMMARY = PASS, WARN, or FAIL. 
For all the downstream processes, only datasets tagged "PASSED" are used. 

#### *1Ab. Aggregating SRR to SRX*
Sequence Read Archive (SRA) Accession Codes are also provided in the quality metrics data. These codes describes the source and type of the dataset.
</br>
Accession code format: aRbx (Examples: SRR#, SRX#, ERR#, DRR#, etc.)</br>
Where:</br>
**First letter (represented by 'a' above)** Represents the source database to which the sample was originally uploaded, before being synchronised with the other 2 databases:</br>

     S – NCBI’s SRA database</br>
     E – EBI’s database</br>
     D – DDBJ database</br>
**Second Letter**</br>

    Always the letter R</br>
**Third letter (represented by 'b' above)** Type of data represented:</br>

    R – Run</br>
    X – Experiment</br>
    S – Sample</br>
    P – Project / study</br>
**x - unique accession code**
</br>
Some datasets have multiple runs (SRR) of the same experiment (SRX).
Multiple runs, if any, were aggregated to its corresponding experiment using the srx_agg function.  

#### *1Ac. Normalisation*
Library size bias was accounted for via TMM normalisation using the `DESeq2` library and then scaled by R's base function `scale()`.

#### *1Ad. Clustering*
After data quality control and normalisation processes, the counts data is clustered using the **hierarchical agglomerative cluster analysis (HAC)** with the [`hclust` function](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust)
	
`hclust` uses **[Agglomerative clustering](https://uc-r.github.io/hc_clustering#algorithms)** also known as AGNES (Agglomerative Nesting) that outputs a tree which can be plotted as a dendrogram. It is good at identifying small clusters.

Spearman correlation was used to compute the distance matrix 

The **[complete (maximum) linkage clustering](https://uc-r.github.io/hc_clustering#algorithms)** was used to measure the dissimilarity between two clusters of observations. This is also referred to as the linkage method.

    Note: Maximum or complete linkage clustering computes all pairwise dissimilarities between the elements in cluster 1 and the elements in cluster 2, and considers the largest value (i.e., maximum value) of these dissimilarities as the distance between the two clusters. It tends to produce more compact clusters which is ideal in determining close relationships between gene expressions. (Reference is the link above).
	
    Note2: The more tight, dense are clusters inside and the less density is outside of them (or the wider apart are the clusters) - the greater is the internal validity.

[Reference](https://stats.stackexchange.com/questions/195456/how-to-select-a-clustering-method-how-to-validate-a-cluster-solution-to-warran/195481#195481)

A series of cluster sizes (total clusters or K) from 50 to 2000 is passed to the `cl_lengthCut` function which outputs a list all the gene IDs tagged to belong to a cluster. 
</br>

### 1B. GO Annotation Data
The [Gene Ontology Annotation File (GAF) format](http://current.geneontology.org/annotations/sgd.gaf.gz) of yeast was downloaded from geneontology.org. 

#### *1Ba. Binary clasiffication*
A binary matrix was constructed with the Entrez Gene IDs (NCBI gene IDs) as rownames and GO IDs as the column names (wide format). If a GO ID is associated with a Gene ID, the cell will equal to 1, otherwise it will be zero. 
</br>
There were Gene IDs present in the GAF that is not present in the counts data. These were removed. 

#### *1Bb. Training and Testing Data Split*
Ten per cent of the gene IDs from the binary GO matrix will be set aside as testing data and the remaining 90 per cent will be used as training data for the algorithm.
</br>
</br>

## 2. Imputation
**Cluster Level Analysis**

>	This explanation focuses in one cluster. This algorithm will be repeated for the rest of the clusters.

### 2A. Correlation per cluster
The correlation of the count data for all Gene IDs belonging to cluster 1 is computed using the Spearman method. This produces an n x n matrix with rows and columns as Gene IDs. 

### 2B. Assigning weight to GO terms in the cluster
The correlation matrix is turned into an edge list using the `igraph` library. An edge list contains pairs of genes with their corresponding correlation value. Each gene ID is filtered from this list along with the gene IDs it has paired with. The correlation value will be treated as the weight of the gene to gene connection resulting in a network for each cluster.

The GO terms for the cluster is filtered and these values will be multiplied across the 1s and 0s of each GO ID column. 

### 2C. Assigning a threshold
After assigning the weights, the values will averaged (*still figuring out the optimal method*) and a threshold will determine a cutoff value for the imputation process.

### 2D. Imputation process
All values euqal to and above the threshold value will be assigned as 1 (meaning the gene is annotated with the GO ID) and all values below will be assigned 0 (meaning the gene is does not belong to the GO ID). This will be repeated with all the gene IDs until a new GO annotation matrix with the same dimensions as the original one but with imputed values is created. 
</br>
</br>

## 3. Optimisation
The imputation process will then be optimised using the techniques below.

### 3A. Parameters to be optimised
A series of cluster sizes and the threshold values is used to find the optimal parameters for imputation.

### 3B. Kfold validation
Kfold cross validation at K = 10 is used to see if the parameters used yield optimal values. The GO annotation for a cluster is divided into 10 parts. Part 1 of the gene IDs will have their annotations set to zero and this dataset will be used to run the first fold of the cross validation. Parts 2 to 10 will be used in the same manner and the process is repeated until it reaches 10 folds.

### 3C. Prediction Scores and the Confusion Matrix
To measure the performance of each parameter pair, Accuracy, Precision, Recall, and the F1 score will be measured using a confusion matrix. 

#### *Confusion Matrix*
- Actual values = the binary GO matrix filtered for the cluster 
- Predicted values = the GO matrix with the imputed values after running the algorithm.
</br>

  |             | **Predicted 1**         | **Predicted 0**    |    
  |:-----------:|:-----------:            |:-----------:       |
  |**Actual 1** | True Postivie (TP)      | False Negative (FN)| 
  |**Actual 0** | False Positive (FP)     | True Negative (TN) |
  
  
  To solve for TP and TN the following equation is used:</br>
  TP or TN = Actual + Predicted

  | **Addition**     | **Predicted 1**         | **Predicted 0**         |
  | :--------------: | :-------------:         | :-------------:         |
  | **Actual 1**     | $\color{red}{\text{2}}$ | 1                       |
  | **Actual 0**     | 1                       | $\color{red}{\text{0}}$ |


  To solve for FP and FN the following equation is used:</br>
  FP or FN = Actual - Predicted

  | **Subtraction** | **Predicted 1**          | **Predicted 0**         |
  |:---------------:|:---------------:         |:---------------:        |
  | **Actual 1**    | 0                        | $\color{red}{\text{1}}$ |
  | **Actual 0**    | $\color{red}{\text{-1}}$ | 0                       |


  Therefore, via matrix addition and subtraction: 

    ```
        diff <- Actual - Predicted 
        sum <- Actual  + Predicted
        
        TP <- length(which(sum == 2)) #True positive
        TN <- length(which(sum == 0)) #True negative
        FP <- length(which(diff == -1)) #False positive
        FN <- length(which(diff == 1)) #False negative
        
        # Precision/Positive Predictive Value (PPV)
        PPV <- TP/(TP+FP)
        # Accuracy (ACC)
        ACC <- (TP+TN)/(TP+TN+FP+FN)
        # F1 score (is the harmonic mean of precision and sensitivity)
        F1 <- (2*TP)/((2*TP)+FP+FN)
        # Recall
        Recall <- TP/(TP+FN)
    ```

## 4. Functions 
### 4A. Data preparation functions

**srx_agg** - Multiple runs, if any, are aggregated to its corresponding experiment. Function from https://github.com/markziemann/getDEE2/issues/7#issuecomment-639290203 

</br>
**cl_lengthCut** - This function's output is a list  grouped by total clusters which gives (1) a table of GeneIDs assigned to a cluster number, (2) the cutree denominator to achieve that total cluster, and (3) the tally of the total Gene IDs belonging to a cluster number.

- Input:
    - hr = heirarchical clustering object from `hclust()` of the stats library
    - min and max = minimum and maximum value of the denominator used in the cutree function to determine total clusters
    - interval = intervals for the vector of values between min and max used as input for          different cluster totals

### 4B. Imputation functions

**wcorr_cluster** - assigns "weights" to the GO matrix of a cluster by multiplying all the correlation values of a gene from the `edge_list()` output; it then takes a weighted average of the matrix to be imputed later with a threshold value. This function's output is a nested list of genes belonging to a cluster with a tally of how many correlation values passed the threshold
that connects the gene to a GO Term.

- Input:
    - gene_list = list of genes to be correlated 
    - edgeList = an edge list (from igraph Library) from `edge_list()` 
    - cl_GO = a list of GO terms for a cluster
    - thresh = threshold from 0 to 1

</br>
**impute** - takes the output from `wcorr_cluster()`, applies the threshold, and assigns 1 to
values equal to or above the cut and zero if otherwise. This produces a binary matrix similar to the original GO annotation matrix but the assignment of 1s and 0s is created by the algorithm: genes that belongs to a GO ID is marked 1 and 0 if not.
This function's output is a nested list of binary matrices for (1) the original GO annotations of genes for a cluster, (2) imputed GO binary matrix, (3) blinded GO matrix from `cross_val()`, and (4) a matrix from matrix 4 subtracted by matrix 1

- Input:
    - GOterms_perCl = GO terms per cluster from `GO_per_cl()`
    - cl_GOall = GO terms per cluster from `GO_per_cl_blinded()`
    - corr_clAll = correlation values grouped by cluster
    - clust_total = total number of clusters
    - cor_edge_list = an edge list (from igraph Library) from the correlation values per               cluster
    - thresh = threshold value from 0 to 1
- Dependencies:
    - `wcorr_cluster`


### 4C. Optimisation functions

**optimise_impute** - takes a series of cluster sizes and thresholds and applies each pair to a cross validation process. This function outputs a nested list of prediction scores from a 10-fold validation process. Scores are grouped by a parameter pair consisting of a total cluster and a threshold value. Each threshold value will be applied to each total cluster value.

- Inputs:
    - cl_list = a list the target total clusters (ex: c(20,100,15000))
    - thresh_list = a numerical list denoting the target threshold value from 0 to 1
    - cuttree_values = the output of `cl_lengthCut()` which gives a table of GeneIDs assigned         to clusters
    - counts = normalized gene counts from RNA seq
    - GO_annot = the Gene Ontology binary matrix for all GeneIDs
- Dependencies:
    - `corr_per_clust`
    - `GO_per_cl`
    - `GO_per_cl_list`
    - `cross_val`

</br>
**corr_per_clust** - takes the counts data of all genes in a cluster to make a correlation matrix. This function's output is a nested list of the genes grouped per cluster and their corresponding correlation matrices. 

- Inputs:
    - x = counts (normalized gene counts from RNA seq)
    - y = matrix of genes w/ cluster number from `cl_lengthCut()`
    - clust_total = total number of clusters

</br>
**GO_per_cl** - subsets the bigger GO binary annotation matrix into a smaller matrix that only contains gene IDs and corresponding GO IDs of a cluster. This function's output is a list of data frames containing all GO terms associated with the genes belonging to a particular cluster.

- Input:
    - x = binary matrix of genes belonging to a GO term
    - y = y = matrix of genes w/ cluster number from `cl_lengthCut()`
    - clust_total = total number of clusters

</br>
**GO_per_cl_list** - lists down all the GO IDs of a cluster; filters out all columns (GO IDs) that sums up to zero to reduce computing load. This function's output is a nested list of GO terms (column names) grouped per cluster from `GO_per_cl()`

- Input:
    - x = list object containing GO terms per cluster from `GO_per_cl()`
    - clust_total = total number of clusters

</br>
**mean_Csweep** - This function's output is a dataframe of the average value of a specific measure of performance from all cross validation folds

- Inputs:
    - kfold_list = a daframe containing all the measures of performance for each fold from                     `cross_val()`
    - stat_type = an index number from the `cross_val()` output list indicating the measure                    of performance to be averaged from 2 to 11 in the ff order: Total Positive                   (TP), Total Negative (TN), False Positive (FP), False Negative (FN),                         Sensitivity (TPR), Specificity (TNR), Precision (PPV), F1 Score (F1), and                    Recall
    
 </br>   
**summary_Csweep** - gets all the average values from `mean_Csweep` and creates a summary. This function's output is a dataframe.

- Inputs:
    - kfold_list = a daframe containing all the measures of performance for each fold from                     `cross_val()`
- Dependencies:
    - `mean_Csweep` 

</br>
**cross_val** - performs a kfold cross validation. This function's output is a list of the overall performance measures of all folds from `stats_all()` and a dataframe containing the genes assigned to each fold.

- Input: 
    - n = number of folds
    - GO_annot = original GO annotation matrix
    - clusters = matrix of genes w/ cluster number
    - GOterms_perCl = GO terms per cluster from `GO_per_cl()`
    - GO_list_perCl = a nested list object containing the GO terms (column names) of the           original db (before blinding) in each cluster
    - corr_clAll = correlation values grouped by cluster
    - clust_total = total number of clusters
    - thresh = threshold value from 0 to 1
- Dependencies
    - `GO_per_cl_blinded`
    - `edge_list`
    - `impute`
    - `stats_cl`
    - `stats_all`
    
</br>
**GO_per_cl_blinded** - a derivative of `GO_per_cl`, this function also subsets from the a bigger blinded GO annotation matrix during cross validation. Instead of taking out zero sum columns, this function makes sure that the column names of the original GO matrix subset is equal to the blinded subset for dimensional consistency. 
This function's output is a list of data frames containing all GO terms (from the blinded db) associated with the genes belonging to a cluster.

- Input:
    - x = Annotation matrix after blinding
    - y = clusters (matrix of genes w/ cluster number)
    - GO_list_perCl = a nested list object containing the GO terms (column names) of the           original db (before blinding) in each cluster
    - clust_total = total number of clusters
- Dependencies
    `GO_list_perCl`

</br>
**edge_list** - transforms all the correlation martrices of each cluster into a an edge list.
This function's output is a list that contains a dataframe with the pair of genes that makes up an edge (listed from and to) with a corresponding correlation value treated as the edge 
weight. This list is created using the igraph library.

- Input:
    - corr_allCl = a list of correlation matrices per cluster from `corr_per_clust()`
    - clust_total = total number of clusters
- Dependencies
    - igraph library

</br>
**stats_cl** - calculates the performance per fold during cross validation (see `cross_val`). This function's output is a nested list of performance scores grouped per cluster.

- Inputs:
    - imputed_df = data frame, the imputed binary matrix using blinded data from                             `impute()`
    - clust_total = total number of clusters

</br>
**stats_all** - This function's output is a list of values measuring the performance of the whole imputation using the blinded and original data frames.

- Input:
    - stats_cl = list, statistical measures of performance per cluster
- Dependencies
    - `stats_cl`
