### ------------- Data Prep 

# Sequencing Run aggregate function
srx_agg <- function(x,counts="GeneCounts") {
  IDX=which(names(x) %in% "GeneCounts")
  mds<-x$MetadataSummary
  n=nrow(x[[IDX]])
  SRX_dat <- vapply(X=unique(mds$SRX_accession) ,function(srx) {
    srrs<-rownames(mds)[which(mds$SRX_accession %in% srx)]
    if (length(srrs)>1) {
      rowSums(x[[IDX]][,srrs])
    } else {
      x[[IDX]][,srrs]
    }
  } , numeric(n))
  rownames(SRX_dat) <- rownames(x[[IDX]])
  colnames(SRX_dat) <- unique(mds$SRX_accession)
  SRX_dat
}


# This function's output is a list that contains (1) a dataframe with genes 
# and its corresponding cluster; (2) the number of clusters produced; 
# and (3) the heirarchical clustering object, hr, that can be used to 
# change the number of clusters produced by adjusting cut_hmax
# Input:
# counts = normalized RNA seq data
# cut_hmax = used to optimise cluster size; denominator for h=max

cluster <- function(counts, cut_hmax){
  
  cluster_data <- list()
  
  # Hierarchical Clustering
  cl <-as.dist(1-cor(t(counts), method="spearman"))
  hr <- hclust(cl , method="complete")
  
  # optimizing the cluster size
  mycl <- cutree(hr, h=max(hr$height/cut_hmax))
  # Check the number of clusters. Can be adjusted by changing the h=max denominator
  mycl_length <- length(unique(mycl))
  
  # Prepare Cluster data frame
  clusters <- as.data.frame(mycl)
  colnames(clusters) <- "ClusterNumber"
  
  # make a new column, GeneID, from rownames
  clusters$GeneID <- rownames(clusters)
  
  cluster_data[["Clusters"]] <- clusters
  cluster_data[["Cl_length"]] <- mycl_length
  cluster_data[["hr"]] <- hr
  
  return(cluster_data)
}


# This function's output is a nested list with of the genes grouped per 
# cluster and their corresponding correlation values. 
# Inputs:
# x = logcounts (normalized gene counts from RNA seq)
# y = clusters (matrix of genes w/ cluster number)
# clust_total = total number of clusters

corr_per_clust <- function(x, y, clust_total){
  corr_cl <- list()
  for (i in 1:clust_total){
    clusterX <- y[y$ClusterNumber == i,]
    cluster_list <- as.list(clusterX$GeneID)
    cluster <- x[rownames(x) %in% cluster_list,]
    corr_result <- cor(t(cluster))
    corr_cl[[paste0("Cluster", i)]] <- corr_result
  }
  return(corr_cl)
}


# This function's output is a list of data frames containing all GO 
# terms associated with the genes belonging to a cluster. 
# Input:
# GO_annot = scerevisiae_GO_matrix_wide (matrix of genes belonging to a GO term)
# y = clusters (matrix of genes w/ cluster number)
# clust_total = total number of clusters

GO_per_cl <- function(GO_annot, clusters, clust_total){
  GO_cl <- list()
  for (i in 1:clust_total){
    clusterX <- clusters[clusters$ClusterNumber == i,]
    cluster_list <- as.list(clusterX$GeneID)
    cluster_GOterms <- GO_annot[GO_annot$GeneID %in% cluster_list,]
    rownames(cluster_GOterms)<- cluster_GOterms[,1] 
    cluster_GOterms[,1] <- c()
    cluster_GOterms <- cluster_GOterms[,which(colSums(cluster_GOterms) > 0)]
    GO_cl[[paste0("Cluster", i)]] <- cluster_GOterms
  }
  return(GO_cl)
}


### ------------- Cluster Analysis

# This function's output is a nested list with of the genes grouped 
# per cluster and their corresponding correlation values. 
# Inputs:
# x = logcounts (normalized gene counts from RNA seq)
# y = clusters (matrix of genes w/ cluster number)
# clust_total = total number of clusters

corr_per_clust <- function(x, y, clust_total){
  corr_cl <- list()
  for (i in 1:clust_total){
    clusterX <- y[y$ClusterNumber == i,]
    cluster_list <- as.list(clusterX$GeneID)
    cluster <- x[rownames(x) %in% cluster_list,]
    corr_result <- cor(t(cluster))
    corr_cl[[paste0("Cluster", i)]] <- corr_result
  }
  return(corr_cl)
}


# This function's output is a list of data frames clustered genes 
# with different total clusters. 
# Input:
# hr = heirarchical clustering object
# min and max = minimum and maximum value of the denominator used
#               in the cutree function to determine total clusters
# interval = intervals for the vector of values between min and max
#               used as input for different cluster totals

cl_lengthCut <- function(hr, min, max, interval){
  
  mycl_cuts <- list()
  
  cut <- seq(from = min, to = max, by = interval)
  
  for (i in cut){
    mycl_length <- length(unique(cutree(hr, h=max(hr$height/i))))
    
    # Prepare Cluster data frame
    clusters <- as.data.frame(mycl)
    colnames(clusters) <- "ClusterNumber"
    clusters$GeneID <- rownames(clusters)
    
    # Tally the total GeneIDs assigned per cluster
    a <- as.numeric(mycl_length)
    tally <- list()
    for (j in 1:a) {
      
      tot <- length(which(clusters$ClusterNumber == j))
      tally[[j]] <- tot
    }
    
    tallyDF <- t(as.data.frame(tally))
    rownames(tallyDF) <- c(1:a)
    
    mycl_cuts[[paste0(mycl_length)]][[paste0("GeneID_assignments")]] <- clusters
    mycl_cuts[[paste0(mycl_length)]][[paste0("Cut_value")]] <- i
    mycl_cuts[[paste0(mycl_length)]][["Tally"]] <- tallyDF
  }
  return(mycl_cuts)
}


# This function's output is a list of data frames containing all 
# GO terms associated with the genes belonging to a cluster. 
# Input:
# x = scerevisiae_GO_matrix_wide (matrix of genes belonging to a GO term)
# y = clusters (matrix of genes w/ cluster number)
# clust_total = total number of clusters

GO_per_cl <- function(x,y,clust_total){
  GO_cl <- list()
  for (i in 1:clust_total){
    clusterX <- y[y$ClusterNumber == i,]
    cluster_list <- as.list(clusterX$GeneID)
    cluster_GOterms <- x[x$GeneID %in% cluster_list,]
    rownames(cluster_GOterms)<- cluster_GOterms[,1] 
    cluster_GOterms[,1] <- c()
    cluster_GOterms <- cluster_GOterms[,which(colSums(cluster_GOterms) > 0)]
    GO_cl[[paste0("Cluster", i)]] <- cluster_GOterms
  }
  return(GO_cl)
}


# This function's output is a nested list of GO terms (column names) grouped 
# per cluster from the original input data frame before blinding
# Input:
# x = list object containing the input data frame before blinding
# clust_total = total number of clusters

GO_list_perCl <- function(x,clust_total){
  GO_names <- list()
  
  for (i in 1:clust_total){
    input <- x[[i]][["Input"]]
    GO_list <- colnames(input)
    
    GO_names[[paste0("Cluster", i)]] <- GO_list
  }
  return(GO_names)
}


# This is a derivative of the GO_list_perCl fuction modified to accomodate
# The different input data structure of the blinded genes
# This function's output is a nested list of GO terms (column names) grouped 
# per cluster from the original input data frame before blinding
# Input:
# x = list object containing the input data frame before blinding
# clust_total = total number of clusters

GO_list_perClBlind <- function(x,clust_total){
  GO_names <- list()
  
  for (i in 1:clust_total){
    input <- x[[i]]
    GO_list <- colnames(input)
    
    GO_names[[paste0("Cluster", i)]] <- GO_list
  }
  return(GO_names)
}



# This function's output is a list of data frames containing all 
# GO terms (from the blinded db) associated with the genes belonging 
# to a cluster. 
# Input:
# x = Annotation matrix after blinding
# y = clusters (matrix of genes w/ cluster number)
# GO_list_perCl = a nested list object containing the GO terms (column names)
#               of the original db (before blinding) in each cluster
# clust_total = total number of clusters

GO_per_cl_blinded <- function(x,y,GO_list_perCl,clust_total){
  GO_cl <- list()
  for (i in 1:clust_total){
    clusterX <- y[y$ClusterNumber == i,]
    cluster_list <- as.list(clusterX$GeneID)
    cluster_GOterms <- x[x$GeneID %in% cluster_list,]
    rownames(cluster_GOterms)<- cluster_GOterms[,1] 
    cluster_GOterms[,1] <- c()
    cluster_GOterms <- cluster_GOterms[, colnames(cluster_GOterms) %in% GO_list_perCl[[i]] ]
    GO_cl[[paste0("Cluster", i)]] <- cluster_GOterms
  }
  return(GO_cl)
}


# This function's output is a nested list of genes belonging to a 
# cluster with their corresponding weighted correlations. 
# Input:
# gene_list = list of genes to be correlated 
# corr_cl = correlation values for cluster
# cl_GO = GO terms for cluster
# cl_num = cluster number

wcorr_cluster <- function(gene_list, corr_cl, cl_GO){
  wcorr_result <- list()
  
  for (gene in gene_list){
    corr_value_df <- as.data.frame(corr_cl[gene,])
    weighted_go <- merge(x = corr_value_df, y = cl_GO, by = "row.names", all.x = TRUE)
    rownames(weighted_go) <- weighted_go[,1]
    weighted_go[,1] <- c()
    weighted_go[is.na(weighted_go)] <- 0
    corr_colsum <- sum(weighted_go[,1])
    weighted_go <- weighted_go[,1]*weighted_go[,2:ncol(weighted_go)]
    normalized_weighted_go <- colSums(weighted_go)/corr_colsum
    wcorr_result[[gene]] <- normalized_weighted_go
  }
  setNames(wcorr_result, paste0(gene))
  return(wcorr_result)
}



### ------------- Imputation function

# This function's output is a nested list of gene clusters 
# containing the data generated for cluster analysis. 
# Input:
# cl_GOall = GO terms per cluster
# corr_clAll = correlation values grouped by cluster
# corr_cl = correlations per cluster
# clust_total = total number of clusters
# thresh = threshold value from 0 to 1

impute <- function (cl_GOall, corr_clAll, clust_total, thresh){
  wGO_list <- list()
  
  for (i in 1:clust_total){
    
    GO_list <- rownames(cl_GOall[[i]])
    
    #If all the genes in the cluster have no GOs and the matrix is empty
    if (is_empty(GO_list) == TRUE){
      wGO_list[[paste0("Cluster", i)]][[paste0("Comment")]] <- "No Gene Ontologies found"
      wGO_list[[paste0("Cluster", i)]][[paste0("Input")]] <- 0
      wGO_list[[paste0("Cluster", i)]][[paste0("Output")]] <- 0
      wGO_list[[paste0("Cluster", i)]][[paste0("Diff")]] <- 0
      
      next
    }
    
    GO_list <- GO_list[order(GO_list)]
    
    corr_cl <- corr_clAll[[i]]
    corr_cl <- corr_cl[order(rownames(corr_cl)),]
    corr_cl <- corr_cl[,order(colnames(corr_cl))]
    gene_list <- rownames(corr_cl)
    
    cl_go <- cl_GOall[[i]]
    
    #If all the genes in the cluster have no GOs and the matrix is empty
    if (is_empty(cl_go) == TRUE){
      wGO_list[[paste0("Cluster", i)]][[paste0("Comment")]] <- "No Gene Ontologies found"
      wGO_list[[paste0("Cluster", i)]][[paste0("Input")]] <- 0
      wGO_list[[paste0("Cluster", i)]][[paste0("Output")]] <- 0
      wGO_list[[paste0("Cluster", i)]][[paste0("Diff")]] <- 0
      
      next
    }
    
    # should add gene that are not included in the GO data otherwise merge 
    # will yield weird results
    diff_go_genes <- setdiff(gene_list,rownames(cl_go))
    add <- data.frame(matrix(0,nrow=length(diff_go_genes),ncol=ncol(cl_go)))
    rownames(add) <- diff_go_genes
    colnames(add) <- colnames(cl_go)
    cl_go <- rbind(cl_go,add)
    cl_go <- cl_go[order(rownames(cl_go)),]
    cl_go <- cl_go[,order(colnames(cl_go))]
    
    wGO_cl <- wcorr_cluster(gene_list=gene_list, corr_cl=corr_cl, cl_GO=cl_go)
    wGO_df <- as.data.frame(do.call(rbind, wGO_cl))
    wGO_thresh <- (as.matrix(wGO_df) > thresh)*1 
    wGO_thresh <- wGO_thresh[order(rownames(wGO_thresh)),]
    wGO_thresh <- wGO_thresh[,order(colnames(wGO_thresh))]
    
    cl_subtract <- wGO_thresh - cl_go
    # need to flip some of the values where NAs were replaced with zeros.
    cl_subtract[cl_subtract<0] <- 1
    
    wGO_list[[paste0("Cluster", i)]][[paste0("Input")]] <- cl_go
    wGO_list[[paste0("Cluster", i)]][[paste0("Output")]] <- wGO_thresh
    wGO_list[[paste0("Cluster", i)]][[paste0("Diff")]] <- cl_subtract
    
  }
  return(wGO_list)
}


# This function otputs a nested list of clusters containing heatmaps for 
# input, output, and output-minus-input binary data frames generated by 
# the imputation funtion
# Input:
# clust_total = total number of clusters you want to evaluate
# wGO_allCl = a list object from the impute function containing the 
# input, output, and output-minus-input data frames

impute_viz <- function(clust_total, wGO_allCl){
  
  wGO_heatmaps <- list()
  
  for (i in 1:clust_total){
    
    par(cex.main = 0.5)
    colfunc <- colorRampPalette(c("white", "red"))
    
    # create heatmap for input df
    wGO_heatmaps[[paste0("Cluster", i)]][[paste0("Input")]] <- 
      heatmap.2(as.matrix(wGO_allCl[[i]][["Input"]]), main=paste0("S.cerevisiae Cluster ",
                                                                  i, " GO Terms"), scale="none", col = colfunc(25), trace="none", margins = 
                  c(5,5))
    
    # create heatmap for output df
    wGO_heatmaps[[paste0("Cluster", i)]][[paste0("Output")]] <- 
      heatmap.2(as.matrix(wGO_allCl[[i]][["Output"]]), main=paste0("S.cerevisiae 
                 Cluster ", i, " Imputed GO Terms"), scale="none", col = colfunc(25), 
                trace="none", margins = c(5,5))
    
    # create heatmap for subtraction df
    wGO_heatmaps[[paste0("Cluster", i)]][[paste0("Diff")]] <- 
      heatmap.2(as.matrix(wGO_allCl[[i]][["Diff"]]), main=paste0("S.cerevisiae Cluster ", 
                                                                 i, " GO Terms vs Imputed GO Terms"), scale="none", col = colfunc(25), 
                trace="none", margins = c(5,5))
  }
  return(wGO_heatmaps)
}


### ------------- Measures of performance 

# This function's output is a nested list of vectors that measures the 
# performance of the imputation using the blinded and original data frames
# grouped per cluster.
# Inputs:
# input_df = data frame, the original GO annotation binary matrix
# blinded_df = data frame, the imputed binary matrix using blinded data
# clust_total = total number of clusters

stats_cl <- function(input_df, blinded_df, clust_total){
  stats_list <- list()
  
  for (i in 1:clust_total){
    input <- input_df[[i]][["Input"]]
    blind <- blinded_df[[i]][["Diff"]]
    
    diff <- input - blind 
    sum <- input + blind
    
    TP <- length(which(sum == 2)) #True positive
    TN <- length(which(sum == 0)) #True negative
    FP <- length(which(diff == -1)) #False positive
    FN <- length(which(diff == 1)) #False negative
    
    # False negative rate (FNR)
    FNR <- FN/(FN+TP)
    # Sensitivity/True Positive Rate (TPR)
    TPR <- 1 - FNR
    # False Positive Rate (FPR)
    FPR <- FP/(FP+TN)
    # Specificity/True Negative Rate (TNR)
    TNR <- 1 - FPR
    # Precision/Positive Predictive Value (PPV)
    PPV <- TP/(TP+FP)
    # accuracy (ACC)
    ACC <- (TP+TN)/(TP+TN+FP+FN)
    # F1 score (is the harmonic mean of precision and sensitivity)
    F1 <- (2*TP)/((2*TP)+FP+FN)
    
    stats_list[[paste0("Cluster", i)]][[paste0("TP")]] <- TP
    stats_list[[paste0("Cluster", i)]][[paste0("TN")]] <- TN
    stats_list[[paste0("Cluster", i)]][[paste0("FP")]] <- FP
    stats_list[[paste0("Cluster", i)]][[paste0("FN")]] <- FN
    
    stats_list[[paste0("Cluster", i)]][[paste0("Sensitivity(TPR)")]] <- TPR
    stats_list[[paste0("Cluster", i)]][[paste0("Specificity(TNR)")]] <- TNR
    stats_list[[paste0("Cluster", i)]][[paste0("Precision(PPV)")]] <- PPV
    stats_list[[paste0("Cluster", i)]][[paste0("Accuracy(ACC)")]] <- ACC
    stats_list[[paste0("Cluster", i)]][[paste0("F1_Score")]] <- F1
  }
  return(stats_list)
}


# This function's output is a list of vectors that measures the 
# performance of the whole imputation using the blinded 
# and original data frames.
# Input:
# stats_cl = list, statistical measures of performance per cluster

stats_all <- function(stats_cl){
  
  stats_total <- list()
  
  # Transform nested list to data frame
  stats_df <- as.data.frame(t(as.data.frame(
    lapply(stats_cl, unlist))))
  
  TP <- sum(stats_df$TP) #True positive
  TN <- sum(stats_df$TN) #True negative
  FP <- sum(stats_df$FP) #False positive
  FN <- sum(stats_df$FN) #False negative
  
  # False negative rate (FNR)
  FNR <- FN/(FN+TP)
  # Sensitivity/True Positive Rate (TPR)
  TPR <- 1 - FNR
  # False Positive Rate (FPR)
  FPR <- FP/(FP+TN)
  # Specificity/True Negative Rate (TNR)
  TNR <- 1 - FPR
  # Precision/Positive Predictive Value (PPV)
  PPV <- TP/(TP+FP)
  # Negative Predictive Value (NPV)
  NPV <- TN/(TN+FN)
  # accuracy (ACC)
  ACC <- (TP+TN)/(TP+TN+FP+FN)
  # balanced accuracy (BA)
  BA <- (TPR+TNR)/2
  # F1 score (is the harmonic mean of precision and sensitivity)
  F1 <- (2*TP)/((2*TP)+FP+FN)
  
  
  stats_total[["Stats_df"]] <- stats_df
  
  stats_total[["TP_all"]] <- TP
  stats_total[["TN_all"]] <- TN
  stats_total[["FP_all"]] <- FP
  stats_total[["FN_all"]] <- FN
  
  stats_total[["Sensitivity"]] <- TPR
  stats_total[["Specificity"]] <- TNR
  stats_total[["Precision"]] <- PPV
  stats_total[["Accuracy"]] <- ACC
  stats_total[["F1_Score"]] <- F1
  
  return(stats_total)
}


# This function's output is a list of the measures of 
# performance from the stats_all function applied to a k-fold 
# cross validation process
# Input
# n = number of folds
# GO_annot = original GO annotation matrix
# clusters = matrix of genes w/ cluster number
# GO_list_perCl = a nested list object containing the GO 
#               terms (column names) of the original db
#               (before blinding) in each cluster
# corr_clAll = correlation values grouped by cluster
# wGO_allClusters = A nested list of gene clusters containing 
#               the data generated from cluster analysis using
#               the impute function
# clust_total = total number of clusters
# thresh = threshold value from 0 to 1

cross_val <- function(n, GO_annot, clusters,
                      GO_list_perCl, corr_clAll,
                      wGO_allClusters, clust_total, 
                      thresh){
  
  stats <- list()
  
  # partition the data
  nfolds <- n
  set.seed(42)
  folds <- sample(1:nfolds, nrow(GO_annot), replace = TRUE)
  
  dict_folds <- data.frame(GO_annot$GeneID, folds)
  saveRDS(dict_folds, "dict_folds.rds")
  
  for (i in 1:n){
    which_fold <- i
    
    ind_blinded <- which(dict_folds$folds == which_fold)
    GO_blinded <- GO_annot
    GO_blinded[ind_blinded, 2:ncol(GO_blinded)] = 0
    
    # Get the GO terms by cluster using the blinded data 
    # and the GO list per cluster
    GO_blindedCl <- GO_per_cl_blinded(GO_blinded, clusters, 
                                      GO_list_perCl, clust_total)
    
    # Impute function
    wGO_blinded <- impute(GO_blindedCl, corr_clAll,
                          clust_total, thresh)
    
    # stats
    stats_perCl <- stats_cl(wGO_allClusters, wGO_blinded, clust_total)
    stats_final <- stats_all(stats_perCl)
    
    stats[[paste0("Fold", i)]] <- stats_final
  }
  return(stats)
}


### ------------- Optimisation 

# This function's output is a list of values to use
# as the denominator for hmax in the cutree function 
# to get a desired number of the total clusters
# Inputs:
# hr = object from the hclust function
# cut_min = hmax denominator value for cl_min (can start with 0)
# cl_min, cl_max = min and max of total clusters
# interval = intervals beween min and max total clusters

mycl_opt <- function(hr, cut_min, 
                     cl_min, cl_max, interval){
  
  mycl_cuts <- list()
  cl_total <- seq(from = cl_min, to = cl_max, by = interval)
  
  # optimizing the cluster size
  for (i in cl_total){
    
    cut = cut_min
    mycl_length = 0
    while (mycl_length!=i) {
      mycl_length <- length(unique(cutree(hr, h=max(hr$height/cut))))
      
      if (mycl_length==i){
        mycl_cuts[[paste0(i)]] <- cut
      } else cut = cut + interval           
    }
  }
  return(mycl_cuts)
}