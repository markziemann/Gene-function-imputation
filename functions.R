### ------------- Data Prep 

# Sequencing Run aggregate function
# Note for Mark: kindly let me know how to cite this properly. 
# This code was from your github. 

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

# This function's output is a list  grouped by total clusters 
# which gives (1) a table of GeneIDs assigned to a cluster number, 
# (2) the cutree denominator to achieve that total cluster, and
# (3) the tally of the total Gene IDs belonging to a cluster number.
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
    
    # for complete method
    #mycl <- cutree(hr, h=max(hr$height/i))
    #for ave linking method
    mycl <- cutree(hr, k=i)
    mycl_length <- length(unique(mycl))
    
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

# This function is an improvement for cl_lengthCut().
# It uses the cutreeHybrid() from the dynamicTreeCut library
# to iterative cut a dendrogram based on a minimum cluster size.
# The output is a list  grouped by total clusters 
# which gives (1) a table of GeneIDs assigned to a cluster number, 
# (2) the cutree denominator to achieve that total cluster, and
# (3) the tally of the total Gene IDs belonging to a cluster number.
# Input:
# hr = heirarchical clustering object
# cl = distance matrix used to build hr
# min and max = minimum and maximum value of the a sequence of values used
#               in the minClusterSize attribute of the cutreeHybrid function 
# interval = intervals for the sequence of values between min and max

cl_cut_dynamic <- function(hr, cl, min, max, interval){
  
  mycl_cuts <- list()
  
  cut <- seq(from = min, to = max, by = interval)
  
  for (i in cut){
    
    dyn_tree <- cutreeHybrid(
      # Input data: basic tree cutiing
      dendro=hr, distM=as.matrix(cl),
      # Branch cut criteria and options
      minClusterSize = i, deepSplit = 4,
      # PAM stage options
      pamStage = TRUE, pamRespectsDendro = FALSE,
      respectSmallClusters = TRUE,
      # Various options
      verbose = 1, indent = 0)
    
    dyn_cl <- as.data.frame(dyn_tree$labels)
    colnames(dyn_cl) <- "ClusterNumber"
    dyn_cl$GeneID <- hr$labels
    
    # Tally the total GeneIDs assigned per cluster
    mycl_length <- length(unique(dyn_tree$labels))
    tally <- list()
    for (j in 1:mycl_length) {
      tot <- length(which(dyn_cl$ClusterNumber == j))
      tally[[j]] <- tot
    }
    
    tallyDF <- t(as.data.frame(tally))
    rownames(tallyDF) <- c(1:mycl_length)
    
    
    mycl_cuts[[paste0(mycl_length)]][[paste0("GeneID_assignments")]] <- dyn_cl
    mycl_cuts[[paste0(mycl_length)]][[paste0("min_cl_size")]] <- i
    mycl_cuts[[paste0(mycl_length)]][["Tally"]] <- tallyDF
  }
  return(mycl_cuts)
}


# This function's output is a list that contains a dataframe with 
# the pair of genes that makes up an edge (listed from and to)
# with a corresponding correlation value treated as the edge 
# weight. This list is created using the igraph library
# Input:
# corr_allCl = a list of correlation matrices per cluster
# clust_total = total number of clusters

edge_list <- function(corr_allCl, clust_total){
  e_list <- list()
  
  for(i in 1:clust_total){
    # Create igraph object
    cor_g <- graph_from_adjacency_matrix(corr_allCl[[i]], mode='directed', weighted = 'correlation', 
                                         diag = TRUE)
    # Extract edge list
    cor_edge_list <- get.data.frame(cor_g, 'edges')
    
    e_list[[paste0("Cluster", i)]] <- cor_edge_list
  }
  return(e_list)
}


### ------------- Cluster Analysis

# This function's output is a nested list of the genes grouped 
# per cluster and their corresponding correlation values. 
# Inputs:
# x = counts (normalized gene counts from RNA seq)
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


# This function's output is a list of data frames containing all 
# GO terms associated with the genes belonging to a particular cluster. 
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

# This is a derivative of the GO_list_perCl fuction modified to accomodate
# The different input data structure of the blinded genes
# This function's output is a nested list of GO terms (column names) grouped 
# per cluster from the original input data frame before blinding
# Input:
# x = list object containing the input data frame before blinding
# clust_total = total number of clusters

GO_per_cl_list <- function(x,clust_total){
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


### ------------- Imputation functions

# This function calculates the geometric mean of an array of numbers

gmean <- function(x){
  exp(mean(log(x[x>0])))
}

# This function's output is a nested list of genes belonging to a 
# cluster with a tally of how many correlation values passed the threshold
# that connects the gene to a GO Term. 
# Input:
# gene_list = list of genes to be correlated 
# edgeList = an edge list (from igraph Library) from the correlation 
#             values in a cluster
# cl_go = a list of GO terms for a cluster
# thresh = threshold from 0 to 1


wcorr_cluster <- function(gene_list, edgeList, cl_go, thresh){
  wcorr_result <- list()
  
  for (gene in gene_list){
    corr_value_df <- edgeList[edgeList$from %in% gene, c(2,3)]
    rownames(corr_value_df) <- corr_value_df$to
    weighted_go <- merge(x = corr_value_df, y = cl_go, by = "row.names", all.x = TRUE)
    rownames(weighted_go) <- weighted_go[,1]
    weighted_go[,1:2] <- c()
    # turn NA to zero for GeneIDs with no GO terms
    weighted_go[is.na(weighted_go)] <- 0
    weighted_go <- weighted_go[,1]*weighted_go[,2:ncol(weighted_go)]
    wGO_thresh <- (as.matrix(weighted_go) > thresh)*1 
    imputed_go <- colSums(weighted_go)
    #imputed_go <- colSums(weighted_go)/colSums(cl_go)
    #imputed_go <- apply(weighted_go, 2, gmean)
    #imputed_go <- apply(weighted_go, 2, max)
    wcorr_result[[gene]] <- imputed_go
  }
  setNames(wcorr_result, paste0(gene))
  return(wcorr_result)
}


# This function's output is a nested list of gene clusters 
# containing the data generated for cluster analysis. 
# Input:
# GOterms_perCl = GO terms per cluster from GO_per_cl()
# cl_GOall = GO terms per cluster from GO_per_cl_blinded()
# corr_clAll = correlation values grouped by cluster
# clust_total = total number of clusters
# cor_edge_list = an edge list (from igraph Library) from the correlation 
#           values per cluster
# thresh = threshold value from 0 to 1

impute <- function (GOterms_perCl, cl_GOall=GO_blindedCl, corr_clAll, 
                    clust_total, cor_edge_list, thresh){
  wGO_list <- list()
  
  for (i in 1:clust_total){
    print(i)
    corr_cl <- corr_clAll[[i]]
    cl_go <- cl_GOall[[i]]
    cl_go_orig <- GOterms_perCl[[i]]
    
    #If all the genes in the cluster have no GOs and the matrix is empty
    if (is_empty(corr_cl) == TRUE || isTRUE(corr_cl == 0) || 
        is_empty(cl_go_orig) == TRUE || isTRUE(cl_go_orig == 0) ){
      wGO_list[[paste0("Cluster", i)]][[paste0("Comment")]] <- "No Gene Ontologies found"
      wGO_list[[paste0("Cluster", i)]][[paste0("Input")]] <- 0
      wGO_list[[paste0("Cluster", i)]][[paste0("Output")]] <- 0
      wGO_list[[paste0("Cluster", i)]][[paste0("Blinded_GO")]] <- 0
      wGO_list[[paste0("Cluster", i)]][[paste0("Diff")]] <- 0
      
      next
    }
    
    corr_cl <- corr_cl[order(rownames(corr_cl)),]
    corr_cl <- corr_cl[,order(colnames(corr_cl))]
    gene_list <- rownames(corr_cl)
    
    # should add gene that are not included in the GO data otherwise merge 
    # will yield weird results
    diff_go_genes <- setdiff(gene_list,rownames(cl_go))
    if(length(diff_go_genes) > 0) {
      add <- data.frame(matrix(0,nrow=length(diff_go_genes),ncol=ncol(cl_go)))
      rownames(add) <- diff_go_genes
      colnames(add) <- colnames(cl_go)
      cl_go <- rbind(cl_go,add)
    }
    cl_go <- cl_go[order(rownames(cl_go)),]
    cl_go <- cl_go[,order(colnames(cl_go))]
    
    # Do the same process above for the orginal
    # annotation matrix to use for model validation
    diff_go_genes2 <- setdiff(gene_list,rownames(cl_go_orig))
    # Check and add non-annotated genes
    if(length(diff_go_genes2) > 0) {
      add <- data.frame(matrix(0,nrow=length(diff_go_genes2),ncol=ncol(cl_go_orig)))
      rownames(add) <- diff_go_genes2
      colnames(add) <- colnames(cl_go_orig)
      cl_go_orig <- rbind(cl_go_orig,add)
    }
    # Check and filter GO terms from blinded GO matrix (cl_go)
    diff_go_annot <- setdiff(colnames(cl_go_orig),colnames(cl_go))
    if(length(diff_go_annot) > 0){
      cl_go_orig <- cl_go_orig[, colnames(cl_go_orig) %in% colnames(cl_go)]
    }
    
    cl_go_orig <- cl_go_orig[order(rownames(cl_go_orig)),]
    cl_go_orig <- cl_go_orig[,order(colnames(cl_go_orig))]
    
    edgeList <- cor_edge_list[[i]]
    wGO_cl <- wcorr_cluster(gene_list=gene_list, edgeList=edgeList, 
                            cl_go=cl_go, thresh=thresh)
    wGO_df <- as.data.frame(do.call(rbind, wGO_cl))
    wGO_thresh <- (as.matrix(wGO_df) > 0)*1 
    #wGO_thresh <- (as.matrix(wGO_df) >=thresh)*1 
    wGO_thresh <- wGO_thresh[order(rownames(wGO_thresh)),]
    wGO_thresh <- wGO_thresh[,order(colnames(wGO_thresh))]
    
    cl_subtract <- wGO_thresh - cl_go
    
    wGO_list[[paste0("Cluster", i)]][[paste0("Input")]] <- cl_go_orig
    wGO_list[[paste0("Cluster", i)]][[paste0("Output")]] <- wGO_thresh
    wGO_list[[paste0("Cluster", i)]][[paste0("Blinded_GO")]] <- cl_go
    wGO_list[[paste0("Cluster", i)]][[paste0("Diff")]] <- cl_subtract
    
    
  }
  return(wGO_list)
}

### ------------- Measures of performance 

# This function's output is a nested list of vectors that measures the 
# performance of the imputation using the blinded and original data frames
# grouped per cluster.
# Inputs:
# imputed_df = data frame, the imputed binary matrix from `impute()`
# clust_total = total number of clusters

stats_cl <- function(imputed_df, clust_total){
  stats_list <- list()
  
  for (i in 1:clust_total){
    
    input <- imputed_df[[i]][["Input"]]
    blind <- imputed_df[[i]][["Output"]]
    
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
    # Accuracy (ACC)
    ACC <- (TP+TN)/(TP+TN+FP+FN)
    # F1 score (is the harmonic mean of precision and sensitivity)
    F1 <- (2*TP)/((2*TP)+FP+FN)
    # Recall
    Recall <- TP/(TP+FN)
    
    stats_list[[paste0("Cluster", i)]][[paste0("TP")]] <- TP
    stats_list[[paste0("Cluster", i)]][[paste0("TN")]] <- TN
    stats_list[[paste0("Cluster", i)]][[paste0("FP")]] <- FP
    stats_list[[paste0("Cluster", i)]][[paste0("FN")]] <- FN
    
    stats_list[[paste0("Cluster", i)]][[paste0("Sensitivity(TPR)")]] <- TPR
    stats_list[[paste0("Cluster", i)]][[paste0("Specificity(TNR)")]] <- TNR
    stats_list[[paste0("Cluster", i)]][[paste0("FPR")]] <- FPR
    stats_list[[paste0("Cluster", i)]][[paste0("FNR")]] <- FNR
    stats_list[[paste0("Cluster", i)]][[paste0("Precision(PPV)")]] <- PPV
    stats_list[[paste0("Cluster", i)]][[paste0("Accuracy(ACC)")]] <- ACC
    stats_list[[paste0("Cluster", i)]][[paste0("F1_Score")]] <- F1
    stats_list[[paste0("Cluster", i)]][[paste0("Recall")]] <- Recall
  }
  return(stats_list)
}


# This function's output is a list of values measuring the 
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
  Recall <- TP/(TP+FN)
  
  stats_total[["Stats_df"]] <- stats_df
  
  stats_total[["TP_all"]] <- TP
  stats_total[["TN_all"]] <- TN
  stats_total[["FP_all"]] <- FP
  stats_total[["FN_all"]] <- FN
  
  stats_total[["Sensitivity"]] <- TPR
  stats_total[["Specificity"]] <- TNR
  stats_total[["FPR"]] <- FPR
  stats_total[["FNR"]] <- FNR
  stats_total[["Precision"]] <- PPV
  stats_total[["Accuracy"]] <- ACC
  stats_total[["F1_Score"]] <- F1
  stats_total[["Recall"]] <- Recall
  
  return(stats_total)
}


# This function's output is a list of the measures of 
# performance from the stats_all function applied to a k-fold 
# cross validation process
# Input
# n = number of folds
# GO_annot = original GO annotation matrix
# clusters = matrix of genes w/ cluster number
# GOterms_perCl = GO terms per cluster from GO_per_cl()
# GO_list_perCl = a nested list object containing the GO 
#               terms (column names) of the original db
#               (before blinding) in each cluster
# corr_clAll = correlation values grouped by cluster
# clust_total = total number of clusters
# thresh = threshold value from 0 to 1

cross_val <- function(n, GO_annot, clusters, GOterms_perCl,
                      GO_list_perCl, corr_clAll,
                      clust_total, thresh){
  
  stats <- list()
  
  # partition the data
  nfolds <- n
  set.seed(42)
  folds <- sample(1:nfolds, nrow(GO_annot), replace = TRUE)
  
  dict_folds <- data.frame(GO_annot$GeneID, folds)
  
  for (i in 1:n){
    
    which_fold <- i
    print(paste0("Fold ", i))
    ind_blinded <- which(dict_folds$folds == which_fold)
    GO_blinded <- GO_annot
    GO_blinded[ind_blinded, 2:ncol(GO_blinded)] = 0
    
    # Get the GO terms by cluster using the blinded data 
    # and the GO list per cluster
    GO_blindedCl <- GO_per_cl_blinded(GO_blinded, clusters, 
                                      GO_list_perCl, clust_total)
    
    # Get the edge list per cluster
    cor_edge_list <- edge_list(corr_allCl=corr_clAll, clust_total)
    
    # Impute function
    wGO_blinded <- impute(GOterms_perCl, GO_blindedCl, corr_clAll,
                          clust_total, cor_edge_list, thresh)
    
    # stats
    stats_perCl <- stats_cl(wGO_blinded, clust_total)
    stats_final <- stats_all(stats_perCl)
    
    stats[[paste0("Fold", i)]] <- stats_final
  }
  stats[["Index_folds"]] <- dict_folds
  return(stats)
  
}

### ------------- Optimisation 

# This function's output is a dataframe of the averaged 
# value of a specific measure of performance from all cross 
# validation folds
# Inputs:
# kfold_list = a daframe containing all the measures of performance 
#             for each fold. See stats_all()
# stat_type = a number indicating the measure of performance to be averaged
#           from 1 to 10 in the ff order: Total Positive (TP), Total Negative
#           (TN), False Positive (FP), False Negative (FN), Sensitivity (TPR),
#           Specificity (TNR), Precision (PPV), and F1 Score (F1)

mean_Csweep <- function(kfold_list, stat_type) {
  mean_df <- data.frame()
  
  mean_df <- do.call(rbind.data.frame, lapply(kfold_list[[1]], "[", stat_type))
  colnames(mean_df)[1] <- names(kfold_list[1])
  
  for (i in 2:length(kfold_list)){
    col_name <- names(kfold_list[i])
    mean_df[col_name] <- do.call(rbind.data.frame, lapply(kfold_list[[i]], "[", stat_type))
  }
  
  Mean <- colMeans(mean_df)
  mean_df <- rbind(mean_df, Mean)
  rownames(mean_df)[rownames(mean_df) == "11"] <- "Mean"
  
  return(mean_df)
}


# This function's output is a dataframe of the averaged 
# value of all measures of performances for all folds from stats_all()
# Inputs:
# kfold_list = a daframe containing all the measures of performance 
#             for each fold. See stats_all()


summary_Csweep <- function(kfold_list){
  summary_df <- data.frame()
  summary <- mean_Csweep(kfold_list, stat_type=6)
  row_name <- names(kfold_list[[1]][[1]][6])
  summary_df <- summary[11,]
  rownames(summary_df)[rownames(summary_df) == "Mean"] <- row_name
  
  stat_type_list <- c(7:length(kfold_list[[1]][[1]]))
  
  for (i in stat_type_list){
    summary <- mean_Csweep(kfold_list, stat_type=i)
    row_name <- names(kfold_list[[1]][[1]][i])
    summary_df[nrow(summary_df)+1,] <- summary[11,]
    rownames(summary_df)[rownames(summary_df) == "Mean"] <- row_name
  }
  return(summary_df)
}


# This function outputs a nested list of prediction scores from 
# a 10-fold validation process. Scores are grouped by a parameter pair
# consisting of a total cluster and a threshold value. Each threshold 
# value will be applied to each total cluster value.
# Inputs:
# cl_list = a list the target total clusters (ex: c(20,100,15000))
# thresh_list = a numerical list denoting the target threshold value from 0 to 1
# cuttree_values = the output of the cl_lengthCut()  function which gives a table
#                 of GeneIDs assigned to clusters
# counts = normalized gene counts from RNA seq
# GO_annot = the Gene Onltology matrix for all GeneIDs

optimise_impute <- function(cl_list, thresh_list, cuttree_values, counts, GO_annot){
  
  scores <- list()
  
  for (i in cl_list){
    print(i)
    for (j in thresh_list){
      print(j)
      cl_tot <- as.character(i)
      cl <- cuttree_values[[cl_tot]][["GeneID_assignments"]]
      
      corr_clAll <- corr_per_clust(counts, cl, i)
      # For faster runtime, set cluster 2 to 0
      #corr_clAll[2] <- 0
      
      GOterms_perCl <- GO_per_cl(GO_annot, cl, i)
      GO_list_perCl <- GO_per_cl_list(GOterms_perCl, i)
      
      sc <- cross_val(n=10, GO_annot=GO_annot, clusters=cl,
                      GOterms_perCl=GOterms_perCl,
                      GO_list_perCl=GO_list_perCl,
                      corr_clAll=corr_clAll,
                      clust_total=i, thresh=j)
      
      scores[[paste0(i, "_", j)]] <- sc
    }
  }
  return(scores)
}


# NEEDS EDITING
#This function outputs a nested list of mean performance values 
# for a series of thresholds determined by a set interval
# Input: 
# interval - int; the interval from 0 to 1 used to generate a sequence
# of values that will serve as a threshold

cross_val_thresh <- function(interval){
  
  kfold_gobal <- list()
  
  thresh <- seq(from = 0, to = 1, by = interval)
  
  for (i in thresh){
    
    kfold <- cross_val(n=10, GO_annot=GO_table, clusters=cluster20,
                       GO_list_perCl=GOterms_perCl, 
                       corr_clAll=corr_allClusters,
                       wGO_allClusters=wGO_allClusters,
                       clust_total=20, thresh= i)
    
    kfold_df <- as.data.frame(do.call(rbind, kfold))
    
    Sensitivity <- mean(unlist(kfold_df$Sensitivity))
    Specificity <- mean(unlist(kfold_df$Specificity))
    Precision <- mean(unlist(kfold_df$Precision))
    Accuracy <- mean(unlist(kfold_df$Accuracy))
    F1_Score <- mean(unlist(kfold_df$`F1 Score`))
    
    
    kfold_gobal[[paste0(i)]][["DF"]] <- kfold_df
    kfold_gobal[[paste0(i)]][["Sensitivity"]] <- Sensitivity
    kfold_gobal[[paste0(i)]][["Specificity"]] <- Specificity
    kfold_gobal[[paste0(i)]][["Precision"]] <- Precision
    kfold_gobal[[paste0(i)]][["Accuracy"]] <- Accuracy
    kfold_gobal[[paste0(i)]][["F1_Score"]] <- F1_Score
  }
  return(kfold_gobal)
}

