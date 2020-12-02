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
# corr_cl = correlations per cluster
# clust_total = total number of clusters

impute <- function (cl_GOall, corr_clAll, clust_total, thresh){
  wGO_list <- list()
  
  for (i in 1:clust_total){
    
    GO_list <- rownames(cl_GOall[[i]])
    GO_list <- GO_list[order(GO_list)]
    
    corr_cl <- corr_clAll[[i]]
    corr_cl <- corr_cl[order(rownames(corr_cl)),]
    corr_cl <- corr_cl[,order(colnames(corr_cl))]
    gene_list <- rownames(corr_cl)
    
    cl_go <- cl_GOall[[i]]
    # should add gene that are not included in the GO data otherwise merge will yield weird results
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

