library(recommenderlab)

load(file = "RestaurantRecommenderSystem/bin_ratings_clean.RData")

model <- Recommender(data = bin_ratings_clean, method = "IBCF", parameter = list(method = "jaccard"))

#install.packages("pacman")

#pacman::p_load(pacman, dplyr, GGally, ggplot2, ggthemes, ggivs, httr, lubridate, plotly, rio, rmarkdown, shiny, stringr, tidyr)
# https://gist.github.com/dasilvaa10/3e1cd989c2b16aee6ae39edc9823406e
get_recommendations<-function(restaurant, food, hybrid = FALSE, w = NULL){
  
  zero_init<-rep(0,ncol(bin_ratings_clean)) 
  
  food_inds<-list()
  
  for (i in 1:length(food)){
    
    food_inds[[i]]<-grep(food[i], colnames(bin_ratings_clean))
    
  }  
  
  food_inds <- sort(unlist(food_inds))
  
  rest_inds<-grep(restaurant, colnames(bin_ratings_clean)) 
  
  select_inds<-intersect(food_inds,rest_inds)
  
  ones_in<-rep(1,length(select_inds))
  
  zero_init[select_inds]<-ones_in
  
  newdat<-rbind(colnames(bin_ratings_clean),zero_init)
  
  colnames(newdat)<-colnames(bin_ratings_clean)
  
  newdat<-as.data.frame(newdat)
  
  newdat<-suppressWarnings(apply(newdat,2,function(x) as.numeric(as.character(x))))
  
  newdat<-newdat[-1, , drop =  FALSE]
  
  newdat<-as.matrix(newdat)
  
  newdat_ratings<-as(newdat,"binaryRatingMatrix")
  
  if (hybrid == FALSE) {
    
    model <- Recommender(data = bin_ratings_clean, method = "IBCF", parameter = list(method = "jaccard"))
    
  } else {
    
    model <- HybridRecommender(
      
      Recommender(data = bin_ratings_clean, method = "IBCF", parameter = list(method = "jaccard")),
      
      Recommender(data = bin_ratings_clean, method = "RANDOM"),
      
      weights = w
      
    )
    
  } 
  
  recommendations <- predict(model, newdat_ratings, n = 100)
  
  recs<-as(recommendations, "list")
  
  recs_df<-as.data.frame(do.call("rbind",strsplit(recs$`1`, "_")))
  
  colnames(recs_df)<-c("Item","Restaurant")
  
  recs_splt<-split.data.frame(recs_df, recs_df$Restaurant )
  
  return(lapply(recs_splt, head, 5))
  
}

res1<-get_recommendations(restaurant="murphys", food=c("mussels","salmon"))
str(res1)
res1[2]


