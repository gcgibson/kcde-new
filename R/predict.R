#####################'
#' Predict function
#' Inputs:
#' ts = time series 
#' k = number of neighbors to use
#' h = horizon 
#' num_lags = 4
#' sigma = bandwith
#' Outputs:
####################'
source('R/kernels.R')

#' @rdname predict_seasonal_kcde
#' @export predict_seasonal_kcde
predict_seasonal_kcde <- function(ts,k,h,num_lags,sigma,eta,sd){
  
  ret_list <- do_box_cox(ts)
  lambda <- ret_list[[2]]
  ts <- ret_list[[1]]
  
  
  # create matrix with dimensions length(ts) by number of lags 
  design_matrix <- matrix(NA,nrow=length(ts),ncol=num_lags)
  for (i in 0:(num_lags-1)){
    design_matrix[,(4-i)] <- (Lag(ts,i))
  }
  # subset to remove NA's created when lags
  design_matrix_complete <- design_matrix[complete.cases(design_matrix),]
  
  # get vector to compare to similarity matrix
  test_x <- tail(ts,num_lags)
  
  #covariate similarities
  sim_mat <- apply(design_matrix_complete,1,function(x){rbf(x,test_x,h = sigma)})
  similarities <- colSums(sim_mat)/sum(colSums(sim_mat))
  names(similarities) <- 1:length(similarities)
  similarities <- similarities[1:(length(similarities)-h)]
  
  #periodic similarities
  periodic_similarities<- unlist(lapply(1:length(ts),function(x){periodic(x,length(ts)+1,eta = eta)}))
  periodic_similarities <- periodic_similarities[num_lags:(length(periodic_similarities)-h)]
  
  #product kernel similarities
  product_similarities <- similarities*periodic_similarities
  
  # get the top k nearest neighbors by similarity
  top_k_similar <- tail(sort(product_similarities),k)
  
  # get indices of design matrix that map to most similar
  top_k_similar_ys <-as.numeric(names(top_k_similar))
  
  # get h-step ahead ts values for each k nearest neighbor of design matrix
  # with probabilities correspnding to the similarities
  if (sum(top_k_similar) == 0){
    top_k_similar <- rep(1/length(top_k_similar),length(top_k_similar))
    print (paste("no nearest neighbor for",sigma,eta))
  }
  pred_density <- pmax(rnorm(10000,mean=sample(ts[top_k_similar_ys+num_lags-1+h],10000,prob=top_k_similar,replace=T),sd=sd),0)
  pred_density <- invert_bc_transform(pred_density,lambda,.5)
  return(pred_density)
}



