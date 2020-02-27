score_dist <- function(dist,truth){
  return (max(log(sum(round(dist,1) ==  round(truth,1))/length(dist)),-10))
}

predownload_data <- function(){
  data_obj <- list()
  data_obj_idx <- 1
  for (region in c('al', 'ak', 'az', 'ca', 'co', 'ct','de','fl','ga', 'hi', 'id', 'il','in', 'ia', 'ks', 'ky', 'la', 'me', 'md', 'ma', 'mi', 'mn', 'ms', 'mo', 'mt', 'ne', 'nv', 'nh', 'nj', 'nm', 'ny_minus_jfk', 'nc', 'nd', 'oh', 'ok', 'or','pa', 'ri', 'sc', 'sd', 'tn', 'tx', 'ut', 'vt', 'va', 'wa', 'wv', 'wi', 'wy','as','mp', 'dc', 'gu', 'pr', 'vi', 'ord', 'lax', 'jfk')){
    req <- curl_fetch_memory(paste0("https://delphi.midas.cs.cmu.edu/epidata/api.php?source=fluview&regions=",region,"&epiweeks=201240-201920"))
    req_json <- jsonlite::prettify(rawToChar(req$content))
    data_obj[[data_obj_idx]] <- jsonlite::fromJSON(req_json)$epidata
    data_obj_idx <- data_obj_idx +1
  }
  total_data <- do.call(rbind, data_obj)
  saveRDS(total_data,"data.RDS")
}


get_ili_for_epiweek_and_issue <- function(start_epiweek,end_epiweek,region){
  total_data <- readRDS("data.RDS")
  return (total_data[total_data$region == region & total_data$epiweek <= end_epiweek & total_data$epiweek >= start_epiweek ,])
  
}


plot_predictive_dist <- function(pred_dist,truth,test_idx,method,previous_point,location){
  p <-ggplot(data=data.frame(x = pred_dist),aes(x=x)) + geom_histogram()+geom_vline(xintercept =truth ) +xlim(0,10) + geom_vline(xintercept = previous_point,col='blue')
  ggsave(filename = paste0(location,"-",method,"-",test_idx),plot = p,device = "png")
}

plot_data_till_time_t <- function(data){
  p <-ggplot(data=data.frame(x = 1:length(data),y=data),aes(x=x,y=y)) + geom_line()
  ggsave(filename = paste0("data-",test_idx),plot = p,device = "png")
}

library(quantmod)
library(pander)
library(cdcfluview)
library(sarimaTD)
library(cdcfluutils)
library(ggplot2)

## declare region to be analyzed 
region_str <- c("Alabama","Alaska","Arizona","California","Colorado","Connecticut", "Delaware","Georgia","Hawaii","Idaho", "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana", "Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana", "Nebraska","Nevada","New Hampshire","New Jersey","New Mexico","New York","North Carolina","North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania","Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont", "Virginia","Washington","West Virginia","Wisconsin","Wyoming","As","mp","District of Columbia", "gu","Puerto Rico","Virgin Islands","ord","lax","New York City")

region_abbrv <- c('al', 'ak', 'az', 'ca', 'co', 'ct','de','ga', 'hi', 'id', 'il','in', 'ia', 'ks', 'ky', 'la', 'me', 'md', 'ma', 'mi', 'mn', 'ms', 'mo', 'mt', 'ne', 'nv', 'nh', 'nj', 'nm', 'ny_minus_jfk', 'nc', 'nd', 'oh', 'ok', 'or','pa', 'ri', 'sc', 'sd', 'tn', 'tx', 'ut', 'vt', 'va', 'wa', 'wv', 'wi', 'wy','as','mp', 'dc', 'gu', 'pr', 'vi', 'ord', 'lax', 'jfk')


kcde_score_df <- matrix(NA,nrow=55*68,ncol=3)
sarima_score_df <- matrix(NA,nrow=55*68,ncol=3)

score_idx <-1
for (location in 1:length(region_str)){
  region_str_local <- region_str[location]
  region_abbrv_local <- region_abbrv[location]
  train_data <- get_ili_for_epiweek_and_issue("201640","201840",region_abbrv_local)
  
  test_data <- get_ili_for_epiweek_and_issue("201840","201920",region_abbrv_local)
  
  
  ## GET SARIMA FIT
  
  
  
  sarima_fit_bc_transform <- fit_sarima(
    y =train_data$wili,
    ts_frequency = 52,
    transformation = "box-cox",
    seasonal_difference = F)
  
  print (sarima_fit_bc_transform)
  
  
  test_epiweeks <- c(201740:201752,paste0(2018,"0",1:9),201810:201818)
  
  
  for (test_idx  in 1:(length(test_epiweeks)-h)){
    h <- length(test_epiweeks)-h - test_idx
    epiweek_cutoff <- test_epiweeks[test_idx]
    truth <- tail(get_ili_for_epiweek_and_issue(test_epiweeks[test_idx + h ],test_epiweeks[test_idx + h],region_abbrv_local)$wili,1)
    
    
    
    ### GET SARIMA PREDS
    sarima_pred_dist<- simulate(
      object = sarima_fit_bc_transform,
      nsim = 10000,
      seed = 1,
      newdata = c(train_data$wili,test_data$wili[1:test_idx]),
      h = h
    )

    sarima_pred_dist <- sarima_pred_dist[!is.na(sarima_pred_dist)]
    sarima_pred_dist <- sarima_pred_dist[!is.nan(sarima_pred_dist)]
  
    kcde_pred_dist <- fit_and_predict_jags(data = c(train_data$wili,test_data$wili[1:test_idx]),epiweeks=c(substr(train_data$epiweek,5,7),substr(test_data$epiweek[1:test_idx],5,7)),
                                           params_ar=sarima_fit_bc_transform$arma)
    
  
    create_submission_file(p_samples = kcde_pred_dist,season = 1,epiweek = test_idx,regions = region_abbrv_local,method = "kcde")
    create_submission_file(p_samples = sarima_pred_dist,season = 1,epiweek = test_idx,regions = region_abbrv_local,method = "sarima")
    
  }
}


#### plots



library(ggplot2)
library(dplyr)
data_for_plot <- total_score_df[total_score_df$epiweek %in% 4:30,] %>% dplyr::group_by(region) %>% dplyr::summarize(kcde_ls=mean(kcde_ls),sarima_ls=mean(sarima_ls))

p <- ggplot(data_for_plot,aes(x=region,y=sarima_ls,col='sarimaTD')) + geom_point() +
  geom_point(aes(x=region,y=kcde_ls,col='sarimaTDS')) + ylim(-10,0)
ggsave(filename = "regional_results",plot = p,device = "png")


create_submission_file <- function(p_samples,season,epiweek,regions,method){
  
  
  p_samples <- array(pmax(pmin(100,p_samples),0.0),dim=dim(p_samples))
  
  # load required libraries
  library(cdcfluutils)
  library(predx)
  
  # create variables to pass into predx
  analysis_time_season = season 
  analysis_time_season_week <- 0
  weeks_in_first_season_year <- get_num_MMWR_weeks_in_first_season_year(analysis_time_season)
  last_analysis_time_season_week = 41
  max_prediction_horizon <- 4
  first_analysis_time_season_week = 10
  
  # create empy list to hold all region submission files
  predx_list <- list()
  
  #create region counter
  region_idx <- 1
  
  #iterate through regions to create submission file per region
  for (reg in regions){
    print (reg)
    predx_list[[region_idx]] <- get_predx_forecasts_from_trajectory_samples(trajectory_samples = p_samples, 
                                                                            location = reg, targets = c( paste0(1:4, " wk ahead")), 
                                                                            season = analysis_time_season, analysis_time_season_week = analysis_time_season_week, 
                                                                            first_analysis_time_season_week = first_analysis_time_season_week, 
                                                                            last_analysis_time_season_week = last_analysis_time_season_week, 
                                                                            predx_types = c("Sample", "Bin", "Point"))
    
    region_idx <- region_idx +1
    
  }
  
  #concate to a single object
  pred_to_write <- dplyr::rbind_list(predx_list) 
  # create the submission df
  library(dplyr)
  submission_df <- predx_to_submission_df(pred_to_write, ew = substr(epiweek,5,7), year = substr(epiweek,1,4), team = paste0(model,"-","projected"))  
  
  if (as.numeric(substr(epiweek,5,7)) <= 20){
    current_season_identifier <- substr(season,6,10)
  } else{
    current_season_identifier <- substr(season,1,4)
  }
  write.csv(submission_df,row.names = F, file =paste0(method,"-","EW",substr(epiweek,5,7),"-",regions,"-",current_season_identifier,".csv"))
  
}

