library(quantmod)
library(pander)
library(cdcfluview)
library(sarimaTD)
library(cdcfluutils)
library(ggplot2)

## declare region to be analyzed 
region_str <- c("Alabama","Alaska","Arizona","California","Colorado","Connecticut", "Delaware","Georgia","Hawaii","Idaho", "Illinois","Indiana","Iowa","Kansas","Kentucky","Louisiana", "Maine","Maryland","Massachusetts","Michigan","Minnesota","Mississippi","Missouri","Montana", "Nebraska","Nevada","New Hampshire","New Jersey","New Mexico","New York","North Carolina","North Dakota","Ohio","Oklahoma","Oregon","Pennsylvania","Rhode Island","South Carolina","South Dakota","Tennessee","Texas","Utah","Vermont", "Virginia","Washington","West Virginia","Wisconsin","Wyoming")

region_abbrv <- c('al', 'ak', 'az', 'ca', 'co', 'ct','de','ga', 'hi', 'id', 'il','in', 'ia', 'ks', 'ky', 'la', 'me', 'md', 'ma', 'mi', 'mn', 'ms', 'mo', 'mt', 'ne', 'nv', 'nh', 'nj', 'nm', 'ny_minus_jfk', 'nc', 'nd', 'oh', 'ok', 'or','pa', 'ri', 'sc', 'sd', 'tn', 'tx', 'ut', 'vt', 'va', 'wa', 'wv', 'wi', 'wy')


### first lets collect all the training data
state_data <- matrix(nrow=length(region_str),ncol=105)

for (location in 1:length(region_str)){
  region_str_local <- region_str[location]
  region_abbrv_local <- region_abbrv[location]
  train_data <- get_ili_for_epiweek_and_issue("201640","201840",region_abbrv_local)
  state_data[location,] <- train_data$wili
  
}

### plot state data

ggplot(data=data.frame(y=c(t(state_data[,1:52])),x=rep(1:52,48),group=rep(1:48,each=52)),aes(x=x,y=y,col=as.factor(group)))+ geom_line() + facet_wrap(~group)  

kcde_score_df <- matrix(NA,nrow=55*68,ncol=3)
sarima_score_df <- matrix(NA,nrow=55*68,ncol=3)
h <- 1
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
    seasonal_difference = T)
  
  print (sarima_fit_bc_transform)
  
  
  test_epiweeks <- c(201740:201752,paste0(2018,"0",1:9),201810:201818)
  
  
  for (test_idx  in 1:(length(test_epiweeks)-h)){
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
    sarima_pred_dist <- sarima_pred_dist[,h]
    
    sarima_pred_dist <- sarima_pred_dist[!is.na(sarima_pred_dist)]
    sarima_pred_dist <- sarima_pred_dist[!is.nan(sarima_pred_dist)]
    
    #KCDE model
    decomp <- stl(ts(c(train_data$wili,test_data$wili[1:test_idx]),frequency = 52),s.window = 52,t.window = 4)
    
    kcde_pred_dist<- simulate(
      object = sarima_fit_bc_transform,
      nsim = 10000,
      seed = 1,
      newdata = decomp$time.series[,1] + decomp$time.series[,2],
      h = h
    )
    kcde_pred_dist <- kcde_pred_dist[,h]
    
    kcde_pred_dist <- kcde_pred_dist[!is.na(kcde_pred_dist)]
    kcde_pred_dist <- kcde_pred_dist[!is.nan(kcde_pred_dist)]
    
    
    sarima_score_df[score_idx,] <- c(region_abbrv_local,test_idx,score_dist( sarima_pred_dist,truth) )
    kcde_score_df[score_idx,] <- c(region_abbrv_local,test_idx,score_dist( kcde_pred_dist,truth) )
    plot_predictive_dist(pred_dist =sarima_pred_dist,truth = truth,test_idx = test_idx,method="sarima" ,previous_point = tail(test_data$wili,1),location=location)
    plot_predictive_dist(pred_dist =kcde_pred_dist,truth=truth,test_idx = test_idx,method="kcde",previous_point = tail(test_data$wili,1),location=location) 
    ### SCORE DATAFRAME COUNTER
    score_idx <- score_idx +1
    
    
    
  }
}

kcde_score_df <- data.frame(kcde_score_df)
colnames(kcde_score_df) <- c("region","epiweek","kcde_ls")
saveRDS(object = kcde_score_df,file="kcde_score_df_2016-2017_h_2")

sarima_score_df <- data.frame(sarima_score_df)
colnames(sarima_score_df) <- c("region","epiweek","sarima_ls")
saveRDS(object = sarima_score_df,file="sarima_score_df_2016-2017_h_2")


print (mean(as.numeric(as.character(kcde_score_df$kcde_ls)),na.rm=T))
print (mean(as.numeric(as.character(sarima_score_df$sarima_ls)),na.rm=T))

print (mean(as.numeric(as.character(kcde_score_df$kcde_ls)),na.rm=T))-print (mean(as.numeric(as.character(sarima_score_df$sarima_ls)),na.rm=T))
