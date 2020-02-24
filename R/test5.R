library(GAS)

fit_and_predict_decomp <- function(data,sarima_fit_bc_transform){
  data_for_sarima <- data$wili 
  decomp <- stl(ts(data_for_sarima,frequency = 52),t.window = 4,s.window = 52)
  reconstructed_series <- decomp$time.series[,1] + decomp$time.series[,2]
  
  sarima_pred_dist<- simulate(
      object = sarima_fit_bc_transform,
      nsim = 10000,
      seed = 1,
      newdata = c(reconstructed_series),
      h = h
    )
    
    
  return (sarima_pred_dist)
}