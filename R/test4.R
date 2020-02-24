#utility function to extract season

extract_season <- function(date){
  if (as.numeric(substr(date,5,7)) <= 20){
    return (as.numeric(substr(date,1,4))-1)
  } else{
    return (as.numeric(substr(date,1,4)))
  }
}
do_box_cox <- function(y){
  est_bc_params <- car::powerTransform(y + bc_gamma, family = "bcPower")
  est_bc_params <- list(
    lambda = est_bc_params$lambda,
    gamma = bc_gamma)
  ret_list <- list()
}


#utility function to move forward h weeks

move_k_week_ahead <- function(epiweek,k){
  
  if (k==1){
    if (substr(epiweek,5,7) == 52){
      new_week <- "01"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (as.numeric(substr(epiweek,5,7)) < 9){
      new_week <- paste0("0",as.numeric(substr(epiweek,5,7))+1)
      new_year <- substr(epiweek,1,4)
    }else {
      new_week <- as.numeric(substr(epiweek,5,7)) + 1
      new_year <- substr(epiweek,1,4)
    }
  } else if (k==2){
    if (substr(epiweek,5,7) == 51){
      new_week <- "01"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 52){
      new_week <- "02"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    }else if (as.numeric(substr(epiweek,5,7)) < 8){
      new_week <- paste0("0",as.numeric(substr(epiweek,5,7))+2)
      new_year <- substr(epiweek,1,4)
    }else {
      new_week <- as.numeric(substr(epiweek,5,7)) + 2
      new_year <- substr(epiweek,1,4)
    }
  } else if (k==3){
    if (substr(epiweek,5,7) == 50){
      new_week <- "01"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 51){
      new_week <- "02"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 52){
      new_week <- "03"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (as.numeric(substr(epiweek,5,7)) < 7){
      new_week <- paste0("0",as.numeric(substr(epiweek,5,7))+3)
      new_year <- substr(epiweek,1,4)
    }else {
      new_week <- as.numeric(substr(epiweek,5,7)) + 3
      new_year <- substr(epiweek,1,4)
    }
  } else if (k==4){
    if (substr(epiweek,5,7) == 49){
      new_week <- "01"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 50){
      new_week <- "02"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 51){
      new_week <- "03"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (substr(epiweek,5,7) == 52){
      new_week <- "04"
      new_year <- as.numeric(substr(epiweek,1,4)) +1
    } else if (as.numeric(substr(epiweek,5,7)) < 6){
      new_week <- paste0("0",as.numeric(substr(epiweek,5,7))+4)
      new_year <- substr(epiweek,1,4)
    }else {
      new_week <- as.numeric(substr(epiweek,5,7)) + 4
      new_year <- substr(epiweek,1,4)
    }
  }
  
  
  return (paste0(new_year,new_week))
}

fit_and_predict_decomp <- function(train_data,sarima_pred_dist){
  decomp <- stl(x=ts(train_data$wili,frequency = 52),s.window = 52,t.window = 5)
  plot(decomp)
  seasonal <- decomp$time.series[,1]
  trend <- decomp$time.series[,2]
  residuals <- decomp$time.series[,3]
  
  trend_arima <- auto.arima(trend)
  seasonal_arima <- auto.arima(seasonal)
  residuals_arima <- auto.arima(residuals)
  original_data_scale_model <- auto.arima(train_data$wili)
  

  train_data$week <- as.numeric(substr(train_data$epiweek,5,7))
  empirical_vars <- train_data %>% group_by(week) %>% summarize(var=var(wili))
  current_week <- as.numeric(substr(move_k_week_ahead(tail(train_data$epiweek,1),h),5,7))
  var_at_week <- empirical_vars[empirical_vars$week == current_week,]$var
  return (rnorm(10000,mean = median.default(sarima_pred_dist),sd=sqrt(var_at_week)))
}