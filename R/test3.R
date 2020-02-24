#utility function to extract season

extract_season <- function(date){
  if (as.numeric(substr(date,5,7)) <= 20){
    return (as.numeric(substr(date,1,4))-1)
  } else{
    return (as.numeric(substr(date,1,4)))
  }
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

fit_and_predict_decomp <- function(train_data){
  decomp <- stl(x=ts(train_data$wili,frequency = 52),s.window = 20,t.window = 5)
  #plot(decomp)
  seasonal <- decomp$time.series[,1]
  trend <- decomp$time.series[,2]
  residuals <- decomp$time.series[,3]
  
  ## seasonal model
  seasonal_df <- data.frame(y = as.vector(seasonal), x= substr(train_data$epiweek,5,7),season=unlist(lapply(train_data$epiweek,extract_season)))

  library(ggplot2)
  #ggplot(seasonal_df,aes(x=x,y=y,col=as.factor(season))) + geom_point()
  
  ## trend model 
  trend_df <- data.frame(y = as.vector(trend),x= substr(train_data$epiweek,5,7),season=unlist(lapply(train_data$epiweek,extract_season)))
  #ggplot(trend_df,aes(x=x,y=y,col=as.factor(season))) + geom_point()
  
  
  
  
  ## h= 1 forecast
  
  model.loc = ("rw_intercept.txt")
  jagsscript = cat("
  model {  
     mu ~ dnorm(0, 0.01); 
     tau.pro ~ dgamma(0.001,0.001); 
     sd.pro <- 1/sqrt(tau.pro);
     
     # Week specific RW
     week_rw[1] <- 0;
     for (i in 2:53){
        week_rw[i] ~ dnorm(week_rw[i-1],100);
     }
  
     # week/season rw
     week_season_rw[1] <- 0;
     for (i in 2:N){
      week_season_rw[i] ~ dnorm(week_season_rw[i-1],100);
     }
     for(i in 1:N) {
        Y[i] ~ dnorm(week_season_rw[i] + week_rw[week_idx[i]], 10000);
     }
  }  
  ", 
                   file = model.loc)
  
  week_for_h_ahead <- as.numeric(substr(move_k_week_ahead(paste0("2018",tail(trend_df$x,1)),h),5,7))
  jags.data = list(Y = c(trend_df$y,NA),week_idx=c(as.numeric(trend_df$x),week_for_h_ahead),N=nrow(trend_df)+1)
  jags.params = c("week_season_rw","week_rw","Y")
  mod_rw_intercept_trend = jags(jags.data, parameters.to.save = jags.params, 
                          model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 1, 
                          n.iter = 10000, DIC = TRUE)
  
  
  
  
  
  model.loc = ("rw_intercept_season.txt")
  jagsscript = cat("
                   model {  
                   mu ~ dnorm(0, 0.01); 
                   tau.pro ~ dgamma(0.001,0.001); 
                   sd.pro <- 1/sqrt(tau.pro);
                   
                   # Week specific RW
                   week_rw[1] <- 0;
                   for (i in 2:53){
                    week_rw[i] ~ dnorm(week_rw[i-1],100);
                   }
                   
                   # week/season rw
                   week_season_rw[1] <- 0;
                   for (i in 2:N){
                   week_season_rw[i] ~ dnorm(week_season_rw[i-1],1000);
                   }
                   for(i in 1:N) {
                   Y[i] ~ dnorm(week_season_rw[i] + week_rw[week_idx[i]], 100000);
                   }
                   }  
                   ", 
                   file = model.loc)
  
  
  jags.data = list(Y = c(seasonal_df$y,NA),week_idx=c(as.numeric(seasonal_df$x),week_for_h_ahead),N=nrow(seasonal_df)+h)
  jags.params = c("week_season_rw","week_rw","Y")
  mod_rw_intercept_season = jags(jags.data, parameters.to.save = jags.params, 
                          model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 1, 
                          n.iter = 10000, DIC = TRUE)
  
  pred_dist <- mod_rw_intercept_season$BUGSoutput$sims.matrix[,paste0("Y[",nrow(seasonal_df)+h,"]")] + mod_rw_intercept_trend$BUGSoutput$sims.matrix[,paste0("Y[",nrow(seasonal_df)+h,"]")]
  
  return (pred_dist)
}