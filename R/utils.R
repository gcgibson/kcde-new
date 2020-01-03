

do_box_cox <- function(y) {
  bc_params <- data.frame(gamma=.5)
  bc_params$lambda <- car::powerTransform(y + .1, family = "bcPower")$lambda
  transformed_y <- car::bcPower(
    U = y + bc_params$gamma,
    lambda = bc_params$lambda)
  
  ret_list <- list()
  ret_list[[1]] <- transformed_y
  ret_list[[2]] <- bc_params$lambda 
  return (ret_list)
}


invert_bc_transform <- function(b, lambda, gamma) {
  ## Two steps: 1) undo box-cox 2) subtract offset gamma
  
  ## 1) undo box-cox
  if(abs(lambda) <= 1e-10) {
    z <- exp(b)
  } else {
    z <- (lambda * b + 1)^(1 / lambda)
  }
  
  ## 2) Subtract gamma
  return(z - gamma)
}
