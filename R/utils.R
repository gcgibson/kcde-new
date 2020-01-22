#' @rdname score_dist
#' @export score_dist
score_dist <- function(dist,truth){
  return (max(-10,log(sum(round(dist,1) == round(truth,1))/length(dist))))
}


#' @rdname do_box_cox
#' @export do_box_cox
do_box_cox <- function(y) {
  bc_params <- data.frame(gamma=.05)
  bc_params$lambda <- car::powerTransform(y + .05, family = "bcPower")$lambda
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
