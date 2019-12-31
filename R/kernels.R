#####################'
#' Set of Kernel Functions
#' Currently implemented 
#' 1) Gaussian
#' 2) periodic
####################'

########'
#' Gussian Kernel 
#' Inputs:
#' x = variable 1
#' y = variable 2
#' h = inverse of kernel bandwith
#' Outputs:
#' similarity between x and y
########'

#' @rdname rbf
#' @export rbf
rbf <- function(x, y, h = 1)
{
  exp(- h * (x - y) ^ 2)
}

########'
#' Periodic Kernel 
#' Inputs:
#' x = variable 1
#' y = variable 2
#' h = inverse of kernel bandwith
#' rho = periodicity of kernel
#' eta = smoothing parameter
#' Outputs:
#' similarity between x and y
########'

#' @rdname periodic
#' @export periodic
periodic <- function(x, y,h=1, rho=pi/52,eta=1e10)
{
  exp(- (sin((rho^h)*(x-y))^2)/(2*(eta^h)))
}