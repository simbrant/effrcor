rcorrmatrix <- function(d, alphad=1, tol=10**(-10), maxiter_cg=10**4) {
  
  d<-as.integer(d)
  if(d<=0 || !is.integer(d)) {
    stop("The dimension 'd' should be a positive integer!\n")
  }
  if(alphad<=0) {
    stop("'alphad' should be positive!\n")
  }
  
  # handling of d=1 and d=2
  if(d==1){
    rr <- matrix(1,1,1); return(rr)
  }
  else if(d==2) {
    rho <- runif(1,-1,1)
    rr <- matrix(c(1,rho,rho,1),2,2); return(rr) 
  } else {
    
    # call C++ routine if d >= 3
    rr <- .rcorrmatrix_cpp_(alphad, d)
    
    return(rr)
  }
}
