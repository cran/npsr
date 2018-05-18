##' @title nps.valid
##' @param Q Histogram of dataset (l*m*n vector)
##' @param l |Z|
##' @param m |X|
##' @param n |Y|
##' @param N Number of Repetitions for Nested Sampling
##' @param S Number of Starting Points for Nested Sampling
##' @description Calculates M_Valid
nps.valid = function(Q, l, m, n, N = sum(Q), S = sum(Q)){
  Rx = expand.grid(rep(list(1:m),l)) # m^l x l matrix
  Ry = expand.grid(rep(list(1:n),m )) #n^m x m matrix
  Rxy = cbind(Rx[rep(1:nrow(Rx), nrow(Ry)),],
              Ry[rep(1:nrow(Ry), each = nrow(Rx)),]) # (n^m)*(m^l) x l+m matrix
  # Create probability matrix for P(X = xi, Y = yi | zi)
  P = Create_P(l,m,n, Rxy)

  # Log-likelihood function for estimating the integral of theta_Z
  llf_z = function(t){
    return(log(Z_product(rep(unlist(t),m*n),Q)))
  }
  # Log-likelihood function for estimating the integral of theta_Z
  llf_xy = function(t){
    return(log(XY_product(unlist(t),P,Q)))
  }

  Int_z = estimate_integral(N=N,S=S,d=l,llf=llf_z, sample_theta = sample_theta)
  Int_xy = estimate_integral(N=N,S=S,d=(m^l)*(n^m),llf=llf_xy, sample_theta = sample_theta)
  return (Int_z*Int_xy)
}

#theta_z should be repeated to be l*m*n dimensions
Z_product = function(theta_z, Q){
  powered = unlist(theta_z)^Q
  product = prod(powered)
  return (product)
}

sample_theta = function(d){
  t = runif(d)
  t = t/sum(t) #make sure that all probablities add up one
  return (list(t))
}


