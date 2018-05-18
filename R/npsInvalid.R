##' @title M_Invalid
##' @param  Q List of unique observations, should be l\*m\*n length
##' @param l |Z|
##' @param m |X|
##' @param n |Y|
##' @param N Number of Repetitions for Nested Sampling
##' @param S Number of Starting Points for Nested Sampling
##' @description Calculates the ML_Invalid
nps.invalid = function(Q,l,m,n, N = sum(Q), S = sum(Q)){
  air = M_air(Q,l,m,n)
  excl = M_excl(Q,l,m,n, N, S)
  air_excl = M_air_excl(Q,l,m,n)
  max = max(c(air,excl,air_excl))

  return(list(air = air,excl = excl ,air_excl = air_excl, max = max))
}
##' @description Calculates the marginal likelihood of M_excl
##' @title M_excl
##' @param Q Histogram of dataset (l\*m\*n vector)
##' @param l |Z|
##' @param m |X|
##' @param n |Y|
##' @param N Number of Repetitions for Nested Sampling
##' @param S Number of Starting Points for Nested Sampling
##' @return The probability that the observations where created from a model which violates the exclusion criterion but not the as-if-randomness criterion
M_excl = function(Q,l,m,n,N = sum(Q), S = sum(Q)){
  # Create all Rxy which violates the exclusion criterion
  h = create_excl_vio_rxy(l,m,n)
  #valid Rxy
  Rx = expand.grid(rep(list(1:m),l)) # m^l x l matrix
  Ry = expand.grid(rep(list(1:n),m )) #n^m x m matrix
  Ry = matrix(unlist(rep(Ry,2)),n^m,l*m) # repeat columns to allow acces dependent on z

  ZXY = expand.grid((1:l), (1:m), (1:n)) # l*m*n x 3 matrix

  max_I = 0;
  sample_theta_xy = function(d){
    t = rep(0,d)
    while(t[d] == 0){ #ensures that model is invalid with the exclusion criterion
      t = runif(d)
      t = t/sum(t) #make sure that all probablities add up to one
    }
    return (list(t))
  }
  sample_theta_z = function(d){
    t = runif(d)
    t = t/sum(t) #make sure that all probablities add up to one
    return (list(t))
  }

  integers = apply(h, 1, function(hi){
    Ry_hi = rbind(Ry, t(data.matrix(hi, rownames.force = NA)));
    Rxy = cbind(Rx[rep(1:nrow(Rx), nrow(Ry_hi)),],
                Ry_hi[rep(1:nrow(Ry_hi), each = nrow(Rx)),]) # ((n^m)+1)*(m^l) x l+m*l matrix
    # Vector of all possible Observations (Qi and ZXYi have to relate to the same obervations)
    # we add the invalid response function (hi) for y to valid response functions
    # P is the probability matrix for the realization of an observation given
    # a pair of response function for x and y (rows are observation, columns are repsponse functions)
    P = Create_P(l,m,n,Rxy, y_zx_depenendent = TRUE)

    # Log-likelihood function for estimating the integral of theta_Z
    llf_z = function(t){
      return(log(Z_product(rep(unlist(t),m*n),Q)))
    }
    # Log-likelihood function for estimating the integral of theta_Z
    llf_xy = function(t){
      return(log(XY_product(unlist(t),P,Q))) #TODO fix for changed T
    }
    Int_z = estimate_integral(N=N,S=S,d=l,llf=llf_z, sample_theta_z)
    Int_xy = estimate_integral(N=N,S=S,d=(m^l)*((n^m)+1),llf=llf_xy,sample_theta_xy)
    return(Int_z * Int_xy)

  })
  return (max(integers))
}
create_excl_vio_rxy = function(l,m,n){
  #our column convention: var1 = (z0,x0), var2 = (z1,x0) ... var(l+1) = (z0,x1) and so on
  h = expand.grid(rep(list(1:n),l*m))
  #logical matrix for all zi, representing wether or not h is independent on zi
  h_indep_zi = sapply(c(1:l), function(i) {
    ind = rowSums(h[,i]==h[,seq(i,l*m,l)]) == l
    return (ind)
  })
  #logical vector which represents which hi is actually dependent on z
  h_dep_z = rowSums((h_indep_zi[,1]&h_indep_zi[,1:l])) != l
  h = h[h_dep_z,]
  return (h)
}

##' @description Calculates the marginal likelihood M_air
##' @title M_air
##' @param Q Histogram of dataset (l\*m\*n vector)
##' @param l |Z|
##' @param m |X|
##' @param n |Y|
##' @return The probability that the observations where created from
##' a model which violates the as-if-randomness criterion but not the exclusion criterion
M_air = function(Q, l,m,n){
  b = (l*(m^l)*((n)^m)) / length(Q)
  result = I_optimized(Q,b)
}
##' @description Calculates the marginal likelihood M_air_excel
##' @title m_air
##' @param Q Histogram of dataset (l\*m\*n vector)
##' @param l |Z|
##' @param m |X|
##' @param n |Y|
##' @return The probability that the observations where created from
##' a model which violates the as-if-randomness criterion but not the exclusion criterion
M_air_excl = function(Q, l,m,n){
  b = (l*(m^l)*(n^(m*l))) / length(Q)
  result = I_naive(Q,b)
}

I_naive = function(Q,b){
  a = length(q)
  b = as.bigz(b)
  numerator = prod(apply(as.matrix(Q+b),1, gamma))
  denom_gsum = factorial((sum(Q+as.bigz(4))))
  denom_gQ = gamma(as.bigz(b))^(a);
  result = numerator/(denom_gsum * denom_gQ);
  return (as.numeric(result));
}
I_optimized = function(Q,b){
  q_max = max(Q)
  qi_max = which(q_max == Q)[1]
  Q_num = rep(Q);
  Q_num[qi_max] = -b + 1; # => gamma(Qj+b) = 1, removing element from product

  numerator = prod(unlist(lapply((Q_num+b), gamma)))

  d_sum_upper = sum(Q+4)
  d_sum_lower = q_max+b
  d_gb = gamma(b)^length(Q)

  result = product_fraction(unlist(lapply((Q_num+b), gamma)),c(d_sum_lower:d_sum_upper, d_gb))
  return (result)
}
##' Reduces out factors of fraction of products and calculates the fraction
##' Analog to prod(num)/prod(den)
##' @param num vector of factors of the numerator
##' @param den vector of factors of the denominator
product_fraction = function(num, den){
  rnum = sort(num, decreasing = TRUE)
  rden = sort(den, decreasing = TRUE)
  num_i = 1 #index to iterate over the numerator
  den_i = 1 #index to iterate over the denominator
  fraction = 1 # fraction we will calculate
  while(num_i <= length(rnum) || den_i <= length(rden)){
    if(num_i > length(rnum) || fraction > 1 && den_i <= length(rden)){
      fraction = fraction/rden[den_i]
      den_i = den_i+1
    }
    else{
      fraction = fraction*rnum[num_i]
      num_i = num_i+1
    }
  }
  return(fraction)
}
