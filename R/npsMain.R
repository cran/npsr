##' Main function of the package.
##' @param df Dataframe with columns z,x and y
##' @param l Number of bins used to discretize Z
##' @param m Number of bins used to discretize X
##' @param n Number of bins used to discretize Y
##' @param N Number of Repetitions for Nested Sampling
##' @param S Number of Starting Points for Nested Sampling
##' @return result object of the test including the fields: nt, valid, invalid, ratio
##' @export
##' @import stats infotheo MASS gmp
##' @examples
##' nps.test(data.frame(x = runif(3), y = runif(3), z = runif(3)),2,2,2, 3, 3)
nps.test = function(df,l,m,n, N, S){
  # Discretize values with given dimensions
  dx = unname(unlist(discretize(df$x,nbins=m, disc="equalwidth")))
  dy = unname(unlist(discretize(df$y,nbins=m, disc="equalwidth")))
  dz = unname(unlist(discretize(df$z,nbins=m, disc="equalwidth")))

  ddf = data.frame(z = dz, x = dx, y = dy)

  nt = nps.necessary(ddf);
  if(!nt$passed) return(FALSE)

  ZXY = expand.grid((1:l), (1:m), (1:n))
  Q = apply(ZXY,1,function(row) {
     qi = sum(ddf$z == row[1] & ddf$x == row[2] & ddf$y == row[3])
     return (qi)
  })
  if(missing(N)){
    N = sum(Q)
  }
  if(missing(S)){
    S = sum(Q)
  }
  invalid = nps.invalid(Q,l,m,n)
  valid = list(valid = nps.valid(Q,l,m,n, N, S))
  ratio = list(ratio = valid$valid/invalid$max)

  result = merge(nt, merge(valid, merge(invalid, ratio)))
  return (result)
}
