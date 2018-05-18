##' @title testIc
##' @description Tests the instrumental constraints on the given dataframe using entropy
##' @param df Dataframe with z, x and y
##' @return FALSE if the data violates the constraints otherwise TRUE
nps.necessary = function(df){
  Iyzx =- condinformation(df$y,df$z,df$x)
  Ixz = mutinformation(df$x,df$z)
  Hx = entropy(df$x)
  passed = (Iyzx + Ixz) <= Hx
  result = list(Iyzx = Iyzx, Ixz = Ixz, Hx = Hx, passed=passed)
  return(result)
}

