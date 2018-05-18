# Taken from as we need a function that is not exported: https://github.com/bnikolic/RNested/blob/master/R/RNested.R
# Bojan Nikolic <bojan@bnikolic.co.uk>
# Nested sampling, implemented in R

#author bnikolic
sset.box <- function(box, nss,
                     llfn)
  {
    p <- apply(box,
               2,
               function(r) {
                 runif(nss, min=r[[1]], max=r[[2]])
               })
    # Also calculate the likelihood and set prior probability for all
    # the points in the starting set
    ll <- apply(p,
                1,
                llfn)
    lpr <- rep(0.0,
               nss)
    data.frame(p=I(p), ll=ll, lpr=lpr)
  }
#author bnikolic
mkPriorSamplePred <- function(s)
{
  function(ll, lp)
    {
      if(ll < s$ll)
        {
          return (FALSE)
        }
      else
        {
          if(lp >= s$lpr)
            {
              return (TRUE);
            }
          else
            {
              if ( exp(lp-s$lpr) > runif(1))
                {
                  return (TRUE);
                }
              else
                {
                  return (FALSE);
                }
            }
        }
    }
}
# bnikolic
rectOffseter <- function(scale)
  {
    function()
      {
        rnorm(length(scale), sd=scale)
      }
  }

#author bnikolic
randomEl <- function(cs)
  {
    N <- dim(cs)[1]
    return (cs[as.integer(runif(1, min=1, max=N)),]    )
  }

mkFixedRectProp <- function(scales)
  {
    off  <- rectOffseter(scales)
    function (x)
      {
        x+off()
      }
  }
#author bnikolic
CPChain <- function(s,
                    proposer,
                    n,
                    llf, lpf,
                    cs)
  {
    pred <- mkPriorSamplePred(s)
    pcurr <- randomEl(cs)$p
    ll <- 0
    lp <- 0
    r <- sapply(1:n, function(x) {
      pnew <- proposer(pcurr)
      llnew <- llf(pnew)
      lpnew <- lpf(pnew)
      if (pred(llnew, lpnew))
        {
          pcurr <<- pnew
          ll <<- llnew
          lp <<- lpnew
          return (TRUE);
        }
      else
        {
          return (FALSE);
        }
    })
    if( any(r) )
      {
        return (list(p=pcurr, ll=ll, lpr=lp));
      }
    else
      {
      return (FALSE);
    }
  }

mkSimplestPSampler <- function(s)
  {
    proposer <- mkFixedRectProp(c(s, s));
    function(worst,llf,lpf,cs) { CPChain(worst,
                                         proposer,
                                         100,
                                         llf,
                                         lpf,
                                         cs)
                               }
  }
#author bnikolic
mkCovarianceSampler <- function(s=1.0)
  {
    cvm <- 0
    proposer <- function(x)
      {
        N <- dim(cvm)[[1]]
        x+mvrnorm(n=1,
                  mu=rep(0,N),
                  Sigma=cvm*s)
      }
    function(worst,llf,lpf,cs) {
      cvm <<- cov(cs$p)
      CPChain(worst,
              proposer,
              100,
              llf,
              lpf,
              cs)
    }
  }

#author bnikolic
boxp <- function(box)
  {
    ff <- function (x)
      {
        for (i in 1:length(x))
            if (x[[i]] < box[[2*i-1]] || x[[i]] > box[[2*i]] )
              {
                return( -998)
              }
        return (0);
      }
    return(ff)
  }
#author bnikolic
nested.step <- function(cs,
                        llf, lpf,
                        psampler)
  {
    worsti <- which.min(cs$ll)
    worst <- cs[worsti,]
    newp <- psampler(worst, llf, lpf, cs)
    if (identical(newp,FALSE))
      {
        return (newp);
      }
    else
      {
      cs[worsti,] <- newp
      return (list(cs, worst));
      }
  }
#author bnikolic
nested.sample <- function(cs,
                          llf, lpf,
                          psampler,
                          cout=rbind(),
                          N=1)
  {
    for (i in 1:N)
      {
        r <- nested.step(cs, llf, lpf, psampler)
        if (identical(r,FALSE))
          break;
        cs <- r[[1]]
        nsamples <- dim(cout)[[1]]
        if (is.null(nsamples))
          nsamples <- 0;
        nlive <-    dim(cs)[[1]]
        weight <- exp(-1.0*nsamples/nlive) - exp(-1.0*(nsamples+1)/nlive)
        cout <- rbind(cout, data.frame(c(r[[2]][1,], w=weight)))
      }
    return (list(cs=cs, cout=data.frame(cout,
                          row.names=NULL)));
  }


