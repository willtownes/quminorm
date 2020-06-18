# functions for fitting Poisson-lognormal distributions

#' @importFrom sads fitpoilog
#' @importFrom bbmle coef logLik
poilog_mle1<-function(x,...){
  #the '1' means it is only applied to one column instead of a matrix
  tryCatch({
    fit<-fitpoilog(x,trunc=NULL,method="BFGS",skip.hessian=TRUE,...)
    res<-coef(fit) #fit$par
    c(res,loglik=as.numeric(logLik(fit))) #fit$logLval
  },
  error=function(e){
    rep(NA,3)
  })
}

poilog_mle_dense<-function(m, mc.cores=1){
  #m a matrix with samples in the columns
  #returns a data frame with nrows=ncols(m)
  #result includes mu, sig params, log likelihood, and BIC
  if(mc.cores>1){ #parallel apply
    res<-colapply_dense_parallel(m, poilog_mle1, mc.cores=mc.cores)
    res<-simplify2array(res)
  } else { #serial apply
    res<-apply(m,2,poilog_mle1)
  }
  res<-as.data.frame(t(res))
  res$bic<- -2*res$loglik + log(nrow(m))*2
  res
}

poilog_mle1_nz<-function(xnz,n){
  #thin wrapper if only the nonzero values were provided (xnz)
  #n: the length of the data vector if zeros were included
  poilog_mle1(c(xnz,rep.int(0,n-length(xnz))))
}

poilog_mle_sparse<-function(m, mc.cores=1){
  if(mc.cores>1){
    res<-colapply_sparse_nonzero(m, poilog_mle1_nz, n=nrow(m),
                                 mc.cores=mc.cores)
  } else {
    res<-colapply_sparse_full(m,poilog_mle1)
  }
  res<-as.data.frame(do.call(rbind,res))
  res$bic<- -2*res$loglik + log(nrow(m))*2
  res
}

#' Estimate Poisson-lognormal scale parameter from zero fraction
#'
#' Given a fixed sigma parameter, use the fraction of zeros to
#' estimate the mu parameter of a Poisson-lognormal
#' distribution via the method of moments.
#'
#' @param lpz a numeric vector representing the natural log of the empirical
#'   fraction of zeros from count data.
#' @param sig either a positive scalar or a vector of same length as lpz
#'   specifying the fixed sigma parameter (see \code{\link[sads]{dpoilog}}).
#' @param lims interval of possible values for the scale
#'   parameter (mu), passed to \code{\link[stats]{uniroot}}.
#'
#' @return a numeric vector of estimated scale parameters (mu) for each element
#' of lpz.
#'
#' @importFrom sads dpoilog
#' @importFrom stats uniroot
poilog_pzero2mu<-function(lpz,sig=2.5,lims=c(-200,200)){
  #Assuming the data follows a Poisson-lognormal
  #if we fix the 'sig' parameter to a specified value
  #we can infer the mu parameter from the fraction of zeros
  #mu != mean of the lognormal, but e^mu is median of lognormal
  #mu can be negative.
  #lpz is a vector of the log of zero fraction for each cell
  #lims are the lower,upper bounds for the mu parameter passed to uniroot
  inner<-function(x,s){
    #x is an element of lpz
    if(is.na(s)){ return(NA) }
    f<-function(mu){
      x-dpoilog(0,mu,sig=s,log=TRUE)
    }
    uniroot(f,lims)$root
  }
  if(length(sig)==1){ #single sig parameter for all cells
    return(vapply(lpz,inner,FUN.VALUE=1.0,s=sig))
  } else { #different tail parameter for each cell
    return(mapply(inner,lpz,sig))
  }
}
