# functions for fitting and quantile normalizing to poisson-lognormal

#' Poisson-lognormal MLEs
#'
#' Compute the maximum likelihood estimates for parameters of the
#' Poisson-lognormal distribution.
#'
#' @param x vector of non-negative integers (the data).
#' @param om optimization method to be used by \code{\link[sads]{fitpoilog}}.
#' @param ... additional arguments passed to \code{\link[sads]{fitpoilog}}.
#'
#' @return a named numeric vector containing the MLEs of mu,sigma (see
#'   \code{\link[sads]{dpoilog}} for details).
#'   The maximized value of the log-likelihood is attached as an attribute.
#'
#' @importFrom sads fitpoilog
#' @importFrom stats coef logLik
#' @export
poilog_mle<-function(x,om="BFGS",...){
  fit<-fitpoilog(x,trunc=NULL,method=om,skip.hessian=TRUE,...)
  mle<-coef(fit) #fit$par
  attr(mle,"loglik")<-as.numeric(logLik(fit)) #fit$logLval
  mle
}

#' Poisson-lognormal MLEs for columns of a matrix
#'
#' Compute the maximum likelihood estimates for parameters of the
#' Poisson-lognormal distribution for each column of a matrix.
#'
#' @param m a matrix or sparse Matrix of non-negative integers (the data).
#' @param ... additional arguments passed to \code{\link{poilog_mle}}.
#'
#' @return a data frame whose rows correspond to the columns of m.
#'   The columns contain the MLEs for the parameters, the log-likelihood, and
#'   the \code{\link[stats]{BIC}} values.
#'
#' @importFrom methods is
#' @export
poilog_mle_matrix<-function(m,...){
  #m a matrix with samples in the columns
  #returns a data frame with nrows=ncols(m)
  #result includes mu,sigma params, log likelihood, and BIC
  mle_func<-function(x){
    tryCatch({
      mle<-poilog_mle(x,...)
      c(mle,loglik=attr(mle,"loglik"))
    },
    error=function(e){
      rep(NA,3)
    })
  }
  if(is(m,"sparseMatrix")){
    res<-as.data.frame(do.call(rbind,colapply_full(m,mle_func)))
  } else {
    res<-as.data.frame(t(apply(m,2,mle_func)))
  }
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
