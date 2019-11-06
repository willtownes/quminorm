# functions for fitting and quantile normalizing to poisson-lognormal

#' @importFrom stats dpois
dpoilog<-function(x,mu,sig,log=FALSE,quadpts=1000){
  #Compute PMF of Poisson-lognormal
  #same functionality as poilog::dpoilog
  #poilog::dpoilog is faster, but this func should work better for
  #tiny probabilities (better numerical stability)
  #mu,sig are same as in poilog::dpoilog
  #if quadpts is a number, the quadrature points are recomputed based on
  #current parameter values
  #if quadpts is a list, we assume it is the result of a call to
  #quadpts<-statmod::gauss.quad.prob(m, dist="normal", mu=0, sigma=1)
  #where m is the number of quadrature points
  #ie that the quadpts are pre-computed based on the Normal(0,1) distribution
  if(is.list(quadpts)){
    q<-quadpts
    lam<-exp(mu+sig*q$nodes)
  } else {
    q<-statmod::gauss.quad.prob(quadpts,dist="normal",mu=mu,sigma=sig)
    lam<-exp(q$nodes)
  }
  w<-log(q$weights)
  f<-function(xi){dpois(xi,lam,log=TRUE)}
  #approximate the log of the integral
  g<-function(xi){matrixStats::logSumExp(w+f(xi))}
  #log PMF for all data points x and all quadrature points t
  #matrix with nrow=quadpts, ncol=length(x)
  lpmf<-vapply(x,g,FUN.VALUE=0.0)
  if(log==TRUE){ return(lpmf) } else { return(exp(lpmf)) }
}

#' @title Log-log curve for Poisson-lognormal PMF
#' @description Draw the probability mass function of a Poisson-lognormal
#'   distribution with both horizontal and vertical axes log transformed.
#' @name llcurve_poilog
#'
#' @param xmax the maximum quantile at which to evaluate the PMF.
#' @param lpar numeric vector containing the parameters mu,sig
#'   see \code{\link[sads]{dpoilog}} for details.
#' @param add logical, should the curve be added to the existing plot window?
#' @param use_sads logical, should the sads package be used to compute
#'   \code{\link[sads]{dpoilog}}? If not, an internal custom implementation is used.
#'   The custom implementation is less numerically stable than sads.
#' @param quadpts positive integer number of quadrature points, increase for
#'   better accuracy
#' @param ... additional arguments passed to \code{\link[graphics]{curve}}
#'
#' @return a list with components x and y of the points that were drawn is
#'   returned invisibly.
#'
#' @export
llcurve_poilog<-function(xmax,lpar,add=TRUE,use_sads=TRUE,quadpts=1000,...){
  #Draw the PMF curve on log-log axes for a Poisson-lognormal distribution
  #Curve goes from zero to xmax
  #lpar are the mu,sigma parameters (log scale, default for sads and poilog packages)
  if(use_sads){
    f<-function(t){
      sads::dpoilog(floor(expm1(t)),mu=lpar[1],sig=lpar[2],log=TRUE)
    }
  } else {
    f<-function(t){
      dpoilog(floor(expm1(t)),mu=lpar[1],sig=lpar[2],log=TRUE,quadpts=quadpts)
    }
  }
  graphics::curve(f,from=0,to=log1p(xmax),add=add,...)
}

#' @title Log-log curve for negative binomial PMF
#' @description Draw the probability mass function of a negative binomial
#'   distribution with both horizontal and vertical axes log transformed.
#' @name llcurve_nb
#'
#' @param xmax the maximum quantile at which to evaluate the PMF.
#' @param lpar numeric vector containing the parameters size,mu
#'   see \code{\link[stats]{dnbinom}} for details.
#' @param add logical, should the curve be added to the existing plot window?
#' @param ... additional arguments passed to \code{\link[graphics]{curve}}
#'
#' @return a list with components x and y of the points that were drawn is
#'   returned invisibly.
#'
#' @export
llcurve_nb<-function(xmax,lpar,add=TRUE,...){
  #Draw the PMF curve on log-log axes for a negative binomial distribution
  #Curve goes from zero to xmax
  #lpar are the size and mu parameters
  f<-function(t){
    stats::dnbinom(floor(expm1(t)),size=lpar[1],mu=lpar[2],log=TRUE)
  }
  graphics::curve(f,from=0,to=log1p(xmax),add=add,...)
}

#' @title Poisson-lognormal MLEs
#' @description Compute the maximum likelihood estimates for parameters of the
#'   Poisson-lognormal distribution.
#' @name poilog_mle
#'
#' @param x vector of non-negative integers (the data).
#' @param om optimization method to be used by \code{\link[sads]{fitpoilog}}.
#' @param ... additional arguments passed to \code{\link[sads]{fitpoilog}}.
#'
#' @return a named numeric vector containing the MLEs of mu,sigma (see
#'   \code{\link[sads]{dpoilog}} for details).
#'   The maximized value of the log-likelihood is attached as an attribute.
#'
#' @export
poilog_mle<-function(x,om="BFGS",...){
  #fit<-poilog::poilogMLE(x,startVals=st,zTrunc=FALSE,method=om,...)
  #mle<-fit$par
  #sads more stable than poilog
  fit<-sads::fitpoilog(x,trunc=NULL,method=om,skip.hessian=TRUE,...)
  mle<-stats::coef(fit) #fit$par
  attr(mle,"loglik")<-as.numeric(stats::logLik(fit)) #fit$logLval
  mle
}

#' @title Negative binomial MLEs
#' @description Compute the maximum likelihood estimates for parameters of the
#'   negative binomial distribution.
#' @name nb_mle
#'
#' @param x vector of non-negative integers (the data).
#' @param ... additional arguments passed to \code{\link[fitdistrplus]{fitdist}}.
#'
#' @return a named numeric vector containing the MLEs of size,mu (see
#'   \code{\link[stats]{dnbinom}} for details).
#'   The maximized value of the log-likelihood is attached as an attribute.
#'
#' @export
nb_mle<-function(x,...){
  fit<-fitdistrplus::fitdist(x,"nbinom",keepdata=FALSE,...)
  mle<-stats::coef(fit)
  attr(mle,"loglik")<-stats::logLik(fit)
  mle
}

mle_matrix<-function(m,lik=c("poilog","nb"),...){
  #m a matrix with samples in the columns
  #returns a data frame with nrows=ncols(m)
  #result includes MLEs for each model, log likelihood, and BIC
  lik<-match.arg(lik)
  mle_funcs<-list(poilog=poilog_mle,nb=nb_mle)
  f<-mle_funcs[[lik]]
  mle_func<-function(x){
    tryCatch({
      mle<-f(x,...)
      c(mle,loglik=attr(mle,"loglik"))
    },
    error=function(e){
      rep(NA,3)
    })
  }
  if(methods::is(m,"sparseMatrix")){
    apply_func<-function(m){
      m<-slam::as.simple_triplet_matrix(m)
      res<-slam::colapply_simple_triplet_matrix(m,mle_func)
      as.data.frame(do.call(rbind,res))
    }
  } else {
    apply_func<-function(m){ as.data.frame(t(apply(m,2,mle_func))) }
  }
  res<-apply_func(m)
  res$bic<- -2*res$loglik + log(nrow(m))*2
  res
}

#' @title Poisson-lognormal MLEs for columns of a matrix
#' @description Compute the maximum likelihood estimates for parameters of the
#'   Poisson-lognormal distribution for each column of a matrix.
#' @name poilog_mle_matrix
#'
#' @param m a matrix or sparse Matrix of non-negative integers (the data).
#' @param ... additional arguments passed to \code{\link{poilog_mle}}.
#'
#' @return a data frame whose rows correspond to the columns of m.
#'   The columns contain the MLEs for the parameters, the log-likelihood, and
#'   the \code{\link[stats]{BIC}} values.
#'
#' @export
poilog_mle_matrix<-function(m,...){
  #mle_matrix for poisson-lognormal
  #result includes mu,sigma params
  mle_matrix(m,"poilog",...)
}

#' @title Negative binomial MLEs for columns of a matrix
#' @description Compute the maximum likelihood estimates for parameters of the
#'   negative binomial distribution for each column of a matrix.
#' @name nb_mle_matrix
#'
#' @param m a matrix or sparse Matrix of non-negative integers (the data).
#' @param ... additional arguments passed to \code{\link{nb_mle}}.
#'
#' @return a data frame whose rows correspond to the columns of m.
#'   The columns contain the MLEs for the parameters, the log-likelihood, and
#'   the \code{\link[stats]{BIC}} values.
#'
#' @export
nb_mle_matrix<-function(m,...){
  mle_matrix(m,"nb",...)
}

#' @title Estimate negative binomial mean parameter from zero fraction
#' @description Given a fixed size parameter, use the fraction of zeros to
#'   estimate the mean parameter of a negative binomial
#'   distribution via the method of moments.
#' @name nb_pzero2mu
#'
#' @param lpz a numeric vector representing the natural log of the empirical
#'   fraction of zeros from count data.
#' @param size either a positive scalar or a vector of same length as lpz
#'   specifying the fixed size parameter (see \code{\link[stats]{dnbinom}}).
#'
#' @return a numeric vector of estimated scale parameters for each element of
#'   lpz.
#'
nb_pzero2mu<-function(lpz,size){
  #Assuming the data follows a negative binomial distribution
  #if we fix the size parameter to a specified value
  #we can infer the mean (mu) parameter from the fraction of zeros
  #lpz is a vector of the log of zero fraction for each cell
  #size can be a scalar or vector
  size*expm1(-lpz/size)
}

#' @title Estimate Poisson-lognormal scale parameter from zero fraction
#' @description Given a fixed sigma parameter, use the fraction of zeros to
#'   estimate the mu parameter of a Poisson-lognormal
#'   distribution via the method of moments.
#' @name poilog_pzero2mu
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
poilog_pzero2mu<-function(lpz,sig=2.5,lims=c(-100,100)){
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
      x-sads::dpoilog(0,mu,sig=s,log=TRUE)
    }
    stats::uniroot(f,lims)$root
  }
  if(length(sig)==1){ #single sig parameter for all cells
    return(vapply(lpz,inner,FUN.VALUE=1.0,s=sig))
  } else { #different tail parameter for each cell
    return(mapply(inner,lpz,sig))
  }
}
