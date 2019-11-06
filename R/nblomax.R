# Find MLE of Negative-Binomial Lomax family
# x|mu~NegBinom(shape=phi, mean=mu)
# mu~Lomax(tail=alpha, scale=lambda)
# if alpha known, can use quadrature via equivalent
# z~Beta(1,alpha)
# mu = lambda*z/(1-z)
# if alpha unknown, can use quadrature via equivalent
# z~Exponential(1)
# mu = lambda*(exp(z/alpha)-1)
# Lomax prior is conjugate to neg binom only if lambda=phi, see bnbinom.R

#' @title Negative binomial-Lomax PMF
#' @description Compute probability mass function of a negative binomial
#'   distribution with mean parameter distributed as Lomax.
#' @name dnblomax
#'
#' @param x a vector of (non-negative integer) quantiles.
#' @param tail a positive scalar controlling the heaviness of the Lomax tail.
#' @param scale a positive scalar.
#' @param shape a positive scalar negative binomial parameter; can be Inf.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param quadpts either a single integer determining the number of quadrature
#'   points, or a list resulting from
#'   statmod::gauss.quad.prob(m, dist="uniform", l=0, u=1)
#'   where m is the number of quadrature points.
#'
#' @details This function computes the likelihood the data x was generated
#'   from the model mu~Lomax(tail,scale), x~negative binomial(shape,mu).
#'   If the shape parameter is Inf this is a Poisson-Lomax model. If shape is 1,
#'   it is a geometric-Lomax model. The Lomax mixing distribution induces
#'   a heavy power-law tail with exponent (1+tail). The smaller the tail
#'   parameter, the heavier the tail of the distribution. If tail<=1 there is no
#'   mean, if tail<2 there is no variance. A larger scale parameter means fewer
#'   zeros and more large values. Shape and scale have no effect on the
#'   power-law tail. If quadpts is a number, the quadrature points are
#'   recomputed based on current parameter values. If quadpts is a list, we
#'   assume it is the result of a call to
#'   quadpts<-statmod::gauss.quad.prob(m, dist="uniform", l=0, u=1)
#'   where m is the number of quadrature points, ie that the quadpts are
#'   pre-computed based on the Uniform(0,1) distribution.
#'
#' @return a vector of (log) probabilities
#' @import statmod,matrixStats
#' @export
dnblomax<-function(x, tail=1, scale=1, shape=1, log=FALSE, quadpts=1000){
  if(is.list(quadpts)){
    q<-quadpts
    mu<-scale*((1-q$nodes)^(-1/tail)-1) #inverse CDF of Lomax
  } else {
    q<-statmod::gauss.quad.prob(quadpts,dist="beta",alpha=1,beta=tail)
    mu<-scale*q$nodes/(1-q$nodes)
  }
  w<-log(q$weights)
  if(is.infinite(shape)){ #poisson
    f<-function(xi){dpois(xi,mu,log=TRUE)}
  } else { #negative binomial
    f<-function(xi){dnbinom(xi,size=shape,mu=mu,log=TRUE)}
  }
  #approximate the log of the integral
  g<-function(xi){matrixStats::logSumExp(w+f(xi))}
  #log PMF for all data points x and all quadrature points t
  #matrix with nrow=quadpts, ncol=length(x)
  lpmf<-vapply(x,g,FUN.VALUE=0.0)
  if(log==TRUE){ return(lpmf) } else { return(exp(lpmf)) }
}

dglomax<-function(x,...){ dnblomax(x,shape=1,...) }

dplomax<-function(x,...){ dnblomax(x,shape=Inf,...) }

#' @title Log-log curve for negative binomial-Lomax PMF
#' @description Draw the probability mass function of a negative binomial-Lomax
#'   distribution with both horizontal and vertical axes log transformed.
#' @name llcurve_lomax
#'
#' @param xmax the maximum quantile at which to evaluate the PMF.
#' @param lpar numeric vector containing the tail, scale, and shape parameters
#'   see \link{dnblomax} for details.
#' @param lik choice of negative binomial, geometric, or Poisson likelihood
#' @param q positive integer number of quadrature points, increase for better
#'   accuracy
#' @param add logical, should the curve be added to the existing plot window?
#' @param ... additional arguments passed to \link[graphics]{curve}
#'
#' @return a list with components x and y of the points that were drawn is
#'   returned invisibly.
#' @export
llcurve_lomax<-function(xmax, lpar=c(1,1,1), lik=c("nb","geom","poi"), q=1000,
                        add=TRUE, ...){
  #Draw the PMF curve on log-log axes for a Discrete-Lomax distribution
  #Curve goes from zero to xmax
  #lpar are the tail,scale,and shape parameters
  #q is number of quadrature points, increase for better accuracy
  lik<-match.arg(lik)
  if(lik=="geom"){ lpar[3]<-1 }
  if(lik=="poi"){ lpar[3]<-Inf }
  f<-function(t){
    dnblomax(floor(expm1(t)),tail=lpar[1],scale=lpar[2],shape=lpar[3],quadpts=q,log=TRUE)
  }
  curve(f,from=0,to=log1p(xmax),add=add,...)
}

#' @import statmod
nblomax_score<-function(x,tail=1,scale=1,shape=1,quadpts=1000,fix_shape=FALSE){
  #the gradient of the log likelihood (data=x) at each data point
  #returns matrix with length(x) rows and num_params cols
  #note if fix_shape=TRUE and shape=Inf this is a poisson-lomax model
  #can also force eg a Geometric by setting shape=1, fix_shape=TRUE
  r<-shape
  if(is.infinite(r)){ fix_shape<-TRUE } #poisson special case
  if(is.list(quadpts)){ #assumes quadpts based on Uniform(0,1)
    q<-quadpts
  } else { #recompute quadpts based on Uniform(0,1)
    q<-statmod::gauss.quad.prob(quadpts,dist="uniform",l=0,h=1)
    stop("error no quadpts provided!")
  }
  w<-log(q$weights) #vector of length m
  m<-length(w) #number of quadrature points
  mu<-scale*((1-q$nodes)^(-1/tail)-1) #vector of length m
  #ll a vector of length n representing marginal probability log p(x_i)
  ll<-dnblomax(x,tail=tail,scale=scale,shape=r,log=TRUE,quadpts=q)
  k<-if(fix_shape){ 2 } else { 3 }
  G<-matrix(0,nrow=m,ncol=k) #tail,scale,shape in cols
  vmu<-(mu+mu^2/r)
  G[,1]<- -log1p(mu/scale)*(mu+scale)/tail #dmu/dtail
  G[,2]<- mu/scale #dmu/dscale
  if(!fix_shape){ #part of gradient with respect to shape param
    G[,3]<- -digamma(r)+1-log1p(r*mu)
  }
  if(is.infinite(r)){ #poisson
    pxmu<-function(xi){dpois(xi,mu,log=TRUE)}
  } else { #neg binom
    pxmu<-function(xi){dnbinom(xi,size=r,mu=mu,log=TRUE)}
  }
  func<-function(xi,lli,fix_shape){ #returns a vector of length k
    #note this implicitly creates local copy of G
    G[,c(1,2)]<- G[,c(1,2)]*(xi-mu)/vmu #*standard GLM score function
    if(!fix_shape){
      G[,3]<- G[,3]+digamma(r+xi)-exp(log1p(r*xi)-log1p(r*mu))
    }
    #recycle p(x|mu)*p(mu)/p(x) importance weights
    #pxmu(xi) and w are vectors of length m, lli is a scalar
    #colSums(sign(G)*exp(log(abs(G))+pxmu(xi) - lli + w))
    colSums(G*exp(pxmu(xi) - lli + w)) #average over quadrature points
  }
  res<-t(mapply(func,x,ll,fix_shape)) #nxk
  res
}

#' @title Negative binomial-Lomax MLEs
#' @description Compute the maximum likelihood estimates for parameters of the
#'   negative binomial-Lomax distribution.
#' @name nblomax_mle
#'
#' @param x vector of non-negative integers (the data).
#' @param quadpts positive integer, number of quadrature points.
#' @param m numeric vector of length three. The prior means for the log of the
#'   tail, scale, and shape parameters. See \link{dnblomax} for details.
#' @param s numeric vector of length three. The prior standard deviations for
#'   the log of the tail, scale, and shape parameters. If s is finite, an L2
#'   (ridge) penalty shrinks the MLEs toward m. If s is Inf, the unpenalized MLE
#'   is obtained.
#' @param shape positive scalar. If NULL the shape parameter is estimated from
#'   data. If not NULL, the shape parameter is fixed at the specified value
#'   and only tail and scale are estimated. Shape=Inf corresponds to
#'   Poisson-Lomax.
#' @param om optimization method to be used by \link[stats]{optim}.
#' @param ... additional arguments passed to \link[stats]{optim}
#'
#' @return a named numeric vector containing the MLEs of tail, scale, and shape
#'   the maximized value of the log-likelihood is attached as an attribute.
#'
#' @import statmod
#' @export
nblomax_mle<-function(x, quadpts=1000, m=rep(0,3), s=rep(Inf,3), shape=NULL,
                      om="BFGS", ...){
  xt<- table(x)/length(x) #fraction of counts for each unique value of x
  xu<- as.integer(names(xt)) #unique values of x
  xt<- as.numeric(xt)
  q<-statmod::gauss.quad.prob(quadpts,dist="uniform",l=0,u=1)
  znames<-c("tail","scale","shape")
  if(!is.null(shape)){
    stopifnot(shape>0)
    ph<-TRUE
    k<-2
    getshape<-function(ez){shape}
    m<-m[1:2]; s<-s[1:2]; znames<-znames[1:2]
  } else {
    ph<-FALSE
    k<-3
    getshape<-function(ez){ez[3]}
  }
  #z = log of tail,scale,shape (initialize with all pars = 1)
  z<-rnorm(k,0,.5)
  names(z)<-znames #<-names(m)<-names(s)
  penalty<-1/(length(x)*s^2)
  f<-function(z){
    ez<-exp(z)
    res<-sum(xt*dnblomax(xu,tail=ez[1],scale=ez[2],shape=getshape(ez),log=TRUE,quadpts=q))
    res - sum(penalty*(z-m)^2)/2
  }
  g<-function(z){
    #gradient of f
    ez<-exp(z)
    score<-nblomax_score(xu,tail=ez[1],scale=ez[2],shape=getshape(ez),quadpts=q,fix_shape=ph)
    #score is length(xu) x k
    score<- colSums(xt*score)*ez #chain rule for logged parameters
    #score is a k-vector
    score - sum(penalty*(z-m))
  }
  res<-optim(z,f,g,control=list(fnscale=-1),method=om,...) #fnscale makes it maximize
  mle<-rep(NA,k)
  names(mle)<-znames
  if(res$convergence==0){
    mle<- exp(res$par)
    #attr(mle,"logscale_hessian")<-res$hessian
  }
  attr(mle,"loglik")<-length(x)*res$value
  mle
}

plomax_mle<-function(x,...){ nblomax_mle(x,shape=Inf,...) }

glomax_mle<-function(x,...){ nblomax_mle(x,shape=1,...) }

#' @title Negative binomial-Lomax MLEs for columns of a matrix
#' @description Compute the maximum likelihood estimates for parameters of the
#'   negative binomial-Lomax distribution for each column of a matrix
#' @name nblomax_mle_matrix
#'
#' @param m a matrix or sparse Matrix of non-negative integers (the data).
#' @param quadpts positive integer, number of quadrature points.
#' @param maxtry positive integer, number of times to retry fitting if initial
#'   fit fails due to numerical problems. Each retry increases quadpts by 500.
#' @param verbose logical, should progress messages be printed?
#' @param ... additional arguments passed to \link{nblomax_mle}
#'
#' @return a data frame whose rows correspond to the columns of m.
#'   The columns contain the MLEs for the parameters, the log-likelihood, and
#'   the \link[stats]{BIC} values.
#'
#' @import slam
#' @export
nblomax_mle_matrix<-function(m,quadpts=1000,maxtry=10,verbose=FALSE,...){
  #m a matrix of counts (0,1,2,...)
  #Computes the maximum likelihood estimates for neg binom-lomax fit
  #for each column of X
  #maxtry number of times to attempt fitting the MLE
  #... additional args passed to nblomax_mle
  mle_func<-function(x){
    mle<-nblomax_mle(x,quadpts=quadpts,...)
    c(mle,loglik=attr(mle,"loglik"))
  }
  if(is(m,"sparseMatrix")){
    apply_func<-function(m){
      m<-slam::as.simple_triplet_matrix(m)
      res<-slam::colapply_simple_triplet_matrix(m,mle_func)
      as.data.frame(do.call(rbind,res))
    }
  } else {
    apply_func<-function(m){ as.data.frame(t(apply(m,2,mle_func))) }
  }
  res<-apply_func(m)
  t<-1
  bad<-is.na(res[,1])
  while(t<maxtry && any(bad)){
    nbad<-sum(bad)
    if(verbose){ print(paste0("iteration: ",t,"; number bad entries: ",nbad)) }
    quadpts<-quadpts+500
    if(nbad>1){
      res2<-apply_func(m[,bad])
    } else {
      res2<-mle_func(m[,bad])
    }
    res[bad,]<-res2
    t<-t+1
    bad<-is.na(res[,1])
  }
  #res<-as.data.frame(res)
  res$bic<- -2*res$loglik + log(nrow(m))*(ncol(res)-1)
  res
}

plomax_mle_matrix<-function(m,...){ nblomax_mle_matrix(m,shape=Inf,...) }

glomax_mle_matrix<-function(m,...){ nblomax_mle_matrix(m,shape=1,...) }

#' @title Estimate Poisson- or geometric-Lomax scale parameter from zero
#'   fraction
#' @description Given a fixed tail parameter, use the fraction of zeros to
#'   estimate the scale parameter of a Poisson-Lomax or geometric-Lomax
#'   distribution via the method of moments.
#' @name pglomax_pzero2scale
#'
#' @param lpz a numeric vector representing the natural log of the empirical
#'   fraction of zeros from count data.
#' @param lik string specifying either the Poisson or geometric likelihood.
#' @param tail either a positive scalar or a vector of same length as lpz
#'   specifying the fixed tail parameter(s).
#' @param quadpts positive integer, number of quadrature points.
#' @param lims interval of possible values for the scale
#'   parameter, passed to \link[stats]{uniroot}.
#'
#' @return a numeric vector of estimated scale parameters for each element of
#'   lpz.
#'
#' @export
pglomax_pzero2scale<-function(lpz, lik=c("poisson","geometric"), tail=0.8,
                              quadpts=1000, lims=c(1e-6,100)){
  #Assuming the data follows a Poisson-Lomax
  #or Geometric-Lomax distribution,
  #if we fix the tail parameter to a specified value
  #we can infer the scale parameter from the fraction of zeros
  #lpz is a vector of the log of zero fraction for each cell
  #quadpts is number of quadrature points to improve accuracy of estimating the PMF
  #lims are the lower,upper bounds for the scale parameter passed to uniroot
  lik<-match.arg(lik)
  dfunc<-if(lik=="poisson"){ dplomax } else { dglomax }
  inner<-function(x,t){
    #x is an element of lpz
    f<-function(s){
      #s is the estimated scale parameter
      x-dfunc(0,tail=t,scale=s,log=TRUE,quadpts=quadpts)
    }
    uniroot(f,lims)$root
  }
  if(length(tail)==1){ #single tail parameter for all cells
    return(vapply(lpz,inner,FUN.VALUE=1.0,t=tail))
  } else { #different tail parameter for each cell
    return(mapply(inner,lpz,tail))
  }
}

################### simulation functions #########################

rburr<-function(n,lambda,k=1.01,shape=1,is_median=TRUE,logscale=FALSE){
  #returns random variates from the Burr or Singh-Maddala Distribution
  #special case of shape=1 corresponds to Lomax distribution
  #these are a special case of the beta prime distribution
  #lambda=median if is_median=TRUE
  #if is_median=FALSE lambda is the scale parameter (larger=larger median)
  #k=power law tail parameter. smaller values=heavier tail
  #shape=concentration parameter. Larger values= peakier distribution
  #for c<=1 the mode is zero. Otherwise mode is nonzero
  #the mean only exists if k>1, variance only exists if k>2
  log_pareto<-log(expm1(rexp(n,k)))
  log_lambda<-log(lambda)
  if(is_median) log_lambda<-log_lambda-log(expm1(log(2)/k))/shape
  lx<-log_lambda+log_pareto/shape
  if(logscale) return(lx)
  return(exp(lx))
}
rlomax<-function(n,tail=1.01,scale=1,logscale=FALSE,cut=Inf){
  x<-rburr(n,scale,k=tail,shape=1,is_median=FALSE,logscale=logscale)
  #x<-VGAM::rlomax(n,scale=scale,shape3.q=tail)
  #if(logscale){ x<-log(x) }
  if(is.infinite(cut)){
    return(x)
  } else {
    #exponential cutoff via competing risks
    #set intersection of survival curves to "cut" value, eg 10^6
    y<-rexp(n,rate=tail*log1p(cut/scale)/cut)
    if(logscale) y<-log(y)
    return(pmin(x,y))
  }
}
rplomax<-function(n,tail=1.01,scale=1,cut=Inf){
  #poisson with lomax prior
  mu<-rlomax(n,tail=tail,scale=scale,logscale=FALSE,cut=cut)
  rpois(n,mu)
}
rglomax<-function(n,tail=1.01,scale=1,cut=Inf){
  lmu<-rlomax(n,tail=tail,scale=scale,logscale=TRUE,cut=cut)
  rgeom(n,plogis(lmu,lower.tail=FALSE))
}
rnblomax<-function(n,tail=1.01,scale=1,shape=1,cut=Inf){
  mu<-rlomax(n,tail=tail,scale=scale,logscale=FALSE,cut=cut)
  rnbinom(n,size=shape,mu=mu)
}
