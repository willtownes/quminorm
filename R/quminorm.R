#Quasi-UMIs:
#estimate scale parameter using only the fraction of zeros
#quantile normalization to a Poisson-Lomax or Geometric-Lomax

#' @importFrom utils tail
make_cdf_nz<-function(thresh,dfunc,maxval=1e6){
  #dfunc=some log-pmf function accepting a single argument (the data)
  #let cdf be the cdf corresponding to dfunc
  #thresh is a probability value, where if cdf>thresh,...
  #...it implies 1-cdf_nz is less than 1/number of nonzero data points
  #beyond this point the quantile function resolves to less than one data point
  #so no point in computing cdf above threshold
  #VALUE: cdf_nz: the zero-truncated version of cdf
  #cdf_nz can then be used to quantile normalize data to the distribution of dfunc
  lthresh<-log(thresh)
  lo<-0; hi<-100
  lcdf<-NULL
  lcdf_max<- -Inf
  ctr<-0 #counter
  while(lcdf_max<=lthresh && hi<maxval){
    ctr<-ctr+1
    z<-seq.int(lo,hi-1) #0:99, 100:199, 200:399, ...
    lpmf<-dfunc(z)
    #update the first pmf value by adding the max cdf to it
    lpmf[1]<-logspace_add(lpmf[1],lcdf_max)
    lcdf<-c(lcdf,log_cumsum_exp(lpmf))
    lcdf_tail<-tail(lcdf,2) #make sure not to change the "2", it is crucial!
    if(diff(lcdf_tail)==0){
      stop("CDF is not increasing, PMF function may have numerical problems!")
    }
    lcdf_max<-lcdf_tail[2] #lcdf_tail must be length 2!!
    lo<-hi #100
    hi<-2*hi #200
  }
  if(hi>maxval){
    stop("Exceeded max value, PMF function may have numerical problems!")
  }
  #remove extra elements of cdf that extend beyond the threshold
  k<-min(which(lcdf>lthresh,arr.ind=TRUE))
  #renormalize CDF for nonzero values only
  #cdf[1] is P(x=0)
  pmf0<-exp(lcdf[1])
  (exp(lcdf[2:k])-pmf0)/(1-pmf0)
}

quminorm_inner<-function(xnz,cdf_nz,nnz=length(xnz)){
  #xnz a vector of positive integers
  #cdf_nz a CDF for some zero-truncated distribution
  #value: the quantile normalized version of xnz
  rk_th<-ceiling(cdf_nz*nnz) #theoretical ranks, increasing order
  dups<-duplicated(rk_th) #TRUE occurs for first duplicate, then all FALSE
  targets<-seq_along(rk_th)[!dups] #the actual qumi values we map to.
  rk_th<-rk_th[!dups]
  xrk<-rank(xnz,ties.method="min") #convert data to empirical ranks
  #xmap gives indices within targets where the quantile function points to
  xmap<-findInterval(xrk,rk_th,left.open=TRUE)+1 #intervals are (lo,hi]
  targets[xmap]
}

#' @importFrom sads dpoilog
quminorm_poilog_nz<-function(xnz,n,shape,sc=NULL,err2na=TRUE){
  #xnz: a data vector with only the positive counts
  #n: the length of the data vector if the zero counts were also included
  #the number of zeros is implied to be n-length(xnz)
  #shape=sig, sc=mu
  #sig,mu are params of lognormal as in sads::dpoilog
  #returns a quantile normalized version of the positive data xnz
  nnz<-length(xnz)
  #stopifnot(n>nnz)
  if(is.null(sc)){
    #assumes at least one gene is a zero count
    lpz<-log(n-nnz)-log(n) #log(fraction of zeros)
    sc<-poilog_pzero2mu(lpz,sig=shape)
  }
  dfunc<-function(x){ dpoilog(x,mu=sc,sig=shape,log=TRUE) }
  pmf0<-exp(dfunc(0))
  #threshold for cdf on the regular scale such that
  #zero truncated cdf (cdf_nz) extends beyond the 1-1/nnz threshold
  thresh<-(1-1/nnz)*(1-pmf0)+pmf0
  if(err2na){
    cdf_nz<-tryCatch(make_cdf_nz(thresh,dfunc),error=function(e){NULL})
  } else {
    cdf_nz<-make_cdf_nz(thresh,dfunc)
  }
  if(is.null(cdf_nz)){
    return(rep(NA,nnz))
  } else {
    return(quminorm_inner(xnz,cdf_nz,nnz))
  }
}

quminorm_poilog<-function(x,shape,...){
  #x a data vector with at least one zero value and the rest positive counts
  #... additional arguments passed to quminorm_poilog_nz
  #returns a quantile normalized version of the data x with zeros unchanged
  nzi<-which(x>0)
  x[nzi]<-quminorm_poilog_nz(x[nzi],length(x),shape,...)
  x
}

quminorm_dense<-function(m,shape,mc.cores=1){
  #this function is best when m is dense and doesn't contain many zeros
  if(mc.cores==1){
    for(i in seq_len(ncol(m))){
      m[,i]<-quminorm_poilog(m[,i],shape)
    }
  } else {
    res<-colapply_dense_parallel(m,quminorm_poilog,shape,mc.cores=mc.cores)
    m[,]<-as.numeric(unlist(res))
  }
  m
}

#' @importFrom methods as
#' @importFrom Matrix drop0
quminorm_sparse<-function(m,shape,mc.cores=1){
  #if m is sparse, it's better to operate only on nonzero elements
  m<-as(drop0(m), "CsparseMatrix") #make sure it is column-oriented sparsity
  qmlist<-colapply_sparse_nonzero(m, quminorm_poilog_nz, n=nrow(m),
                           shape=shape, mc.cores=mc.cores)
  #replace all nonzero elements with stacked cols vector
  m@x<-as.numeric(unlist(qmlist))
  m #return sparse matrix
}
