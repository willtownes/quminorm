####### Miscellaneous utility functions ##########

####### logsumexp functions #######

#x,y must be scalars!
#' @importFrom matrixStats logSumExp
logspace_add<-function(x,y){ logSumExp(c(x,y)) }

log_cumsum_exp<-function(lp){
    #lp a vector of log probabilities
    #returns the log of the cumulative sum of the probabilities
    Reduce(logspace_add,lp,accumulate=TRUE)
}

####### sparse Matrix colapply functions #######

#' @importFrom slam as.simple_triplet_matrix colapply_simple_triplet_matrix
colapply_full<-function(X,FUN,...){
    #apply a function f to each column of sparse Matrix X
    #FUN is applied to the entire column, including the zero values
    #for a faster alternative that only operates on nonzero values, see _____
    #if FUN returns a scalar the result is a vector of length ncol(m)
    #otherwise the result is a list of length ncol(m)
    #... additional args passed to f
    colapply_simple_triplet_matrix(as.simple_triplet_matrix(X),FUN,...)
}

#' @importFrom Matrix nnzero
listCols<-function(m){
    #converts a sparse Matrix into a list of its columns
    #each list item contains only the nonzero elements of the column
    res<-split(m@x, findInterval(seq_len(nnzero(m)), m@p, left.open=TRUE))
    names(res)<-colnames(m)
    res
}

#' @importFrom parallel mclapply
colapply_nonzero<-function(X,FUN,...,mc.cores=1){
    #apply a function FUN to NONZERO elements of each column of sparse Matrix X
    #for an alternative that operates on all values, see ____
    #mc: should parallel processing be used? Only recommended if FUN is slow
    #... additional args passed to mclapply or to FUN
    #this function always returns a list of length ncol(X)
    if(mc.cores==1){
        return(lapply(listCols(X),FUN,...))
    } else {
        return(mclapply(listCols(X),FUN,...,mc.cores=mc.cores))
    }
}

