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

####### dense matrix parallel colapply functions #######

listCols_dense<-function(X){
    #https://stackoverflow.com/questions/6819804/how-to-convert-a-matrix-to-a-list-of-column-vectors-in-r
    split(X,col(X,as.factor=TRUE))
}

#' @importFrom parallel mclapply
colapply_dense_parallel<-function(X,FUN,...,mc.cores=2){
    mclapply(listCols_dense(X),FUN,...,mc.cores=mc.cores)
}

####### sparse Matrix colapply functions #######

#' @importFrom slam as.simple_triplet_matrix colapply_simple_triplet_matrix
colapply_sparse_full<-function(X,FUN,...){
    #apply a function f to each column of sparse Matrix X
    #FUN is applied to the entire column, including the zero values
    #for a faster alternative that only operates on nonzero values,...
    #...see colapply_nonzero
    #if FUN returns a scalar the result is a vector of length ncol(m)
    #otherwise the result is a list of length ncol(m)
    #... additional args passed to f
    colapply_simple_triplet_matrix(as.simple_triplet_matrix(X),FUN,...)
}

#' @importFrom methods as
#' @importFrom Matrix nnzero
listCols_sparse<-function(X){
    #converts a sparse Matrix into a list of its columns
    #each list item contains only the nonzero elements of the column
    X<-as(X,"CsparseMatrix")
    res<-split(X@x, findInterval(seq_len(nnzero(X)), X@p, left.open=TRUE))
    names(res)<-colnames(X)
    res
}

#' @importFrom parallel mclapply
colapply_sparse_nonzero<-function(X,FUN,...,mc.cores=1){
    #apply a function FUN to NONZERO elements of each column of sparse Matrix X
    #for an alternative that operates on all values, see colapply_full
    #mc: should parallel processing be used? Only recommended if FUN is slow
    #... additional args passed to mclapply or to FUN
    #this function always returns a list of length ncol(X)
    if(mc.cores>1){
        return(mclapply(listCols_sparse(X),FUN,...,mc.cores=mc.cores))
    } else {
        return(lapply(listCols_sparse(X),FUN,...))
    }
}

