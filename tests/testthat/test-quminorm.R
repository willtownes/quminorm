context("test quminorm and poilog functions")
library(Matrix)
library(SingleCellExperiment)
set.seed(1234)

#set up some test data
ncells<-3
ngene<-2000
#mean_tot_reads<-1e5
mean_tot_umi<-1000
lambda<-matrix(exp(rnorm(ncells*2000,0,2)),ncol=ncells)
n<-rpois(ncells,mean_tot_umi)
lambda<-t(t(lambda)*n/colSums(lambda))
#lambda<-lambda*rep(mu_cell,each=ngene)
m0<-matrix(rpois(length(lambda),lambda),ncol=ncells)
pcr<-1+rgamma(length(m0),2,rate=2) # mean=1+1=2, minimum=1
m<-round(m0*pcr)
tpm<-t(t(m)/colSums(m))*1e6
m2<-Matrix(m,sparse=TRUE)
tpm2<-Matrix(tpm,sparse=TRUE)
se<-SummarizedExperiment(assays=list(readcounts=m))
sce<-SingleCellExperiment(assays=list(readcounts=m,readcounts_sparse=m2,
                                      tpm=tpm,tpm_sparse=tpm2))

test_that("poilog mle can be estimated from data",{
    fit1<-poilog_mle(m0)
    expect_gt(min(fit1$sig),0) #sigma param must be positive
    expect_equal(colnames(fit1),c("mu","sig","loglik","bic"))
    expect_is(fit1,"data.frame")
    expect_equal(dim(fit1),c(ncol(m0),4))
    m0s<-Matrix(m0,sparse=TRUE)
    expect_equal(fit1,poilog_mle(m0s)) #check sparse and dense give same result
})

test_that("quminorm preserves sparsity and rank-ordering of features", {
    q0<-quminorm(m)
    expect_gte(min(q0),0) #no negative values
    #preserve rank-ordering
    ranks_original<-rank(m[,1],ties.method="min")
    ranks_qumi<-rank(q0[,1],ties.method="min")
    expect_true(all(ranks_qumi<=ranks_original))
    expect_equal(q0>0, m>0) #preserve sparsity
})

test_that("quminorm handles different input types properly", {
    shp<-1.5
    aname<-paste0("qumi_poilog_",shp)
    #check no error with different object types
    q1<-quminorm(m,shape=shp)
    expect_is(q1,"matrix") #check correct output type
    q2<-quminorm(m2,shape=shp)
    expect_s4_class(q2,"CsparseMatrix") #check correct output type
    expect_equivalent(q1,as.matrix(q2)) #check consistent output
    se2<-quminorm(se,"readcounts",shape=shp)
    expect_s4_class(se2,"SummarizedExperiment")
    expect_equal(assay(se2,aname),q1)
    sce2<-quminorm(sce,"tpm_sparse",shape=shp)
    expect_s4_class(sce2,"SingleCellExperiment")
    expect_equal(assay(sce2,aname),q2)
})
