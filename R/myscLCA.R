#' LCA scRNA analysis
#'
#' This function loads a gene by cell transcript count matarix (row by column).
#'
#' @param datmatrix data matrix, gene by cell
#' @param clust.max maximum number of clusters specified by user, default: 10
#' @param trainingSetSize the number of cells used in training set; default: 1000; maximum: the total number of cells in data matrix
#' @param datBatch a vector of batch group ID for each cell, default: NULL
#' @return an integer vector showing clustering membership
#' @examples
#' # Load example dataset
#' data(myscExampleData)
#' 
#' # Both data matrix and true labels of cells are provided
#' names(myscExampleData)
#' 
#' # With 14,074 genes and 250 cells, 83% of entries are zero
#' dim(myscExampleData$datamatrix)
#' 
#' # Three types of cells
#' table(myscExampleData$truelabel)
#' 
#' # Start my scRNAseq LCA analysis
#' myclust.res <- myscLCA(myscExampleData$datamatrix)
#' 
#' # Top results are provide in a list
#' length(myclust.res)
#' 
#' # The result ranked as the best, compared with the true labels:
#' table(myclust.res[[1]],myscExampleData$truelabel)
#' 
#' # The result ranked as the second best:
#' table(myclust.res[[2]],myscExampleData$truelabel)
#' 
#' # The result ranked as the third best:
#' table(myclust.res[[3]],myscExampleData$truelabel)
#' 
#' @export 

myscLCA <- function(datmatrix, clust.max=10, trainingSetSize=1000, datBatch=NULL){
  myRes <- list()
  myResSilScore <- c()
  spec <- T
  if( ncol(datmatrix) > trainingSetSize){
    mytrainset <- sample(ncol(datmatrix),trainingSetSize)
    TPM <- datmatrix[,mytrainset]
    myBatch <- datBatch[mytrainset]
  }else{
    TPM <- datmatrix
    myBatch <- datBatch
  }
  TPMselect <- which(rowSums(TPM) > 0)
  TPM <- TPM[TPMselect,]
  UMI <- colSums(TPM)
  TPM.log <- log2(TPM / matrix(UMI, ncol=ncol(TPM), nrow=nrow(TPM), byrow=T) + 2.5e-5)
  pca <- prcomp(TPM.log, scale=T)
  pcaT <- prcomp(t(TPM.log), scale=T)
  nTWT <- which(TWtest(pcaT$sdev)$pvalue <= 0.05)
  nTWT <- nTWT[which(abs(1-cosDist.part(rbind(log2(UMI) - mean(log2(UMI)), t(pcaT$x[,nTWT])), idx1 = 1, idx2 = -1)) < 0.5)]
  if(!is.null(datBatch)){
    nTWT <- nTWT[ which( p.adjust(apply(pcaT$x[,nTWT],2,function(x){kruskal.test(x,myBatch)$p.value}),"fdr") > 0.05) ]
  }
  nTW <- which(TWtest(pca$sdev)$pvalue <= 0.05)
  nTW <- nTW[which(abs(1-cosDist.part(rbind(log2(UMI), t(pca$rotation[,nTW])), idx1 = 1, idx2 = -1)) < 0.5)]
  if(!is.null(datBatch)){
    nTW <- nTW[ which( p.adjust(apply(pca$rotation[,nTW],2,function(x){kruskal.test(x,myBatch)$p.value}),"fdr") > 0.05) ]
  }

  pcaT.dis <- cosDist(pcaT$x[,nTWT])
  factors.best <- nTW
  cutoff <- quantile(apply(pca$rotation[,factors.best], 1, sumsq), 0.975)
  outliers = which(apply(pca$rotation[,factors.best], 1, sumsq) > cutoff)
  rot2 = pca$rotation[-outliers,factors.best]
  myDist.best = cosDist(pca$rotation[,factors.best]) 
  if (spec) {
    A = 1 - myDist.best / 2
    D = diag(apply(A, 1, sum))
    L <- (D %^% (-1/2)) %*% A %*% (D %^% (-1/2))
    evL <- eigen(L, symmetric=TRUE)
    nF = clust.max
    if (nF < 10) nF = 10
    if (nF > 50) nF = 50
    X = evL$vectors[-outliers,1:nF]
    X = X/sqrt(apply(X, 1, sumsq))
    myHClust2 = hclust(dist(X), method="average")
  }
  myHClust = hclust(as.dist(myDist.best[-outliers, -outliers]), method="average")
  sil.best = NA
  sil.best2 = rep(NA, clust.max)
  clusters = matrix(NA, nrow=clust.max, ncol=ncol(TPM)) 
  for (nClust in 2: clust.max) {
    clust1 = cutree(myHClust, nClust) # inital seed
    for (z in 1:1000) {
      my.center = matrix(0, nrow=nClust, ncol=length(factors.best))
      for (i in 1:nClust)
        if (sum(clust1 == i) > 1) {
          my.center[i,] = colMeans(rot2[which(clust1==i),])
        } else if (sum(clust1 == i) == 1) {
          my.center[i,] = rot2[which(clust1==i),]
        }
      dist2 = cosDist.part(rbind(rot2, my.center), idx1 = length(clust1) + 1:nClust, idx2 = 1:length(clust1))
      clust2 = apply(dist2, 2, which.min)
      if (sum(clust1 != clust2) == 0) break
      clust1 = clust2
    }
    if (spec) {
      clust1 = cutree(myHClust2, nClust)# inital seed
      mycenter = matrix(NA, nrow=nClust, ncol=ncol(X))
      for (i in 1:nClust) {
        if (sum(clust1 == i) > 1) {
          mycenter[i,] = colMeans(X[which(clust1==i),])
        } else {
          mycenter[i,] = X[which(clust1==i),]
        }
      }
      clust3 = kmeans(X, mycenter)$cluster
      if (mean(silhouette(clust3, myDist.best[-outliers,-outliers])[,3]) > mean(silhouette(clust2, myDist.best[-outliers,-outliers])[,3])) clust2 = clust3
    }
    clust1 = array(NA, ncol(pca$rotation))
    clust1[-outliers] = clust2
    clust2 = clust1
    for (i in outliers) {
      tmp = array(NA, nClust)
      for (j in 1:nClust) tmp[j] = mean(myDist.best[i, which(clust2 == j)])
      clust1[i] = which.min(tmp)
    }
    sil.best2[nClust] = mean(silhouette(clust1, pcaT.dis)[,3])
    clusters[nClust,] = clust1
    if (is.na(sil.best) | sil.best2[nClust] > sil.best) {
      sil.best = sil.best2[nClust]
      clust.best = clust1
    }
  }
  if(sil.best <0){
    mymembership <- rep(1,ncol(datmatrix))
    myRes[[1]] <- mymembership
    myResSilScore[1] <-sil.best
  }else{
    myTop3 <- max( maxN(sil.best2[-1],3), 0)
    myTop3 <- which(sil.best2 >= myTop3)
    for( xi in 1:length(myTop3)){
      clust.best <- clusters[myTop3[xi],]
      pval = matrix(NA, nrow=myTop3[xi], ncol=max(nTW))
      
      for (j in 1:nrow(pval)) for (i in nTW) pval[j,i] = wilcox.test(pca$rotation[which(clust.best==j), i], pca$rotation[which(clust.best!=j), i])$p.value
      factors.best = which(apply(pval, 2, min) < 0.05/(ncol(pval)-1)/nrow(pval))
      nClust = nrow(pval)
      if (length(factors.best) >= 2) {
        cutoff = quantile(apply(pca$rotation[,factors.best], 1, sumsq), 0.975)
        outliers = which(apply(pca$rotation[,factors.best], 1, sumsq) > cutoff)
        rot2 = pca$rotation[-outliers,factors.best]
        myDist.best = cosDist(pca$rotation[,factors.best]) 
        if (spec) {
          A = 1 - myDist.best / 2
          D = diag(apply(A, 1, sum))
          L <- (D %^% (-1/2)) %*% A %*% (D %^% (-1/2))
          evL <- eigen(L, symmetric=TRUE)
          nF = length(factors.best)
          if (nF > 50) nF = 50
          X = evL$vectors[-outliers,1:nF]
          X = X/sqrt(apply(X, 1, sumsq))
          myHClust2 = hclust(dist(X), method="average")
        }
        myHClust = hclust(as.dist(myDist.best[-outliers, -outliers]), method="average")
        clust1 = cutree(myHClust, nClust) # inital seed
        for (z in 1:1000) {
          my.center = matrix(0, nrow=nClust, ncol=length(factors.best))
          for (i in 1:nClust)
            if (sum(clust1 == i) > 1) {
              my.center[i,] = colMeans(rot2[which(clust1==i),])
            } else if (sum(clust1 == i) == 1) {
              my.center[i,] = rot2[which(clust1==i),]
            }
          dist2 = cosDist.part(rbind(rot2, my.center), idx1 = length(clust1) + 1:nClust, idx2 = 1:length(clust1))
          clust2 = apply(dist2, 2, which.min)
          if (sum(clust1 != clust2) == 0) break
          clust1 = clust2
        }
        if (spec) {
          clust1 = cutree(myHClust2, nClust)# inital seed
          mycenter = matrix(NA, nrow=nClust, ncol=ncol(X))
          for (i in 1:nClust) {
            if (sum(clust1 == i) > 1) {
              mycenter[i,] = colMeans(X[which(clust1==i),])
            } else {
              mycenter[i,] = X[which(clust1==i),]
            }
          }
          clust3 = kmeans(X, mycenter)$cluster
          if (mean(silhouette(clust3, myDist.best[-outliers,-outliers])[,3]) > mean(silhouette(clust2, myDist.best[-outliers,-outliers])[,3])) clust2 = clust3
        }
        clust.best = array(NA, ncol(pca$rotation))
        clust.best[-outliers] = clust2
        clust2 = clust.best
        for (i in outliers) {
          tmp = array(NA, nClust)
          for (j in 1:nClust) tmp[j] = mean(myDist.best[i, which(clust2 == j)])
          clust.best[i] = which.min(tmp)
        }
      } else if (length(factors.best) == 1) {
        myDist.best = as.matrix(dist(pca$rotation[,factors.best]))
        myHClust = hclust(dist(pca$rotation[-outliers,factors.best]), method="average")
        X = pca$rotation[-outliers,factors.best]
        clust1 = cutree(myHClust, nClust)# inital seed
        mycenter = array(NA, nClust)
        for (i in 1:nClust) mycenter[i] = mean(X[which(clust1 == i)])
        clust1 = kmeans(X, mycenter)$cluster
        clust.best = array(NA, ncol(pca$rotation))
        clust.best[-outliers] = clust1
        clust2 = clust.best
        for (i in outliers) {
          tmp = array(NA, nClust)
          for (j in 1:nClust) tmp[j] = mean(myDist.best[i, which(clust2 == j)])
          clust.best[i] = which.min(tmp)
        }
        
      }
      
      if( ncol(datmatrix) > trainingSetSize){
        mymembership <- clust.best
        
        myNumClust <- length(unique(clust.best))
        if(length(factors.best) > 1){
          center.best <- t(sapply(1:myNumClust,function(x){colMeans(pca$rotation[which(mymembership==x),factors.best])}))
        }else{
          center.best <- matrix( tapply(pca$rotation[,factors.best],mymembership,mean), ncol=1)
        }
        mysvd.d <- length(TPMselect-1)*pca$sdev[factors.best]^2
        mytrans.mat <- (1/mysvd.d)*t(pca$x[,factors.best])
        TPM <- datmatrix[TPMselect,]
        UMI <- colSums(TPM)
        TPM.log <- log2(TPM / matrix(UMI, ncol=ncol(TPM), nrow=nrow(TPM), byrow=T) + 2.5e-5)
        TPM.scale <- scale(TPM.log)
        myQsize <- ncol(TPM.scale)
        myDsize <- nrow(center.best)
        myQueries <- mytrans.mat %*% TPM.scale
        myassignment <- cosDist.part(rbind(t(myQueries),center.best),idx1=myQsize+1:myDsize,idx2=1:myQsize)
        myQmembership <- apply(myassignment[,1:myQsize],2,which.min)
        
        mymembership <- myQmembership
      }else{
        mymembership <- clust.best
      }
      myRes[[xi]] <- mymembership
      myResSilScore[[xi]]=sil.best2[myTop3[xi]]
    }#xi
  }
  myRes <- myRes[order(myResSilScore,decreasing = T)]
  return(myRes);
}
