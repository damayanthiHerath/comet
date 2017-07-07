library("fpc")
library("factoextra")

#The datadf should be inthe form, data.frame(contigId,GC content,log(coverage),contiglength,tetranulceotidefrequencyvalues)
#The eps1 and minpnts1 refers to the parameters for identifying core,border, noise points. 
#dispparamtier2 referes to the minimum distance of the points in a cluster in DPGMM.
#choices refers to the number of principal components, tetranucleotide frequencies are being reduced to. 
binTwoTiereddpanddbscan<- function(datadf,eps1,minpnts1,dispparamtier2,choices){
  
  #Tier 1 binning considering GC Content and Coverage
  dbscan::kNNdistplot(datadf[,2:3], k =  minpnts1)
  db <- fpc::dbscan(datadf[,2:3], eps = eps1, MinPts = minpnts1,scale=TRUE)
  print(fviz_cluster(db, datadf[,2:3], stand = FALSE, geom = "point"))
  classification.1<- data.frame('id'=data$id,'label'=db$cluster)
  
  tier1labels<- unique(classification.1$label)
  
  tier2summary= data.frame()
  
  disp.param2<- dispparamtier2
  for (i in 1:length(tier1labels)) {
   
  (v1.class <- tier1labels[i])
   toclust.subset <- data[classification.1$label == v1.class,]
    #skip the outliers
    if(v1.class==0) {
      ids = toclust.subset[,1]
      classification.2<- data.frame('id'=ids,'label-2'='0')
   }
    

    else{
      
      features <- 6:(ncol(toclust.subset)-1)
      toclust <- toclust.subset[,features]
      pca <- prcomp(toclust,scale=T)
     
      toclust.pca <- pca$x[,choices]
      ids = toclust.subset[,1]
      dp.update <- dirichletClusters(toclust.pca, disp.param2)
      classification.2<- data.frame('id'=ids,'label-2'=dp.update$cluster)
      plot(toclust.pca[,1],toclust.pca[,2],col=dp.update$cluster,xlab = 'PC1',ylab='PC2',pch=19)
    }
    tier2summary<- rbind(tier2summary,classification.2)
    
  }
  ans<- merge(classification.1,tier2summary,by.x = "id",by.y = "id")
  ans$bin <- paste(ans$label,ans$label.2,sep = "")
 
  return(ans)
  
}
#The DPGMM clustering method was adopted from the R implementation
#https://statistical-research.com/index.php/2013/04/07/dirichlet-process-infinite-mixture-models-and-clustering/

dirichletClusters <- function(orig.data, disp.param = NULL, max.iter = 100, tolerance = .001)
{
  n <- nrow( orig.data )
  
  data <- as.matrix( orig.data )
  pick.clusters <- rep(1, n)
  k <- 1
  
  mu <- matrix( apply(data,2,mean), nrow=1, ncol=ncol(data) )
  
  is.converged <- FALSE
  iteration <- 0
  
  ss.old <- Inf
  ss.curr <- Inf
  
  while ( !is.converged & iteration < max.iter ) { # Iterate until convergence
    iteration <- iteration + 1
    
    for( i in 1:n ) { # Iterate over each observation and measure the distance each observation' from its mean center's squared distance from its mean
      distances <- rep(NA, k)
      
      for( j in 1:k ){
        distances[j] <- sum( (data[i, ] - mu[j, ])^2 ) # Distance formula.
      }
      
      if( min(distances) > disp.param ) { # If the dispersion parameter is still less than the minimum distance then create a new cluster
        k <- k + 1
        pick.clusters[i] <- k
        mu <- rbind(mu, data[i, ])
      } else {
        pick.clusters[i] <- which(distances == min(distances))
      }
      
    }
    
    ##
    # Calculate new cluster means
    ##
    for( j in 1:k ) {
      if( length(pick.clusters == j) > 0 ) {
        mu[j, ] <- apply(subset(data,pick.clusters == j), 2, mean)
      }
    }
    
    ##
    # Test for convergence
    ##
    ss.curr <- 0
    for( i in 1:n ) {
      ss.curr <- ss.curr +
        sum( (data[i, ] - mu[pick.clusters[i], ])^2 )
    }
    ss.diff <- ss.old - ss.curr
    ss.old <- ss.curr
    if( !is.nan( ss.diff ) & ss.diff < tolerance ) {
      is.converged <- TRUE
    }
    
  }
  
  centers <- data.frame(mu)
  ret.val <- list("centers" = centers, "cluster" = factor(pick.clusters),
                  "k" = k, "iterations" = iteration)
  
  return(ret.val)
}
