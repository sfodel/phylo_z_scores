phylo_z_scores <- function(phy, OTU, ncores, nreps, nfactors, taxonomy) {
  ##phy = a phylogenetic tree in the form of a phylo object (e.g. imported .nwk file with read.tree of package ape)
  ##OTU = a "species" (e.g. Sequence Variants) by communities table with raw sequence counts or cell counts. Row names and column names necessary for "species" (e.g. SV hashes) and communities, respectively
  ##ncores = the number of desired cores to use for parallelization. Strongly suggested to use at least 20 cores, because the standard deviation of the null distribution depends on that.
  ##nreps = the number of permutations to use when creating the null distributions for the calculation of the z-scores. A minimum of 100 is recommended.
  ##nfactors = the number of phylogenetic factors to search for in phylofactorization. If the p-values become non-significant before reaching this number, the script stops automatically, reporting up to that number of factors
  ##taxonomy = a two column table with "species" names (e.g. SV hashes) in the first column and taxonomic assignments in the second.
  
  ###Output is an S4 object that has the following objects collated:
    ### 1.@species_cluster: A data.table of the individual species scores and taxonomy
    ### 2.@constant : A data.table containing only species with z-scores < -2
    ### 3.@taxa.weights : A collapsed data.table of @constant per unique taxonomies
    ### 4.@phylofactorization : The phylofactor object performed on the sum_score of @species_cluster for nfactors
    ### 5.@phylo_tree : The ggtree visualization of the factors of @phylofactorization
  
  library(phylofactor)
  library(parallel)
  library(picante)
  library(ggtree)
  
  timestamp()
  print("Loading inputs")
  OTU <- apply(OTU,2,FUN=function(x) x/sum(x))
  species <- rownames(OTU)
  phydist <- cophenetic(phy)
  pdist <- phydist[species,species]
  rm(phydist)
  taxonomy <- as.data.table(taxonomy)
  names(taxonomy) <- c("species", "taxonomy")
  
  ##Make data.table ix with rows like the upper diagonal of your (community*community) matrix
  timestamp()
  print("Making initial community-community data.table")
  ix <- data.table('i'=rep(1:ncol(OTU),times=ncol(OTU):1-1))
  ix[,j:=(i+1):ncol(OTU),by=i]
  
  ix[,i:=colnames(OTU)[i]] #give names to community i in the table
  ix[,j:=colnames(OTU)[j]] #give names to community j in the table
  ###Optional - Remove community pairs from column name field that are separated by underscore
  #ix=ix[paste(sapply(strsplit(i,'_'),getElement,1), sapply(strsplit(i,'_'),getElement,2))!=paste(sapply(strsplit(j,'_'),getElement,1), sapply(strsplit(j,'_'),getElement,2))]
  ##For us, this removes community pairs from within the same patch (6*19)+4 = 118 community pairs because we lack one replicate in one stream
  
  ##Load sub-functions
  
  NTD_calc <- function(i,j,spp,y,pdist,OTU){
    
    if (y[i]>0 & y[j]==0){
      other.species <- setdiff(colnames(pdist),spp)
      other.species <- intersect(other.species,rownames(OTU)[which(OTU[,j]>0)])
      min.dist <- min(pdist[spp,other.species])
      a <- y[i]
    } else if (y[j]>0 & y[i]==0){
      other.species <- setdiff(colnames(pdist),spp)
      other.species <- intersect(other.species,rownames(OTU)[which(OTU[,i]>0)])
      min.dist <- min(pdist[spp,other.species])
      a=y[j]
    }
    names(a) <- NULL
    return(c('abundance'=a,'distance'=min.dist))
  }
  
  getNTDvecs <- function(spp,OTU,ix,pdist.=pdist){
    y <- OTU[spp,]
    #limit community-community comparisons to where species is present in one, absent in another
    signed_abunds <- sign(as.numeric(y>0)-0.5)
    names(signed_abunds) <- names(y)
    contributing_pairs <- which(signed_abunds[ix$i]*signed_abunds[ix$j]<0)
    ix <- ix[contributing_pairs,]
    
    output <- matrix(NA,nrow=nrow(ix),ncol=2)
    for (ii in 1:nrow(ix)){
      output[ii,]<-unlist(NTD_calc(ix$i[ii],ix$j[ii],spp,y,pdist,OTU))
    }
    colnames(output) <- c('abundances','distances')
    output <- as.data.table(output)
    output <- cbind(output,ix)
    output[,pres_in_i:=y[i]>0]
    return(output)
  }
  
 
  cl <- makeCluster(ncores)
  clusterEvalQ(cl,library(data.table))
  clusterExport(cl,varlist=c('pdist','getNTDvecs','NTD_calc','ix'), envir= environment())
  timestamp()
  print("Making initial species table")
  species_NTD_effects <- parLapply(cl,species,getNTDvecs,OTU=OTU,ix=ix)
  stopCluster(cl)
  rm('cl')
  
  names(species_NTD_effects) <- species
  NTD_data <- species_NTD_effects %>%  rbindlist(use.names = T,idcol = T)
  names(NTD_data)[1] <- 'species'
  
  
  NTD_data[,community_index:=paste(i,j,sep='_')]
  NTD_data[,community_N:=.N,by=community_index]
  setkey(NTD_data,species,community_index)
  
  rDist <- NTD_data[,c('species','i','j')]
  rDist[,lagdist:=0]
  rDist[,dist:=0]
  rDist[,ldist:=0]
  rDist[,sq_dist:=0]
  rDist[,sq_ldist:=0]
  rDist[,sd_dist:=0]
  rDist[,sd_ldist:=0]
  
  mindist <- function(spp,i,j,dist=NULL,OTU.=OTU){
    y=OTU[spp,]
    if (y[i]>0 & y[j]==0){
      other.species <- setdiff(colnames(pdist),spp)
      other.species <- intersect(other.species,rownames(OTU)[which(OTU[,j]>0)])
      min.dist <- min(dist[spp,other.species])
      
    } else if (y[j]>0 & y[i]==0){
      other.species <- setdiff(colnames(pdist),spp)
      other.species <- intersect(other.species,rownames(OTU)[which(OTU[,i]>0)])
      min.dist <- min(dist[spp,other.species])
      
    }
    return(min.dist)
  }
  
  rdist_calc <- function(spp,i,j,OTU.=OTU,pdist.=pdist) sapply(spp,mindist,i[1],j[1],taxaShuffle(pdist))
  
  rdists <- function(nreps,rDist.=rDist,OTU.=OTU,pdist.=pdist){
    for (z in 1:nreps){
      rDist[,lagdist:=dist]
      rDist[,dist:=dist+rdist_calc(species,i,j),by=c('i','j')]
      rDist[,ldist:=ldist+log(dist-lagdist)]
      rDist[,sq_dist:=sq_dist+(dist-lagdist)^2/nreps]
      rDist[,sq_ldist:=sq_ldist+log(dist-lagdist)^2/nreps]
    }
    rDist[,dist:=dist/nreps]
    rDist[,ldist:=ldist/nreps]
    rDist[,sd_dist:=sqrt(sq_dist-dist^2)]
    rDist[,sd_ldist:=sqrt(sq_ldist-ldist^2)]
    rDist$sq_dist <- NULL
    rDist$sq_ldist <- NULL
    return(rDist)
  }
  
  cl <- makeCluster(ncores)
  clusterEvalQ(cl,library(picante))
  clusterEvalQ(cl,library(data.table))
  clusterExport(cl,varlist = c('rDist','OTU','pdist','rdist_calc','mindist'), envir= environment())
  
  reps <- rep(ceiling(nreps/ncores),ncores)
  
  timestamp()
  print("Forming null distributions - this may take a few days so chill a bit!")
  rDists <- parLapply(cl,reps,rdists)
  stopCluster(cl)
  rm('cl')
  timestamp()
  
  for (i in 1:length(rDists)){
    rDist$dist <- rDist$dist+rDists[[i]]$dist/ncores
    rDist$ldist <- rDist$ldist+rDists[[i]]$ldist/ncores
    rDist$sd_dist <- rDist$sd_dist+rDists[[i]]$sd_dist/ncores
    rDist$sd_ldist <- rDist$sd_ldist+rDists[[i]]$sd_ldist/ncores
  }
  
  rDist$lagdist <- NULL
  rDist$sq_dist <- NULL
  rDist$sq_ldist <- NULL
  
  names(rDist)[4] <- 'null_dist'
  names(rDist)[5] <- 'null_ldist'
  names(rDist)[6] <- 'null_sd_dist'
  names(rDist)[7] <- 'null_sd_ldist'
  
  setkey(rDist,species,i,j)
  setkey(NTD_data,species,i,j)
  
  NTD_data <- NTD_data[rDist]
  NTD_data[,l_z_score:=(log(distances)-null_ldist)/null_sd_ldist]
  
  ntd_spp_cluster <- NTD_data[,list(avg_score=mean(l_z_score), avg_std_score=mean(l_z_score)*(.N), sum_score=sum(l_z_score), median_score=median(l_z_score), N_comm_pairs = sum(l_z_score)/mean(l_z_score)),
                                            by='species']
  
  setkey(taxonomy,species)
  setkey(ntd_spp_cluster,species)
  ntd_spp_cluster <- taxonomy[ntd_spp_cluster]
  
  #Subsetting species with median z-scores lower than -2
  ntd_spp_cluster_constant <- subset(ntd_spp_cluster, subset = ntd_spp_cluster$median_score < -2)
  
  #Clustering by unique taxonomies
  taxa_weights <- ntd_spp_cluster_constant[,list(sum_score=sum(sum_score), avg_sum_score=mean(sum_score), sum_comm_pairs=sum(N_comm_pairs), avg_comm_pairs=mean(N_comm_pairs)), by='taxonomy']
  taxa_weights[,N:=(sum_score/avg_sum_score)]
  
  # phylofactorization of sum_score ----------------------------------------
  
  print("Performing phylofactorization")
  sum_score_vec <- vector(length = Ntip(phy))
  for (i in 1:length(sum_score_vec)) {
    for (j in 1:length(sum_score_vec)) {
      if (phy$tip.label[i] == ntd_spp_cluster$species[j]) 
        sum_score_vec[i] <- ntd_spp_cluster$sum_score[j]
    }
  } #I hate myself for using for loops but goddamn them they are very intuitive
  
  pf_sum_score <- twoSampleFactor(sum_score_vec, tree = phy, nfactors = nfactors, ncores = ncores) #The user is prompted to check the distribution of sum_score and standardize, if needed, before running this
  pf_sum_score_tree <- pf.tree(pf_sum_score, tree= phy, layout='rectangular')
  pp <- pf_sum_score_tree$ggplot
  
  class(pf_sum_score) <- "list"
  phylo_z <- setClass("phylo_z", slots=list(species_cluster="data.table", constant="data.table", taxa.weights="data.table", phylofactorization="list", phylo_tree="list"))
  output <- phylo_z(species_cluster = ntd_spp_cluster, constant = ntd_spp_cluster_constant, phylofactorization = pf_sum_score, taxa.weights = taxa_weights, phylo_tree = pf_sum_score_tree)
  return(output)
  
}
