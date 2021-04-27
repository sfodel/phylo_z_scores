# phylo_z_scores
Null modelling of phylogenetic nearest taxon distances for individual taxa

This repository contains the R code to identify phylogenetic clades under homogeneous ecological selection (HoS clades), that is monophyletic clades of
taxa (e.g. Sequence Variants in a 16S rRNA gene survey) with phylogenetically closer relatives across communities than expected by chance.
The method consists of the following steps: 

1) For a given pair of communities, j, k, and for each SV, i, that is present in one but not both communities, we calculate its “phyloscore”. The phyloscore is a z-score quantifying how different its nearest phylogenetic distance is to a null expectation in which species are randomly drawn to be present in the community in which SV i is absent. For example, if we examine SV i across communities j and k and i is present in community j and not in community k, we first find the nearest phylogenetic distance di,j,k of i based on the SVs that are present in community k. We then sample a null distribution of M minimum phylogenetic distances {d_(i,j,k,m)^0 }_(m=1)^(m=M) between our focal SV i and the SVs present in community k in which SV i is absent. If there are Nk species present in community k, we randomly sample Nk species other than SV i to be present in the null community, compute the distance to the nearest present taxon to our focal SV, and repeat this process M=100 times to estimate the distribution of nearest phylogenetic distances in our null model. Finally, we calculate the phyloscore as:

 z_(i,j,k)=(log⁡(d_(i,j,k) )-〈log⁡〖d_(i,j,k,m)^0 〗 〉_m)/(σ_(i,j,k)^0 )

where 〈log⁡〖d_(i,j,k,m)^0 〗 〉_m  is the average of the null distribution’s log-transformed nearest taxon distances and σ_(i,j,k)^0 the standard deviation of this distribution.

2) We then calculate for each SV its total phyloscore as the sum of its phyloscores across all community pairs. We use phylofactorization (https://github.com/reptalex/phylofactor, https://esajournals.onlinelibrary.wiley.com/doi/abs/10.1002/ecm.1353)  to identify monophyletic clades of SVs with significantly different total phyloscores compared to the complement set of SVs and to extract the consensus taxonomic classification of the SVs within. Phylofactorization is a graph-partitioning algorithm that sequentially cuts edges in a phylogeny, splitting the tree into disjoint sub-trees with high within-group similarities and between-group differences. At each iteration, phylofactorization cuts the edge that maximizes an objective function quantifying the difference between the sets of SVs on either side of the edge. Here, the objective function was the absolute value of the t-statistic from a two-sample t-test of equal variance on the total phyloscores between the two groups of SVs on each side of the edge. The output of phylofactorization will be phylogenetic clades containing SVs with distinctly different total phyloscores compared to outgroups. This second step is particularly important for distinguishing between niche-related patterns and dispersal. Dispersal would result in uniform phyloscores across all the phylogeny -high scores in the case of dispersal limitation and low scores in the case of homogenizing dispersal- whereas clade-specific high niche occupancy would result in monophyletic groups with lower/higher phyloscores, respectively, compared to outgroups. 
Importantly, because the phylogenetic distance pool is preserved across all permutations, the phyloscores for each SV are determined by presence-absence patterns alone and are independent of the branch lengths and patterns of speciation in the phylogeny. Because of that, we need to ensure that all potential sources of bias to the presence/absence matrix are excluded prior to the calculation of the phyloscores. To our perception the potential sources of bias can be either sequencing error biases or inadequate sampling biases. We describe below how we treated both these potential biases, and we recommend similar assessment in studies using our framework.
a) Sequencing error biases. Sequencing errors can skew the distribution of presences/absences by inputting false positives, i.e., non-existent presences. Because these errors are more likely to happen in more abundant sequences, these false positives might tend to cluster around abundant SVs in the dataset, artificially decreasing the phyloscores of abundant phylogenetic clades. To treat these potential biases we excluded SVs observed in only one replicate sample. Taking into consideration the error correction implemented the usual microbiome analyses tools such as in DADA2, we have no reason to believe that any residual errors will be differentially distributed between HoS and non-HoS clades. 
b) Inadequate sampling biases. Inversely to sequencing errors, inadequate sampling biases can introduce false negatives (non-existent absences). In other words, low sequencing effort can be enough to capture all the diversity of abundant phylogenetic clades but not that of less abundant ones. In this way some SVs in the latter clades can be left out and the presence/absence matrix can be artificially sparse in these clades. This concerns not only each HoS clade, but more importantly the outgroups against which these clades have distinctly different phyloscores. Thus, we needed to ensure that the non-HoS clades and each HoS clade are adequately sampled. For that we recommend to perform individual rarefactions for each of these phylogenetic groups and to ensure that all such curves saturate, supporting that there is no bias due to inadequate sampling effort for any of the clades in question.


