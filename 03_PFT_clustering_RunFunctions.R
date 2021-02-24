#-------------------------------------------------------------------------------
#  Load the data
# plant functional traits, cPFTs and (brightness-normalized) reflectance
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Data and analyses/Raw data")
data_PFT_spec <- Load_PFT_spec(LUTdata = "_refl_sg", width = 1 )
data_PFT_spec_br <- Load_PFT_spec(LUTdata = "_refl_br_sg", width = 1 )

# plant optical traits
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Data and analyses/PLSR results/b) PLSR best models")
data_pred_plsr_cv <- read.csv("PLSR_best_predictions_cv.csv",sep=";",header=TRUE, row.names=1)

#--------------------------------------------------------------------------------
# Choose dataset and traits to work with

setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Data distribution")
Traits <-c("SLA","LDMC","LNC_mass","Chltotal_mass","LA","Height") # include lifespan?
traits_field <- data_to_cluster_trans(data_PFT_spec$PFT_traits, data_PFT_spec$PFT_meta, Traits, name="traits_field",
                                      trait_log=c("LA"), trait_sqrt = c())

Traits_plsr <- c("SLA_log","LDMC_log","LNC_mass_log","Chltotal_mass_sqrt","LA_log","Height_sqrt")
traits_plsr_cv <- data_to_cluster_backwardtrans(data_pred_plsr_cv, data_PFT_spec$PFT_meta, Traits_plsr, name = "PFT_PLSR_cv", 
                                                trait_log=c("LNC_mass_log","SLA_log","LDMC_log"), trait_sqrt=c("Chltotal_mass_sqrt","Height_sqrt"))
# Keep optical LA log-transformed

colnames(traits_field$datavalues_trans) <- c("SLA","LDMC","Leaf N", "Leaf Chl a+b", "log(Leaf area)","Plant Height")
colnames(traits_plsr_cv$datavalues_backtrans) <- c("SLA","LDMC","Leaf N", "Leaf Chl a+b", "log(Leaf area)","Plant Height")
#--------------------------------------------------------------------------------
# ORDINATION of spectral data

# 1. PCA
# setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/PCA reflectance")
# 
# pca_spec <- pca_explore_spec(data_PFT_spec)
# pca_spec_SF <- pca_explore_spec_SuppFig(data_PFT_spec, traits_field$datavalues, name="Refl")
# 
# pca_spec_br <- pca_explore_spec(data_PFT_spec_br)
# pca_spec_br_SF <- pca_explore_spec_SuppFig(data_PFT_spec_br, traits_field$datavalues, name="Refl_br")
# 
# graphics.off()

#--------------------------------------------------------------------------------
# CORRELATION STRUCTURE

setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Correlation scatterogram")
trait_names <- c("SLA", "LDMC", "Leaf N", "Leaf Chl a+b", "log(Leaf area)", "Plant height")
corr_traits(traits_field, transf = "_trans", std ="", name="field",trait_names)
corr_traits(traits_plsr_cv, transf = "_backtrans", std ="", name="plsr_cv",trait_names)

#--------------------------------------------------------------------------------
# CLUSTERING

# ------- Cluster algorithm: hierarchical and kmeans
# Clustering on functional traits
cluster_PFT <- cluster_traits(traits_field$datavalues_trans, "euc", "Ward", 10,stand=TRUE) 

# Clustering on functional traits predicted from reflectance through PLS regression
cluster_POT_cv_backtrans <- cluster_traits(traits_plsr_cv$datavalues_backtrans, "euc", "Ward", 10,stand=TRUE)

# Clustering on PC's
  # Set number of variables (=PC's) equal to number of traits used above
  # cluster_POT_refl <- cluster_traits(pca_spec$x[,1:6], "euc", "Ward", 10,stand=TRUE) 
  # cluster_POT_refl_br <- cluster_traits(pca_spec_br$x[,1:6], "euc", "Ward", 10,stand=TRUE) 

graphics.off()
    
# Save an overview of cluster assignment per observation
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Overview observations - groups")
overview <- c("cluster_PFT",
              "cluster_POT_cv_backtrans")
              # "cluster_POT_refl","cluster_POT_refl_br")

for (i in 1:length(overview)){
  table <- cbind("Observation ID" = labels(get(overview[i])$distancematrix),
                 "group no. DHC" = get(overview[i])$groups_DHC,
                 "group no. NbClust" = get(overview[i])$groups_NbClust,
                 "group no. Kmeans" = get(overview[i])$groups_Kmeans,
                 "group no. Kmeans_alg" = get(overview[i])$groups_Kmeans_alg)
  write.table(table,paste(overview[i],"csv",sep="."), sep=";", col.names = T, row.names=F)
}
# !!! the results of Kmeans and Kmeans_alg are not the same!
# Probably because the algorithm each time produces slightly different results
    
# ------- Visualise clusters on pca
  # Note from FLAMES course:
  #' Both 'prcomp()' and 'princomp()' give the square roots of the eigenvalues ("standard deviations").
  #' Also, both functions do not give the loadings (= matrix A) but the eigenvectors (= matrix U):
  #' sd >< eigenvalues --> square them
  #' Eigenvectors (= matrix U): pca$rotation
  #' Component matrix (= matrix A): pca$rotation %*% diag(pca$sdev)
  #' Or: tcrossprod(pca$rotation, diag(pca$sdev))

  setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Clusters on PCA")
  
  # Clusters on functional traits indicated on PCA of functional traits
  # for plotting purposes, turn axes if wanted
  
  pca_PFT <- prcomp(traits_field$datavalues_trans,scale=TRUE,center=TRUE)
  pca_PFT$loadings <- pca_PFT$rotation %*% diag(pca_PFT$sdev)
  turn=c(-1,-1,-1)
  axes <- c(1, 2, 3)
  for (i in 1:length(axes)){
    pca_PFT$loadings[,i] <- pca_PFT$loadings[,i] * turn[i]
    pca_PFT$x[,i] <- pca_PFT$x[,i] * turn[i] 
    pca_PFT$rotation[,i] <- pca_PFT$rotation[,i] * turn[i]
  }  
  
  summary(pca_PFT)
  clusters_on_pca(pca_PFT, traits_field, cluster_PFT, clusterdecision = "DHC", form =  "trans", name = "traits_field_", colorscheme="Dark2")
  clusters_on_pca_ellipses(pca_PFT, 0.68, group_order_nb=c(7,10,3,4,1,2,8,5,9,6), text="ePFT ",
                           traits_field, cluster_PFT, clusterdecision = "DHC", form =  "trans", name = "traits_field_", colorscheme="Dark2")
  clusters_on_pca_ellipses(pca_PFT, 0.95, group_order_nb=c(7,10,3,4,1,2,8,5,9,6), text="ePFT ",
                           traits_field, cluster_PFT, clusterdecision = "DHC", form =  "trans", name = "traits_field_", colorscheme="Dark2")
  graphics.off()
  
# Clusters on PLSR predicted functional traits indicated on PCA of PLSR predicted functional traits
  pca_POT_PLSR_cv_backtrans <- prcomp(traits_plsr_cv$datavalues_backtrans,scale=TRUE,center=TRUE)
  pca_POT_PLSR_cv_backtrans$loadings <- pca_POT_PLSR_cv_backtrans$rotation %*% diag(pca_POT_PLSR_cv_backtrans$sdev)
  
  turn=c(-1,-1,1)
  axes <- c(1, 2, 3)
  for (i in 1:length(axes)){
    pca_POT_PLSR_cv_backtrans$loadings[,i] <- pca_POT_PLSR_cv_backtrans$loadings[,i] * turn[i]
    pca_POT_PLSR_cv_backtrans$x[,i] <- pca_POT_PLSR_cv_backtrans$x[,i] * turn[i] 
    pca_POT_PLSR_cv_backtrans$rotation[,i] <- pca_POT_PLSR_cv_backtrans$rotation[,i] * turn[i]
  }  
  
  summary(pca_POT_PLSR_cv_backtrans)
  clusters_on_pca(pca_POT_PLSR_cv_backtrans,traits_plsr_cv, cluster_POT_cv_backtrans, clusterdecision = "DHC", form =  "backtrans", name = "traits_PLSR_cv_backtrans_", colorscheme = "Set1")
  clusters_on_pca_ellipses(pca_POT_PLSR_cv_backtrans, 0.68, group_order_nb=c(6,9,3,2,4,8,5,7,1), text="ePOT ",
                           traits_plsr_cv, cluster_POT_cv_backtrans, clusterdecision = "DHC", form =  "backtrans", name = "traits_PLSR_cv_backtrans_", colorscheme="Set1")
  clusters_on_pca_ellipses(pca_POT_PLSR_cv_backtrans, 0.95, group_order_nb=c(6,9,3,2,4,8,5,7,1), text="ePOT ",
                           traits_plsr_cv, cluster_POT_cv_backtrans, clusterdecision = "DHC", form =  "backtrans", name = "traits_PLSR_cv_backtrans_", colorscheme="Set1")
  graphics.off()
  
  # Clusters on refl indicated on PCA of PLSR predicted functional traits
  # clusters_on_pca(pca_POT_PLSR_cv_backtrans,traits_plsr_cv, cluster_POT_refl, clusterdecision = "DHC", form =  "backtrans", name = "refl_traits_PLSR_cv_backtrans_", colorscheme = "Set1")
  # clusters_on_pca(pca_POT_PLSR_cv_backtrans,traits_plsr_cv, cluster_POT_refl_br, clusterdecision = "DHC", form =  "backtrans", name = "refl_br_traits_PLSR_cv_backtrans_", colorscheme = "Set1")
  # graphics.off()
  
 
# ------- Visualise clusters as dendrogram
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Dendrogram")
# Dendrogram of clustering on functional traits
dend_PFT_DHC <- plot_dend("ward.D2", cluster_PFT, clusterdecision = "DHC", name = "dend_PFT_", colorscheme = "Dark2") 

# Dendrogram of clustering on PLSR predicted functional traits
dend_POT_PLSR_cv_backtrans_DHC <- plot_dend("ward.D2", cluster_POT_cv_backtrans, clusterdecision = "DHC", name = "dend_POT_PLSR_cv_backtrans_", colorscheme = "Set1") 
graphics.off()

# Dendrogram of clustering on reflectance
# dend_POT_refl_DHC <- plot_dend("ward.D2", cluster_POT_refl, clusterdecision = "DHC", name = "dend_POT_refl_", colorscheme = "Set1") 
# dend_POT_refl_br_DHC <- plot_dend("ward.D2", cluster_POT_refl_br, clusterdecision = "DHC", name = "dend_POT_refl_br_", colorscheme = "Set1") 
# graphics.off()

# ------- Visualise 2 dendrograms as tanglegram
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Tanglegram")
dend_compare_tanglegram(cluster_PFT, dend_PFT_DHC, dend_POT_PLSR_cv_backtrans_DHC, clusterdecision = "DHC", name = "DHC_PFT_vs_PLSR_cv", colorsheme1="Dark2")
# dend_compare_tanglegram(cluster_PFT, dend_PFT_DHC, dend_POT_refl_DHC, clusterdecision = "DHC", name = "DHC_PFT_vs_refl", colorsheme1="Dark2")
# dend_compare_tanglegram(cluster_PFT, dend_PFT_DHC, dend_POT_refl_br_DHC, clusterdecision = "DHC", name = "DHC_PFT_vs_refl_br", colorsheme1="Dark2")


#--------------------------------------------------------------------------------
# CLUSTERING EVALUATION
  
##### --------- 1. Cluster validity and stability:
cluster_results <- list(cluster_PFT = cluster_PFT,
                        cluster_POT_cv_backtrans = cluster_POT_cv_backtrans)
                        # cluster_POT_refl = cluster_POT_refl,
                        # cluster_POT_refl_br = cluster_POT_refl_br) 
cluster_quality <- matrix(nrow=length(cluster_results), ncol=3)
colnames(cluster_quality) <- c("dataset","cophcorr", "agglcoef")
for (i in 1:length(cluster_results)){
  cluster_quality[i,"dataset"] <- names(cluster_results)[i]
  cluster_quality[i,"cophcorr"] <- round(ClusterQuality(cluster_results[[i]], "DHC", "test")$cophcorr,2)
  cluster_quality[i,"agglcoef"] <- round(ClusterQuality(cluster_results[[i]], "DHC", "test")$agglcoef,2)
}
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Clustering quality")
write.table(cluster_quality, file="Cluster_quality_DHC.csv", sep =";", row.names=F, col.names=T) 
graphics.off()


# Bootstrap Jaccard similarity value <= 0.5: indication of a "dissolved cluster".
# Bootstrap Jaccard similarity value => 0.85: "Highly stable" clusters 
# Silhoutte A value close to 1 indicates that the data point is assigned to a very appropriate cluster.
ClusterStability(cluster_PFT, clusterdecision="DHC")
ClusterStability(cluster_POT_cv_backtrans, clusterdecision="DHC")
graphics.off()
# Stability of clusters is really bad!


##### --------- 2. Describe clusters: species membership
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Species group membership")
membership_PFT <- membership(traits_field, cluster_PFT, "DHC", name="membership_PFT_")
membership_POT_cv_backtrans <- membership(traits_plsr_cv, cluster_POT_cv_backtrans, "DHC", name="membership_POT_cv_backtrans_")
graphics.off()



##### 3. Describe clusters: Standardized trait profile of each group (after log-transformation if necessary)
#        = Fig. 10: from Roth et al. (2016)
# Fig_ProfileGroups(cluster_traits_field, "DHC", traits_field$datavalues_log, name = "ProfileGroups_traits_field_", colorscheme = "Dark2")
# Fig_ProfileGroups(cluster_PFT_PLSR_cal_mass, "DHC", traits_field_plsr_cal$datavalues_log, name = "ProfileGroups_PFT_PLSR_cal_byPFT_PLSR_cal_mass_", colorscheme = "Set1")
# Fig_ProfileGroups(cluster_PFT_PLSR_cv_mass, "DHC", traits_field_plsr_cv$datavalues_log, name = "ProfileGroups_PFT_PLSR_cv_byPFT_PLSR_cv_mass_", colorscheme = "Set1")


# ------------------ Group profiles ordered according to appearance in tanglegram
# group_order_nb specifies the order number of each group starting from 1 (define after plotting tanglegram)
# in other words: from group 1 until last group, indicate it's position/order in the dendrogram 
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Profile groups")

trait_labels <-  c("SLA","LDMC","Leaf N", "Leaf Chl a+b", "log(Leaf area)","Plant Height")
Fig_ProfileGroups_ordered(cluster_PFT, "DHC", traits_field$datavalues_trans, name = "ProfileGroups_ordered_PFT_", colorscheme = "Dark2",
                          group_order_nb=c(7,10,3,4,1,2,8,5,9,6), text="ePFT ",  trait_labels)
Fig_ProfileGroups_ordered(cluster_POT_cv_backtrans, "DHC", traits_plsr_cv$datavalues_backtrans, name = "ProfileGroups_ordered_POT_PLSR_cv_backtrans_", colorscheme = "Set1", 
                          group_order_nb=c(6,9,3,2,4,8,5,7,1), text="ePOT ", trait_labels)

# Profiles of POT described by field-measured traits
Fig_ProfileGroups_ordered(cluster_POT_cv_backtrans, "DHC", traits_field$datavalues_trans, name = "ProfileGroups_ordered_POT_PLSR_cv_backtrans_by_traits_field", colorscheme = "Set1", 
                          group_order_nb=c(6,9,3,2,4,8,5,7,1), text="ePOT ", trait_labels)
graphics.off()


##### 4. Describe clusters: Spectral profile of each group
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Profile groups reflectance")
Fig_ProfileGroupsRefl_ordered(cluster_PFT, "DHC", data_PFT_spec, name = "ProfileGroupsRefl_PFT_", colorscheme = "Dark2",  
                              group_order_nb=c(7,10,3,4,1,2,8,5,9,6))
Fig_ProfileGroupsRefl_ordered(cluster_POT_cv_backtrans, "DHC", data_PFT_spec, name = "ProfileGroupsRefl_POT_cv_backtrans_", colorscheme = "Set1",  
                              group_order_nb=c(6,9,3,2,4,8,5,7,1))


##### --------- 5. Describe clusters: Distribution of groups per trait 
# Kruskal Wallis test: does not assume normality, nor homoscedasticity
wd1 <- "C:/Users/u0091812/Box Sync/04. Optical types/Figures/KruskalWallis"
wd2 <- "C:/Users/u0091812/Box Sync/04. Optical types/Figures/Profile traits_KW exploration"
KW_PFT<-KruskalWallis(traits_field$datavalues_trans, cluster_PFT$groups_DHC,wd1, wd2, name="PFT")
KW_POT_cv_backtrans<-KruskalWallis(traits_plsr_cv$datavalues_backtrans, cluster_POT_cv_backtrans$groups_DHC,wd1, wd2, name="POT_cv_backtrans")
dev.off()

# Post hoc Dunn test
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Profile traits_boxplot_Dunn")
facet_names <-  c("SLA","LDMC","Leaf N", "Leaf Chl a+b", "log(Leaf area)","Plant Height")
Fig_ProfileTraits_Dunn(cluster_PFT, "DHC", traits_field$datavalues_trans, "PFT", "Dark2",
                       group_order_nb=c(7,10,3,4,1,2,8,5,9,6), text="ePFT ", facet_names, trans="")
Fig_ProfileTraits_Dunn(cluster_POT_cv_backtrans, "DHC", traits_plsr_cv$datavalues_backtrans, "POT_cv_backtrans", "Set1",
                      group_order_nb=c(6,9,3,2,4,8,5,7,1), text="ePOT ", facet_names, trans="")
graphics.off()


##### 6. Describe clusters: PERMANOVA: amount of functional variance explained by groups
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Permanovatest")
# Multivariate variation between functional traits explained
traits_structure <- c("Plant Height","log(Leaf area)")
traits_economics<-c("Leaf N","SLA","LDMC","Leaf Chl a+b")
Permanova_PFT_traits_field <- permanova(traits_field$datavalues_trans, cluster_PFT$groups_DHC, 
                                        partitioning = "PFT_traits_field",999, traits_structure, traits_economics, trans="")
Permanova_POT_cv_backtrans_traits_field <- permanova(traits_field$datavalues_trans, cluster_POT_cv_backtrans$groups_DHC,
                                           "POT_cv_backtrans_traits_field",999, traits_structure, traits_economics, trans="")
# Permanova_POT_refl_traits_field <- permanova(traits_field$datavalues_trans, cluster_POT_refl$groups_DHC, 
#                                              "POT_refl_traits_field",999, traits_structure, traits_economics, trans="")
# Permanova_POT_refl_br_traits_field <- permanova(traits_field$datavalues_trans, cluster_POT_refl_br$groups_DHC, 
#                                                 "POT_refl_br_traits_field",999, traits_structure, traits_economics, trans="")



# Multivariate variation between optical traits explained
traits_structure_PLSR <- c("Plant Height","log(Leaf area)")
traits_economics_PLSR<-c("Leaf N","SLA","LDMC","Leaf Chl a+b")
Permanova_POT_cv_backtrans_traits_plsr_cv_backtrans <- permanova(traits_plsr_cv$datavalues_backtrans, cluster_POT_cv_backtrans$groups_DHC, 
                                                                 "POT_cv_backtrans_traits_plsr_cv_backtrans",999, traits_structure_PLSR, traits_economics_PLSR, trans="")
# Permanova_POT_refl_traits_plsr_cv_backtrans <- permanova(traits_plsr_cv$datavalues_backtrans, cluster_POT_refl$groups_DHC, 
#                                                "POT_refl_traits_plsr_cv_backtrans",999, traits_structure_PLSR, traits_economics_PLSR, trans="")
# Permanova_POT_refl_br_traits_plsr_cv_backtrans <- permanova(traits_plsr_cv$datavalues_backtrans, cluster_POT_refl_br$groups_DHC, 
#                                                          "POT_refl_br_traits_plsr_cv_backtrans",999, traits_structure_PLSR, traits_economics_PLSR, trans="")
Permanova_PFT_traits_plsr_cv_backtrans <- permanova(traits_plsr_cv$datavalues_backtrans, cluster_PFT$groups_DHC,
                                                            "PFT_traits_plsr_cv_backtrans",999, traits_structure_PLSR, traits_economics_PLSR, trans="")

Permanova_traits_field <- rbind(Permanova_PFT_traits_field, Permanova_POT_cv_backtrans_traits_field)
Permanova_traits_plrs_cv_backtrans <- rbind(Permanova_PFT_traits_plsr_cv_backtrans, Permanova_POT_cv_backtrans_traits_plsr_cv_backtrans)
# Permanova_traits_field <- rbind(Permanova_PFT_traits_field, Permanova_POT_cv_backtrans_traits_field, Permanova_POT_cal_backtrans_traits_field, Permanova_POT_refl_traits_field, Permanova_POT_refl_br_traits_field)
# Permanova_traits_plrs_cv_backtrans <- rbind(Permanova_PFT_traits_plsr_cv_backtrans, Permanova_POT_cv_backtrans_traits_plsr_cv_backtrans, Permanova_POT_refl_traits_plsr_cv_backtrans, Permanova_POT_refl_br_traits_plsr_cv_backtrans)

conv_groups <- c("Plant_growth_form_LEDA", "Plant_lifespan_LEDA", "Plant_growth_habit", "CSR_Hunt")
for (i in 1 : length(conv_groups)){
  # Multivariate variation explained of functional traits 
  Permanova_cPFT_traits_field <- permanova(traits_field$datavalues_trans, traits_field$meta[,conv_groups[i]], 
                                           paste(conv_groups[i],"traits_field",sep="_"),999, traits_structure, traits_economics, trans="")
   Permanova_traits_field <- rbind(Permanova_traits_field, Permanova_cPFT_traits_field)
  # Permanova_traits_field_notrans <- rbind(Permanova_traits_field_notrans, Permanova_cPFT_traits_field_notrans)
  
  # Multivariate variation explained of functional traits 
   Permanova_cPFT_traits_plsr_cv_backtrans <- permanova(traits_plsr_cv$datavalues_backtrans, traits_field$meta[,conv_groups[i]],
                                                       paste(conv_groups[i],"traits_cv_backtrans",sep="_"),999, traits_structure_PLSR, traits_economics_PLSR, trans="")
  Permanova_traits_plrs_cv_backtrans <- rbind(Permanova_traits_plrs_cv_backtrans, Permanova_cPFT_traits_plsr_cv_backtrans)
}

write.table(Permanova_traits_field, file="Permanova_traits_field.csv", sep =";", col.names=T, row.names = F)
write.table(Permanova_traits_plrs_cv_backtrans, file="Permanova_traits_plrs_cv_backtrans.csv", sep =";", col.names=T, row.names = F)



##### 7. Describe clusters: PERMANOVA: amount of functional variance explained by varying amount and choice of traits
# only applicable to partitionings that were not made by using Functional traits directly
# setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Permanova_loops")
# permanova_loop(traits_field$datavalues_trans, traits_field$meta[,"CSR_Hunt"], "CSR_Hunt",99)


#--------------------------------------------------------------------------------
# CLUSTERING CORRESPONDENCE
# Apply the functions

##### 1. Tanglegram and entanglement of two dendrograms
# see up


##### 2. Cophenetic correlationa and Baker's Gamma Index between dendrograms
dend_compare_metrics(dend_PFT_DHC,dend_POT_PLSR_cv_backtrans_DHC) # dend_POT_refl_DHC, dend_POT_refl_br_DHC)
graphics.off()



##### 3. Adjusted Rand Index and Fowlkes-Mallows Index between clusters
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Cluster correspondence")
corsp_PFT <- cluster_compare_metrics(traits_field, cluster_PFT, "DHC", 
                                     cluster_PFT, "DHC", name = "PFT")
corsp_PFT_POT_cv <- cluster_compare_metrics(traits_field, cluster_POT_cv_backtrans, "DHC",
                                                    cluster_PFT, "DHC", name = "PFT_POT_cv")
# corsp_PFT_POT_relf <- cluster_compare_metrics(traits_field, cluster_POT_refl, "DHC",
#                                               cluster_PFT, "DHC", name = "PFT_POT_refl")
# corsp_PFT_POT_relf_br <- cluster_compare_metrics(traits_field, cluster_POT_refl_br, "DHC",
#                                                  cluster_PFT, "DHC", name = "PFT_POT_refl_br")
corsp_PFT_POT <- cbind(corsp_PFT, corsp_PFT_POT_cv)
colnames(corsp_PFT_POT) <- paste(rep(c("PFT","POT_cv"),each=2),
                                 rep(c("ARI","FMI"),2))
# corsp_PFT_POT <- cbind(corsp_PFT, corsp_PFT_POT_cv, corsp_PFT_POT_cal, corsp_PFT_POT_relf, corsp_PFT_POT_relf_br)
# colnames(corsp_PFT_POT) <- paste(rep(c("PFT","POT_cv","POT_cal","POT_refl","POT_refl_br"),each=2),
#                                  rep(c("ARI","FMI"),5))
write.table(corsp_PFT_POT, file=paste("Cluster_Correspondence_PFT_POTs_overview","csv",sep="."), sep =";", row.names=T, col.names=NA)


##### 4. Procrustes analysis on 2 pca's
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Procrustes")
VI_pca_procrustes_ordered(pca_PFT, pca_POT_PLSR_cv_backtrans, nb_comp = 3, 
                          cluster_PFT, cluster_POT_cv_backtrans,
                          name="PFT_POT_cv_3PCs", colorscheme1="Dark2",
                          group1_order_nb=c(7,10,3,4,1,2,8,5,9,6), group2_order_nb=c(6,9,3,2,4,8,5,7,1))
VI_pca_procrustes_ordered(pca_PFT, pca_POT_PLSR_cv_backtrans, nb_comp = 6, 
                          cluster_PFT, cluster_POT_cv_backtrans,
                          name="PFT_POT_cv_6PCs", colorscheme1="Dark2",
                          group1_order_nb=c(7,10,3,4,1,2,8,5,9,6), group2_order_nb=c(6,9,3,2,4,8,5,7,1))
