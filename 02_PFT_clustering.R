######################################################################
###                                                                ###
###               Cluster functional and optical data              ###
###                                                                ###
###             Statistical and visual descrption of clusters      ###
###                       Comparisons with cPFTs                   ###
###                                                                ###
######################################################################


# In this script we code all the functions that are necessary to perform the analyses and graphic displays of the manuscript titled
# Spectrally defined plant functional types adequately capture multidimensional trait variation in herbaceous communities

# Remark: preparatory analyses that estimate optical traits from reflectance are covered in the script "01_PLSR trait estimation.R"

# Script by Elisa Van Cleemput, 2017-20

#--------------------------------------------------------------------------------

# Clean workspace

rm(list=ls())

#--------------------------------------------------------------------------------

library(xlsx)
# library(devtools)
# install_github("vqv/ggbiplot")
# library(ggbiplot) # installation error. Needed?
library(cluster)
library(gclus) # dmat.color
library(dynamicTreeCut) # cutreeHybrid
library(resemble) # SAM calculation with fDiss function, SID with sid function

setwd("C:/Users/u0091812/Box Sync/Literature/R packages/")
source("biostats.R") # Download this library from https://www.umass.edu/landeco/teaching/ecodata/labs/ecodata_labs.html 

#-------------------------------------------------------------------------------
#  Load the data
library(qdapTools) # lookup
Load_PFT_spec <- function(LUTdata, width) {
  ### 1. Load functional data
  PFT <- read.csv("Species_x_traits.csv",sep=";",header=TRUE)
    target <- which(names(PFT) == 'ID')[1]
    PFT_meta<-PFT[,1:target]
    PFT_traits<-PFT[,-c(1:target)]
    PFT_traits$Chltotal_mass <- PFT_traits$Chltotal
    PFT_traits$LNC_mass <- 10 * PFT_traits$LNC # % dry mass to mg/g
    PFT_traits$LNC_area <- PFT_traits$LNC_mass * 1/PFT_traits$SLA * 10^-1 # mg/cm²

  ### 2. Load data on CSR, taxonomy and growth form
  cPFT<-read.xlsx("Conventional_PFTs.xlsx", sheetName = "Overview for R")
    PFT_meta$CSR_Hunt <- lookup(PFT_meta$Species.code, cPFT[,c("Species.code","CSR_Hunt")])
    PFT_meta$CSR_Pierce <- lookup(PFT_meta$Species.code, cPFT[,c("Species.code","CSR_Pierce")])
    PFT_meta$Family <- lookup(PFT_meta$Species.code, cPFT[,c("Species.code","Family")])
    PFT_meta$Order <- lookup(PFT_meta$Species.code, cPFT[,c("Species.code","Order")])
    PFT_meta$Superorder <- lookup(PFT_meta$Species.code, cPFT[,c("Species.code","Superorder")])
    PFT_meta$Group <- lookup(PFT_meta$Species.code, cPFT[,c("Species.code","Group")])
    PFT_meta$Plant_growth_habit <- lookup(PFT_meta$Species.code, cPFT[,c("Species.code","Plant_growth_habit")])
    PFT_meta$Plant_lifespan_Pierce <- lookup(PFT_meta$Species.code, cPFT[,c("Species.code","Plant_lifespan_Pierce")])
    PFT_meta$Plant_lifespan_LEDA <- lookup(PFT_meta$Species.code, cPFT[,c("Species.code","Plant_lifespan_LEDA")])
    PFT_meta$Plant_growth_form_LEDA <- lookup(PFT_meta$Species.code, cPFT[,c("Species.code","Plant_growth_form_LEDA")])
    
  ### 4. Load spectral data
  spec <- read.csv("Species_x_reflectance_unmixed.csv",sep=";",header=TRUE)
    target <- which(names(spec) == 'nm350')[1]
    spectra<-spec[,target:ncol(spec)]
    meta<-spec[,1:target-1]
    rownames(meta)=meta[,"ID"]
  
  #--------------------------------------------------------------------------------
  # Preprocessing steps: spectral data
  # 1) smoothing
  # 2) Create library for 1st deriv or 2nd deriv spectra
  # 3a) spectral binning (optional) + removal of noise bands
  # 3b) removal of noise bands
  # 4) brightness normalization (Feilhauer et al. 2010)
  
  
  ### 1) Create speclib and Smooth
  library(hsdar)
  wl=as.numeric(gsub("nm","",names(spectra)))
  speclib <- speclib(as.matrix(spectra),wl)
  SI(speclib)<-meta
  idSpeclib(speclib)<-as.character(meta$ID)
  # speclib_smooth<- smoothSpeclib(speclib,method="sgolay", n=51) 
  
  library(stringr) #str_sub
  SI(speclib)$Instrument <- str_sub(SI(speclib)[,1],-3,-1)
  libr_sv<-subset(speclib,Instrument == ".sv")
  libr_asd<-subset(speclib,Instrument == "asd")
  libr_sv<- smoothSpeclib(libr_sv,method="sgolay", n=101)
  libr_asd<- smoothSpeclib(libr_asd,method="sgolay", n=51)
  spectra_sv <- spectra(libr_sv)
  meta_sv <- SI(libr_sv)
  spectra_asd <- spectra(libr_asd)
  meta_asd <- SI(libr_asd)
  spectra_total <- rbind(spectra_asd,spectra_sv)
  wl <- seq(350,2500,1)
  colnames(spectra_total)<-wl
  meta_total <- rbind(meta_asd,meta_sv)
  Refl <- cbind(meta_total, spectra_total)
  Refl <- Refl[order(Refl$ID),] # order the observations according to ID
  meta_ordered <- Refl[,c(1:which(names(Refl)=="Instrument"))]
  speclib_smooth <- speclib(as.matrix(Refl[,-c(1:which(names(Refl)=="Instrument"))]), wl)
  SI(speclib_smooth)<-meta_ordered
  
  
  ### 2) Create library for 1st deriv or 2nd deriv spectra
  if (isTRUE(grepl("1deriv",LUTdata,fixed=T))){
    speclib_smooth <- derivative.speclib(speclib_smooth, m = 1)
  } else if (isTRUE(grepl("2deriv",LUTdata,fixed=T))){
    speclib_smooth <- derivative.speclib(speclib_smooth, m = 2)
  } else {
    speclib_smooth <- speclib_smooth
  }
  
  ### 3a) Spectral binning (optional) + removal of noise bands
  if (width == "10"){
    # function for spectral binning
    wl <- c(seq(400, 1340, 10), seq(1460,1780,10), seq (1970, 2400, 10))
    resamp <- function (spec, wl, wl.out){
      out <- approx (wl, spec, xout=wl.out)$y
    }
    
    # apply function to spectra --> result = smoothed spectra in matrix: binned + chosen spectral regions only
    spectra_smooth <- speclib_smooth@spectra@spectra_ma # smoothed spectra in matrix
    wavelength <- speclib_smooth@wavelength # wavelengths in vector
    spectra <- t (apply (spectra_smooth, 1, resamp, wl=wavelength, wl.out=wl)) 
    mask(speclib_smooth)<-c(349,400,1340,1460,1780,1970,2400,2501) # for visualisation
    speclib_plot<-speclib(spectra,wl) #when visualising binned spectra
    SI(speclib_plot)<-meta_ordered
    
  }  else if (width == "1"){
  ### 3b) Removal of noise bands in case of no spectral binning
    mask(speclib_smooth)<-c(349,400,1340,1460,1780,1970,2400,2501) # for visualisation
    spectra <- speclib_smooth@spectra@spectra_ma # smoothed spectra in matrix
    wl<-speclib_smooth@wavelength
  }
  
  #### 4) Brightness normalization (optional) (only for raw and unmixed spectra, not 1deriv and 2deriv)
  if (isTRUE(grepl("br",LUTdata,fixed=T))){
    # Binned or detailed spectra
    brightness<-sqrt(apply(spectra^2,1,sum))
    spectra<-spectra/brightness
    speclib_smooth<-speclib(spectra,wl)
    SI(speclib_smooth)=meta_ordered
    mask(speclib_smooth)<-c(349,400,1340,1460,1780,1970,2400,2501) # for visualisation, this is binned in case of binning!
    
  }  else {
    spectra<-spectra
  }
  colnames(spectra) <- wl
  rownames(spectra) <- meta_ordered[,"ID"]
  
  output <- list ("PFT" = PFT,
                  "PFT_traits" = PFT_traits,
                  "PFT_meta" = PFT_meta,
                  "spectra" = spectra,
                  "spec_meta" = meta_ordered,
                  "wl" = wl)
  return(output)
}

#--------------------------------------------------------------------------------
# Choose dataset and traits to work with
data_to_cluster_trans <- function(dataset, meta, traitlist, name, trait_log, trait_sqrt){
  dataset_select <- dataset[,traitlist]
  rownames(dataset_select)<-meta[,'ID']
  label<-meta[,"Species.code"]
  
  # ----- Figures without transformations
  png(paste("Distr_dotchart_withoutTransf_",".png",sep=name), width = 12, height = 8, units = 'in', res = 600)
  # x11(width=2000,height=1200)
  par(mfrow=c(3,4))
  for (i in 1:length(traitlist)){
    dotchart(dataset_select[,i], main=colnames(dataset_select)[i])
  }
  mtext(paste(name,"data without transformations", sep=": "), outer = T, side = 3, line = -1.5, cex = 1)
  dev.off()
  
  png(paste("Distr_hist_withoutTransf_",".png",sep=name), width = 12, height = 8, units = 'in', res = 600)
  # x11(width=2000,height=1200)
  par(mfrow=c(3,4))
  for (i in 1:length(traitlist)){
    hist(dataset_select[,i], main=colnames(dataset_select)[i])
  }
  mtext(paste(name,"data without transformations", sep=": "), outer = T, side = 3, line = -1.5, cex = 1)
  dev.off()
  
  # ----- Transformations
  dataset_select_trans <- dataset_select
  # Log-transformations
  if (is.null(trait_log)){
    dataset_select_trans <- dataset_select_trans
  } else {
    for (i in 1:length(trait_log)){
      dataset_select_trans[,trait_log[i]] <- log10(dataset_select[,trait_log[i]])
    }
  }
  # Sqrt-transformations
  if (is.null(trait_sqrt)){
    dataset_select_trans <- dataset_select_trans
  } else {
    for (i in 1:length(trait_sqrt)){
      dataset_select_trans[,trait_sqrt[i]] <- sqrt(dataset_select[,trait_sqrt[i]])
    }
  }
  
  # ----- Figures with transformations
  png(paste("Distr_dotchart_withTransf_",".png",sep=name), width = 12, height = 8, units = 'in', res = 600)
  # x11(width=2000,height=1200)
  par(mfrow=c(3,4))
  for (i in 1:length(traitlist)){
    dotchart(dataset_select_trans[,i], main=colnames(dataset_select_trans)[i])
  }
  mtext(paste(name,"data with transformations", sep=": "), outer = T, side = 3, line = -1.5, cex = 1)
  dev.off()
  
  png(paste("Distr_hist_withTransf_",".png",sep=name), width = 12, height = 8, units = 'in', res = 600)
  # x11(width=2000,height=1200)
  par(mfrow=c(3,4))
  for (i in 1:length(traitlist)){
    hist(dataset_select_trans[,i], main=colnames(dataset_select_trans)[i])
  }
  mtext(paste(name,"data with transformations", sep=": "), outer = T, side = 3, line = -1.5, cex = 1)
  dev.off()
  
 
  list ("datavalues" = dataset_select,
        "datavalues_trans" = dataset_select_trans,
        "meta" = meta, 
        "label" = label)
}
data_to_cluster_backwardtrans <- function(dataset, meta, traitlist, name, trait_log, trait_sqrt){
  dataset_select <- dataset[,traitlist]
  rownames(dataset_select)<-meta[,'ID']
  label<-meta[,"Species.code"]
  
  # ----- Figures without transformations
  png(paste("Distr_dotchart_withoutTransf_",".png",sep=name), width = 12, height = 8, units = 'in', res = 600)
  # x11(width=2000,height=1200)
  par(mfrow=c(3,4))
  for (i in 1:length(traitlist)){
    dotchart(dataset_select[,i], main=colnames(dataset_select)[i])
  }
  mtext(paste(name,"data without transformations", sep=": "), outer = T, side = 3, line = -1.5, cex = 1)
  dev.off()
  
  png(paste("Distr_hist_withoutTransf_",".png",sep=name), width = 12, height = 8, units = 'in', res = 600)
  # x11(width=2000,height=1200)
  par(mfrow=c(3,4))
  for (i in 1:length(traitlist)){
    hist(dataset_select[,i], main=colnames(dataset_select)[i])
  }
  mtext(paste(name,"data without transformations", sep=": "), outer = T, side = 3, line = -1.5, cex = 1)
  dev.off()
  
  # ----- Transformations
  dataset_select_backtrans <- dataset_select
  # Log-transformations
  if (is.null(trait_log)){
    dataset_select_backtrans <- dataset_select_backtrans
  } else {
    for (i in 1:length(trait_log)){
      dataset_select_backtrans[,trait_log[i]] <- 10^(dataset_select[,trait_log[i]])
    }
  }
  # Sqrt-transformations
  if (is.null(trait_sqrt)){
    dataset_select_backtrans <- dataset_select_backtrans
  } else {
    for (i in 1:length(trait_sqrt)){
      dataset_select_backtrans[,trait_sqrt[i]] <- (dataset_select[,trait_sqrt[i]])^2
    }
  }
  
  # ----- Figures with transformations
  png(paste("Distr_dotchart_withBackTransf_",".png",sep=name), width = 12, height = 8, units = 'in', res = 600)
  # x11(width=2000,height=1200)
  par(mfrow=c(3,4))
  for (i in 1:length(traitlist)){
    dotchart(dataset_select_backtrans[,i], main=colnames(dataset_select_backtrans)[i])
  }
  mtext(paste(name,"data with backward transformations", sep=": "), outer = T, side = 3, line = -1.5, cex = 1)
  dev.off()
  
  png(paste("Distr_hist_withBackTransf_",".png",sep=name), width = 12, height = 8, units = 'in', res = 600)
  # x11(width=2000,height=1200)
  par(mfrow=c(3,4))
  for (i in 1:length(traitlist)){
    hist(dataset_select_backtrans[,i], main=colnames(dataset_select_backtrans)[i])
  }
  mtext(paste(name,"data with backward transformations", sep=": "), outer = T, side = 3, line = -1.5, cex = 1)
  dev.off()
  
  
  list ("datavalues" = dataset_select,
        "datavalues_backtrans" = dataset_select_backtrans,
        "meta" = meta, 
        "label" = label)
}

#--------------------------------------------------------------------------------
# ORDINATION
#--------------------
# 1. PCA
library(randomcoloR) # distinctColorPalette
library(gginnards) # move_layers # https://cran.r-project.org/web/packages/gginnards/vignettes/user-guide-2.html

# myggbiplot (based on ggbiplot function) creates PCA figures
library(ggplot2)
library(plyr)
library(scales)
library(grid)
library(ggrepel) #geom_text_repel
library(randomcoloR) # distinctColorPalette
myggbiplot <- function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                        obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                        ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                        alpha = 1, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                        varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                        col_arrows="gray20",var.labels=FALSE, # <- add new arguments
                        ...) {
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  
  # Variable names
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  
  # Variables for text label placement
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  
  # Base plot
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    # g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
    #                                        xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
    #                                                                                               "picas")), color = muted("red"))
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0,
                                           xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
                                                                                                  "picas")), color = col_arrows)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.labels) {
    # g <- g + geom_text(data = df.v, aes(label = varname, 
    #                                     x = xvar, y = yvar, angle = angle, hjust = hjust), 
    #                    color = "darkred", size = varname.size)
    # g <- g + geom_text(data = df.v, aes(label = varname,
    #                                     x = xvar, y = yvar, angle = angle, hjust = hjust),
    #                    color = col_arrows, size = varname.size)
    g <- g + geom_text_repel(data = df.v, aes(label = varname,
                                              x = xvar, y = yvar, angle = angle),
                             color = col_arrows, size = varname.size)
    
  }
  return(g)
}

# the following functions are only necessary when interested in pca on reflectance
pca_explore_spec <- function(datalist) {
  dataset <- datalist$spectra
  label <- datalist$spec_meta$Species.code
  rownames(dataset) <- datalist$spec_meta$ID
  
  # library(stringr) #str_sub
  # label <- as.factor(str_sub(datalist$spec_meta$First_measurement,-3,-1))

  # Loadings of the bands in each PC (eigenvectors) 
  pca<-prcomp(dataset,scale=TRUE,center=TRUE)
  # summary(pca)
  x11(width=1500,height=720)
  par(mfrow=c(1,2))
  screeplot(pca) # idem as plot(pca) # It seems reasonable to withhold 3 PC's (also when looking at proportion of variance they explain)
  plot(pca,type='l')
  abline(h=1, col="blue")
  summary(pca)
  print(pca)
  
  
  # Indicate what should be visualised on the plot:
  # x11(width=1500,height=720)
  par(mfrow=c(1,2))
  g <- myggbiplot(pca, obs.scale = 1, var.scale = 1, #choices = 2:3,
                ellipse = FALSE, circle = FALSE, var.axes=T, labels=label) #,groups=label2) #grouping only works PFT_intra dataset!
  g <- g + scale_color_discrete(name = '') + theme_bw()
  print(g)
  
  palette_color <- rep(distinctColorPalette(nlevels(label)/4),4) # find a good color palette.
  palette_color <- rep(distinctColorPalette(10),4) # find a good color palette.
  palette_shape <- c(rep(15,10),rep(16,10),rep(17,10),rep(18,10))
  
  x11(width=1500,height=720)
  g <- myggbiplot(pca, obs.scale = 1, var.scale = 1, #choices = 2:3,
                  ellipse = F, circle = F, var.axes=F,alpha=0, # specify groups = vis and remove alpha = 0 when no shape differences
                  col_arrows = "dimgray",var.labels=F) + 
    geom_point(aes(colour=label,shape=label), size = 3) +
    scale_color_manual(name='Species',values=palette_color) +
    scale_shape_manual(name='Species',values=palette_shape)+
    theme(legend.direction = 'vertical', 
          legend.position = 'right') +
    theme(legend.position="none")+
    scale_x_continuous(breaks=seq(floor(min(pca$x[,1]))+floor(min(pca$x[,1]))%%2,ceiling(max(pca$x[,1]))+ceiling(max(pca$x[,1]))%%2,2))+
    scale_y_continuous(breaks=seq(floor(min(pca$x[,2]))+floor(min(pca$x[,2]))%%2,ceiling(max(pca$x[,2]))+ceiling(max(pca$x[,2]))%%2,2))+
    # scale_x_continuous(sec.axis = dup_axis()) + scale_y_continuous(sec.axis = dup_axis()) +
    theme_bw()+
    # ggtitle("Mean Optical traits of unmixed Table Green measurements")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  # g_reorder <- move_layers(g,idx=4L,position="bottom")
  # print(g_reorder)
  print(g)
  
  return(pca)
} 
library(Hmisc) # rcorr
library(tibble) # rownames_to_column
library(reshape) # melt
pca_explore_spec_SuppFig <- function(datalist, FunctionalTraits, name){
  dataset <- datalist$spectra
  label <- datalist$spec_meta$Species.code
  Spec_meta <- datalist$spec_meta 
  PFT_meta <- datalist$PFT_meta
  
  pca<-prcomp(dataset,scale=TRUE,center=TRUE)
  summary(pca)
  x11(width=1500,height=720)
  par(mfrow=c(1,2))
  screeplot(pca) # idem as plot(pca) # It seems reasonable to withhold 3 PC's (also when looking at proportion of variance they explain)
  plot(pca,type='l')
  abline(h=1, col="blue")
  
  ### Fig. 3 from Roth et al. 2016 ### 
  # Loadings of the bands in each PC (eigenvectors) 
  PC<-as.matrix(pca$rotation[,1])
  PC<-cbind(PC,as.matrix(pca$rotation[,2]),as.matrix(pca$rotation[,3]),as.matrix(pca$rotation[,4]))
  colnames(PC)<-c('PC1','PC2','PC3','PC4')
  
  # insert NA in the atmospheric windows: (based on hsdar mask)
  # mask(speclib_smooth)<-c(349,400,1340,1460,1780,1970,2400,2501) 
  #   - 350 tem 400 (voor atm window 1); 
  #   - 1340 tem 1460 (tussen PC 981 en 982);
  #   - 1780 tem 1970 (tussen PC 1342 en 1343); 
  #   - 2400 tem 2500 (na PC 1773)
  mask1<-matrix(rep(NA,((400-350))*4),ncol=4)
  mask2<-matrix(rep(NA,((1460-1340)-1)*4),ncol=4)
  mask3<-matrix(rep(NA,((1970-1780)-1)*4),ncol=4)
  mask4<-matrix(rep(NA,((2500-2400))*4),ncol=4)
  PC_Complete<-rbind(mask1,PC[1:941,],mask2,PC[942:1262,],mask3,PC[1263:1693,],mask4)
  PC_Complete<-cbind(c(350:2500),PC_Complete)
  colnames(PC_Complete)<-c('wl','PC1','PC2','PC3','PC4')
  PC_Complete<-as.data.frame(PC_Complete)
  
  # Labels:
  PoV <- (round(100*(pca$sdev^2/sum(pca$sdev^2)),1))[1:4] # Proportion of variance explained (also in summary) =  normalized standard deviations
  labels <- c(paste(paste("PC1 (", PoV[1], sep=""), "% explained var.)"), 
              paste(paste("PC2 (", PoV[2], sep=""), "% explained var.)"),
              paste(paste("PC3 (", PoV[3], sep=""), "% explained var.)"),
              paste(paste("PC4 (", PoV[4], sep=""), "% explained var.)"))
  
  x11(width=720,height=720)
  g <- ggplot(PC_Complete)+
    geom_line(aes(wl,PC4,color="PC4",linetype = 'PC4'),size=1)+
    geom_line(aes(wl,PC3,color="PC3",linetype = 'PC3'),size=1)+
    geom_line(aes(wl,PC2,color="PC2",linetype = 'PC2'),size=1)+
    geom_line(aes(wl,PC1,color="PC1",linetype = 'PC1'),size=1)+
    geom_hline(yintercept=0)+
    labs(x="wavelength (nm)", y="Loading")+
    scale_colour_manual(name = "Guide1", breaks = c("PC1","PC2","PC3","PC4"), values=c("blue", "green4", "red","orange"),
                        labels = labels)+
    scale_linetype_manual(name = "Guide1", breaks = c("PC1","PC2","PC3","PC4"),values=c("solid","dotted","dashed","dotdash"),
                          labels = labels)+
    theme(legend.justification = c(1,1), legend.position = c(0.99,0.99),
          legend.title = element_blank(),legend.key=element_blank(),legend.key.width = unit(3, "line"))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.5))
  print(g)
  ggsave(paste("PC_loadings_",".jpeg",sep=name),width = 15,height = 15,units=("cm"),dpi=600)
  # In this figure loadings are visualised rather than correlations. How can we calculate these? Same procedure as for figure 5.
  
  
  ### Table 4 from Roth et al. 2016 ### 
  # Coordinates of the individual observations on the PC's:
  PC_Coord<-as.matrix(pca$x[,1])
  PC_Coord<-as.data.frame(cbind(PC_Coord,as.matrix(pca$x[,2]),as.matrix(pca$x[,3]),as.matrix(pca$x[,4])))
  # PC_Coord<-cbind(PC_Coord,as.matrix(pca$x[,2]),as.matrix(pca$x[,3]))
  colnames(PC_Coord)<-c('PC1','PC2','PC3','PC4')
  
  # Correlation between each PC and categorical variables 
  # PC1_Meta<-cbind(PC_Coord[,c('PC1'),drop=F],meta_mean[,c('Species','DvpPhase','Patch','Species_Cond','Site','Date')])
  # PC2_Meta<-cbind(PC_Coord[,c('PC2'),drop=F],meta_mean[,c('Species','DvpPhase','Patch','Species_Cond','Site','Date')])
  # PC3_Meta<-cbind(PC_Coord[,c('PC3'),drop=F],meta_mean[,c('Species','DvpPhase','Patch','Species_Cond','Site','Date')])
  PC1_Meta<-cbind(PC_Coord[,c('PC1'),drop=F],Spec_meta[,c('Species.code','Site','Date')],
                  PFT_meta[,c('CSR_Hunt','Family','Order','Superorder','Group','Plant_growth_habit','Plant_lifespan_Pierce','Plant_lifespan_LEDA','Plant_growth_form_LEDA')])
  PC2_Meta<-cbind(PC_Coord[,c('PC2'),drop=F],Spec_meta[,c('Species.code','Site','Date'),drop=F],
                  PFT_meta[,c('CSR_Hunt','Family','Order','Superorder','Group','Plant_growth_habit','Plant_lifespan_Pierce','Plant_lifespan_LEDA','Plant_growth_form_LEDA')])
  PC3_Meta<-cbind(PC_Coord[,c('PC3'),drop=F],Spec_meta[,c('Species.code','Site','Date')],
                  PFT_meta[,c('CSR_Hunt','Family','Order','Superorder','Group','Plant_growth_habit','Plant_lifespan_Pierce','Plant_lifespan_LEDA','Plant_growth_form_LEDA')])
  PC4_Meta<-cbind(PC_Coord[,c('PC4'),drop=F],Spec_meta[,c('Species.code','Site','Date')],
                  PFT_meta[,c('CSR_Hunt','Family','Order','Superorder','Group','Plant_growth_habit','Plant_lifespan_Pierce','Plant_lifespan_LEDA','Plant_growth_form_LEDA')])
  
  PC1_Meta_lm<-sapply(PC1_Meta[,-1], function(x) summary(lm(PC1_Meta$PC1 ~ x)) )
  PC2_Meta_lm<-sapply(PC2_Meta[,-1], function(x) summary(lm(PC2_Meta$PC2 ~ x)) )
  PC3_Meta_lm<-sapply(PC3_Meta[,-1], function(x) summary(lm(PC3_Meta$PC3 ~ x)) )
  PC4_Meta_lm<-sapply(PC4_Meta[,-1], function(x) summary(lm(PC4_Meta$PC4 ~ x)) )
  PC_Meta_R2<-t(as.data.frame(PC1_Meta_lm['r.squared',]))
  PC_Meta_R2<-cbind(PC_Meta_R2,t(as.data.frame(PC2_Meta_lm['r.squared',])),t(as.data.frame(PC3_Meta_lm['r.squared',])),t(as.data.frame(PC4_Meta_lm['r.squared',])))
  PC_Meta_R2<-round(PC_Meta_R2,2)
  colnames(PC_Meta_R2)<-c('PC1','PC2','PC3','PC4')
  PC_Meta_R2adj<-t(as.data.frame(PC1_Meta_lm['adj.r.squared',]))
  PC_Meta_R2adj<-cbind(PC_Meta_R2adj,t(as.data.frame(PC2_Meta_lm['adj.r.squared',])),t(as.data.frame(PC3_Meta_lm['adj.r.squared',])),t(as.data.frame(PC4_Meta_lm['adj.r.squared',])))
  PC_Meta_R2adj<-round(PC_Meta_R2adj,2)
  colnames(PC_Meta_R2adj)<-c('PC1','PC2','PC3','PC4')
  write.table(PC_Meta_R2, file=paste("PC_categ_R2_",".csv",sep=name), sep =";", col.names=NA) 
  write.table(PC_Meta_R2adj, file=paste("PC_categ_R2adj_",".csv",sep=name), sep =";", col.names=NA) 
  
  
  ### Fig. 5 from Roth et al. 2016### 
  # Correlation between each PC and functional traits
  PC_trait<-cbind(PC_Coord,FunctionalTraits)
  corr_PCtrait<-rcorr(as.matrix(PC_trait),type="pearson")
  corr_PCtrait_r<-corr_PCtrait$r[5:nrow(corr_PCtrait$r),1:4]
  
  corr_PCtrait_r<-as.data.frame(t(corr_PCtrait_r))
  corr_PCtrait_r<-rownames_to_column(corr_PCtrait_r,"PC")
  
  df.long<-melt(corr_PCtrait_r,id.vars=c("PC"))
  x11(width=1200,height=720)
  p<- ggplot(df.long,aes(variable,value))+
    geom_col(aes(fill=variable))+
    scale_y_continuous(limits=c(-1,1))+
    scale_fill_brewer(name="Traits", palette="Spectral")+
    facet_wrap(~PC, ncol=4)+
    labs(y="Correlation")+
    theme_bw()+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  print(p)
  ggsave(paste("PC_trait_cor_",".jpeg",sep=name),width = 20,height = 20,units=("cm"),dpi=600)
  
  return(PC_Meta_R2)
}


#--------------------------------------------------------------------------------
# CORRELATION STRUCTURE
#--------------------
panel.cor.custom <- function(x, y, digits = 2, cex.cor,  ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  # txt <- paste("r= ", txt, sep = "")
  if(missing(cex.cor))  cex <- 0.5/strwidth(txt)
  # txt <- text(0.5, 0.5, txt)
  txt <- text(0.5, 0.5, txt,cex = cex*(abs(r)+0.4))
  # txt <- text(0.5, 0.5, txt,cex = cex)
  
  #p-value calculation --> kan je ook weglaten
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  # txt2 <- paste("p = ", txt2, sep = "")
  txt2 <- paste("(", ")", sep = txt2)
  if(p<0.01) txt2 <- paste("(", "0.01)", sep = "< ")
  # text(0.5, 0.2, txt2)
  txt <- text(0.5, 0.2, txt2,cex = cex*(abs(r)+0.2))
}
corr_traits <- function(data, transf, std, name,trait_names) {
  if (transf == ""){
    dataset <- data$datavalues
  } else if (transf == "_trans"){
    dataset <- data$datavalues_trans
  }else if (transf == "_backtrans"){
    dataset <- data$datavalues_backtrans
  }
  
  colnames(dataset) <- trait_names
  
  if (std == "_stand"){
    dataset<-data.stand(dataset,method='standardize',margin='column',plot=F) #test<-scale(dataset)
  }
  
  # Examine correlation among the functional traits: Scatter plot matrix
  dataset.r <- abs(cor(dataset)) # get correlations
  dataset.col <- dmat.color(dataset.r) # get colors
  # reorder variables so those with highest correlation are closest to the diagonal
  dataset.o <- order.single(dataset.r) 
  
  # x11(width=720,height=720)
  # cpairs(dataset, dataset.o, panel.colors=dataset.col, gap=.5,
  #        main="Variables Ordered and Colored by Correlation" , upper.panel=panel.cor)
  # png(paste("Correlation_traits",".png",sep=std), width = 8, height = 8, units = 'in', res = 600)
  # png(paste("Correlation_traits_",".png",sep=paste(transf,std,sep="")), width = 8, height = 8, units = 'in', res = 600)
  png(paste(paste("Correlation_traits_", transf ,sep=name),".png",sep=std), width = 8, height = 8, units = 'in', res = 600)
  
  # cpairs(dataset, dataset.o, panel.colors=dataset.col, gap=.5,
  #        main="Pearson correlation" , upper.panel=panel.cor.custom)
  cpairs(dataset, order=NULL, panel.colors=dataset.col, gap=.5,
         main=NULL , upper.panel=panel.cor.custom)
  dev.off()

  }  # non-transformed data (?)


#--------------------------------------------------------------------------------
# CLUSTERING
#--------------------
# ------- Cluster algorithm
B_spline <- function(x){
  library(fda)
  data <- x
  
  # Create spline for each continuous spectral interval?
  # mask(speclib_smooth)<-c(349,400,1340,1460,1780,1970,2400,2501) # for visualisation
  
  # http://www.psych.mcgill.ca/misc/fda/
  # Book
  # We will create cubic B-splines (order = 4; degree = 3)
  # ----------------- Spectral region 1 (400 - 1340 nm)
  # breaks_1 <-  seq(400,1340,20)
  # breaks_1 <-  c(400,430,460,seq(490,750,15),800)
  breaks_1 <-  c(seq(400,1340,30),1340)
  
  basisobj_1 <- create.bspline.basis(rangeval = c(400,1340), norder = 4, breaks = breaks_1)
  # basismatrix_1 <- eval.basis(breaks_1, basisobj_1)
  # basismatrix_1b <- predict(basisobj_1,breaks_1) # idem to basismatrix_1
  plot(basisobj_1, ylab = "B-spline basis functions B(nm)", xlab="nm", col="black")
  
  data_1 <- data[,as.character(breaks_1)]
  colnames(data_1) <- breaks_1
  fdo_1 <- Data2fd(argvals = breaks_1, y = t(data_1), basisobj = basisobj_1)
  coefs_1 <- fdo_1$coefs
  
  # Graphical evaluation of the fdo:
  # plotfit.fd(y = t(data_1), argvals = breaks_1, fdo_1, lty=1, lwd=2, residual=T)
  # plotfit.fd(y = t(data_1), argvals = breaks_1, fdo_1, lty=1, lwd=2, residual=F)
  
  # ----------------- Spectral region 2 (1460 - 1780 nm)
  breaks_2 <-  c(seq(1460,1780,30),1780)
  
  basisobj_2 <- create.bspline.basis(rangeval = c(1460,1780), norder = 4, breaks = breaks_2)
  
  plot(basisobj_2, ylab = "B-spline basis functions B(nm)", xlab="nm", col="black")
  
  data_2 <- data[,as.character(breaks_2)]
  colnames(data_2) <- breaks_2
  fdo_2 <- Data2fd(argvals = breaks_2, y = t(data_2), basisobj = basisobj_2)
  coefs_2 <- fdo_2$coefs
  
  # plotfit.fd(y = t(data_2), argvals = breaks_2, fdo_2, lty=1, lwd=2, residual=T)
  # plotfit.fd(y = t(data_2), argvals = breaks_2, fdo_2, lty=1, lwd=2, residual=F)
  
  # ----------------- Spectral region 3 (1970 - 2400 nm)
  breaks_3 <-  c(seq(1970,2400,30),2400)
  
  basisobj_3 <- create.bspline.basis(rangeval = c(1970,2400), norder = 4, breaks = breaks_3)
  
  plot(basisobj_3, ylab = "B-spline basis functions B(nm)", xlab="nm", col="black")
  
  data_3 <- data[,as.character(breaks_3)]
  colnames(data_3) <- breaks_3
  fdo_3 <- Data2fd(argvals = breaks_3, y = t(data_3), basisobj = basisobj_3)
  coefs_3 <- fdo_3$coefs
  
  # plotfit.fd(y = t(data_3), argvals = breaks_3, fdo_3, lty=1, lwd=2, residual=T)
  # plotfit.fd(y = t(data_3), argvals = breaks_3, fdo_3, lty=1, lwd=2, residual=F)
  
  # ----------------- Combine spline coefficients of the entire spectral range
  coefs_total <- rbind(coefs_1, coefs_2, coefs_3)
}

library(factoextra)
library(NbClust)
cluster_traits <- function(dataset, distance, method, nbclusters, stand) {
  if (stand == TRUE) {
    dataset.std<-data.stand(dataset,method='standardize',margin='column',plot=F) 
    dataset.range<-data.stand(dataset,method='range',margin='column',plot=F) 
    dataset.scale<-scale(dataset, center=T, scale=T) # idem as dataset.scale
    # check that we get mean of 0 and sd of 1
    colMeans(dataset.std)  # faster version of apply(scaled.dat, 2, mean)
    apply(dataset.std, 2, sd)
  } else if (stand == FALSE){
    dataset.std <- as.data.frame(dataset)
  }
  
  ### a) Agglomerative hierarchical clustering
  # Ward Hierarchical Clustering
  dist_SAM <- function (x){
    y <- fDiss(x, method="cosine",scaled=F,center=F) # use the cosine/SAM function from fDISS package --> different results!!! seem more correct (also used for black table analyses)
    colnames(y) <- rownames(y) <- rownames(x) #this is the only line
    z <- as.dist(y)
  }

  if (distance == "euc") {
    d <- dist(dataset.std, method = "euclidean") # distance matrix; methods: euclidean
  } else if (distance == "Gower") {
    library(proxy)
    # summary(pr_DB) # a list of all metrics available
    # pr_DB$get_entry("Gower") # get more information about a particular one
    d<-dist(dataset.std,method="Gower")#idem to:
    dist_Gower <- function (x){
      y <- dist(x, method="Gower")
    }
    d2<-dist_Gower(dataset.std)
  } else if (distance == "SAM") {
    # pr_DB$get_entry("cosine")
    d2<-dist(as.matrix(dataset.std, method="cosine")) #  does not seem correct -- should we still take the cos of the results (cfr. formula)?
    d<-dist_SAM(as.matrix(dataset.std))
    d3 <- fDiss(as.matrix(dataset.std), method = "cosine",center=F,scaled=F) # idem to d
  } else if (distance == "Manhattan"){
     d<-dist(dataset.std,method="manhattan",diag=T,upper=T)
  } else if (distance == "SID"){
    d<-sid(em,mode="density",center=F,scaled=F)$sid
  } else if (distance == "SAM+MAN"){
    # https://stackoverflow.com/questions/21332959/summing-2-distance-matrices-for-getting-a-third-overall-distance-matrix-ecolo
    # https://datascience.stackexchange.com/questions/27407/clustering-with-multiple-distance-measures
    d_sam<-dist_SAM(as.matrix(dataset.std))
    d_man<-dist(dataset.std,method="manhattan",diag=T,upper=T)
    
    # Extract relevant code form Fuse code (github)
    weights=c(4/6, 2/6) # frequency of traits estimated from brightness-normalized vs. conventional plsr
    dots <- list(d_sam, d_man)
    dots <- lapply(dots, as.dist, diag = FALSE, upper = TRUE)
    D <- do.call("cbind", dots)
    range01 <- function(x){(x-min(x))/(max(x)-min(x))} # different from fuse code: we scale between 0 and 1
    D01 <- apply(D, 2, FUN = range01)
    retval <- rowSums(D01 * weights)
    class(retval) <- "dist"
    attr(retval, "Labels") <- attr(dots[[1]], "Labels")
    attr(retval, "Size") <- attr(dots[[1]], "Size")
    attr(retval, "Diag") <- FALSE
    attr(retval, "Upper") <- FALSE
    attr(retval, "method") <- "fuse"
    attr(retval, "weights") <- weights
    attr(retval, "call") <- match.call()
    d <- retval
    # library(analogue)
    # d <- analogue::fuse(d_sam, d_man, weights=c(5/7, 2/7))
    # d5050 <- analogue::fuse(d_sam, d_man) # 50%-50% SAM and MAN
  } else if (distance == "Bspline"){
    BSpline_coef <- B_spline(dataset)
    d <- dist(t(BSpline_coef), method = "euclidean") # distance matrix; methods: euclidean
  }
  
  if (method == "Ward") {
    fit <- hclust(d, method="ward.D2") # Apparantly Ward.D does not correctly implement the Ward algorithm
    x11(width=720,height=720)
    # hclus.table(fit) #summary
    hclus.scree(fit) # this plot can help to decide upon number of clusters
    nbclusters<-nbclusters # specify the number of desired clusters (k)
    groups <- cutree(fit, k=nbclusters)
    
    # Statistically define or defend no. of clusters
    # http://stat.sys.i.kyoto-u.ac.jp/prog/pvclust/
    ####### Dynamic Hybrid Clustering
    DHC <- cutreeHybrid(dendro = fit, distM = as.matrix(d), minClusterSize = 2)
    # DHC$values returns the same output as:
    # DC <- cutreeDynamic(dendro = fit, distM = as.matrix(d), minClusterSize = 2)
    nbclusters_DHC <- length(unique(DHC$labels))
    
    x11(width=720,height=720)
    plot(fit)
    rect.hclust(fit, k=nbclusters_DHC, border="red") # draw dendogram with red borders around the clusters
    # rect.hclust(fit, h=heights) # draw dendogram with red borders around the clusters 

    # Determine no. of clusters based upon multiple indices
    x11(width=720,height=720)
    multimethod <- NbClust(data = scale(dataset), diss = NULL, distance = "euclidean",
                    min.nc = 2, max.nc = 15, method = "ward.D2")$Best.partition
    x11(width=720,height=720)
    plot(fit)
    rect.hclust(fit, k=max(multimethod), border="blue") # draw dendogram with blue borders around the clusters
    
    # Kmeans clustering instead of hierarchical clustering
    x11(width=720,height=720)
    Kmeans <- NbClust(data = scale(dataset), diss = NULL, distance = "euclidean",
                           min.nc = 2, max.nc = 15, method = "kmeans")$Best.partition
    # Check if the kmeans algorihtm allocates observations to the same group as NbClust
    Kmeans_alg <- kmeans(scale(dataset), max(Kmeans))$cluster
  }
  
  # # Idem to
  # library(dendextend)
  # groups <- dataset.std %>%  # Retrieve groupings
  #   # scale %>% 
  #   dist(method = "euclidean") %>% 
  #   hclust(method = "ward.D2") %>% # Hierarchical clustering 
  #   cutree(k=nbclusters) 
  # fit <- dataset.std %>%   # For more information on clustering performance we need the variable 'fit'
  #   # scale %>% 
  #   dist(method = "euclidean") %>% 
  #   hclust(method = "ward.D2") 

    output <- list ("fit" = fit,
                  "groups" = groups,
                  "groups_DHC" = DHC$labels,
                  "groups_NbClust" = multimethod,
                  "groups_Kmeans" = Kmeans,
                  "groups_Kmeans_alg" = Kmeans_alg,
                  "distancematrix" = d)
    return(output)
}

# ------- Visualise clusters on pca
library(RColorBrewer)
create_color_list <- function(scheme){
  if (scheme == "Dark2"){
    getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
    colors=getPalette(8)
    # Extent the colour list with three more colours:
    colors <- c(colors, "#63B8FF", "#7A378B", "#000080")
  } else if (scheme == "Set1") {
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    colors=getPalette(9) 
    # getPalette2 = colorRampPalette(brewer.pal(8,"Accent"))
    # colors2 = getPalette2(8)
    # the yellow color is too bright so we will change it + add extra colours
    # colors <- c(colors[1:5],colors[7:9], "#F0027F",colors2) #, "#ffc414"
    colors <- c(colors[1:5],colors[7:9],"#FFC125", "#00008B", "#20B2AA","#E9967A", "#228B22", "#8B7500")
  }
  
  return(colors)
}


library(pca3d)  # pca3d
# require(installr)
# library(magick)
library(rgl) # legend3d

clusters_on_pca <- function(pca, datalist, cluster_output, clusterdecision, form, name,colorscheme) {
  
  if (form == "orig"){
    dataset <- datalist$datavalues
    label <- datalist$label
  } else if (form == "trans"){
    dataset <- datalist$datavalues_trans
    label <- datalist$label
  } else if (form == "backtrans"){
    dataset <- datalist$datavalues_backtrans
    label <- datalist$label
  } else if (form == "spectra"){
    dataset <- datalist$spectra
    label <- datalist$spec_meta$Species.code
    rownames(dataset) <- datalist$spec_meta$ID
  } else if (form == "Bspline"){
    dataset <- t(B_spline(datalist$spectra))
    label <- datalist$spec_meta$Species.code
  }
  
  if (clusterdecision == "DHC"){
    groups<-cluster_output$groups_DHC
  } else if (clusterdecision == "nbclusters"){
    groups<-cluster_output$groups
  } else if (clusterdecision == "NbClust"){
    groups<-cluster_output$groups_NbClust
  } else if (clusterdecision == "Kmeans"){
    groups<-cluster_output$groups_Kmeans
  } else if (clusterdecision == "Kmeans_alg"){
    groups<-cluster_output$groups_Kmeans_alg
  }
  
  if(form != "spectra"){
    write.table(as.data.frame(round(pca$loadings,2)), file=paste("PCA_loadings_",".csv",sep=paste(name,form,sep="")), sep =";", col.names=NA)
  }
  
  colors <- create_color_list(colorscheme)
  
  # PC 1 and 2 with and without variable names 
  x11(width=1500,height=720)
  par(mfrow=c(1,2))
  g <- myggbiplot(pca, obs.scale = 1, var.scale = 1, #choices = 2:3,
                  ellipse = FALSE, circle = FALSE, var.axes=T, labels=label,groups=as.factor(groups), #as.factor(membership$groups)
                  col_arrows = "gray18",var.labels=T) #grouping only works PFT_intra dataset!
  print(g)  
  g <- g + scale_color_manual(name='Species',values=colors) +
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  g <- g + theme(legend.position = "none")
  g_reorder <- move_layers(g,idx=2L,position="bottom")
  print(g_reorder)
  ggsave(paste(paste("ClusterOnPCA",name,sep="_"),"_PC12.jpeg",sep=clusterdecision),width = 20,height = 15,units=("cm"),dpi=600)
  
  x11(width=1500,height=720)
  par(mfrow=c(1,2))
  g <- myggbiplot(pca, obs.scale = 1, var.scale = 1, #choices = 2:3,
                  ellipse = FALSE, circle = FALSE, var.axes=T, labels=label,groups=as.factor(groups), #as.factor(membership$groups)
                  col_arrows = "gray18",var.labels=F) #grouping only works PFT_intra dataset!
  g <- g + scale_color_manual(name='Species',values=colors) +
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  g <- g + theme(legend.position = "none")
  g_reorder <- move_layers(g,idx=2L,position="bottom")
  print(g_reorder)
  ggsave(paste(paste("ClusterOnPCA",name,sep="_"),"_PC12_withoutVarnames.jpeg",sep=clusterdecision),width = 20,height = 15,units=("cm"),dpi=600)
} 

clusters_on_pca_ellipses <- function(pca, CI,group_order_nb, text, datalist, cluster_output, clusterdecision, form, name,colorscheme) {
  
  if (form == "orig"){
    dataset <- datalist$datavalues
    label <- datalist$label
  } else if (form == "trans"){
    dataset <- datalist$datavalues_trans
    label <- datalist$label
  } else if (form == "backtrans"){
    dataset <- datalist$datavalues_backtrans
    label <- datalist$label
  } else if (form == "spectra"){
    dataset <- datalist$spectra
    label <- datalist$spec_meta$Species.code
    rownames(dataset) <- datalist$spec_meta$ID
  } else if (form == "Bspline"){
    dataset <- t(B_spline(datalist$spectra))
    label <- datalist$spec_meta$Species.code
  }
  
  if (clusterdecision == "DHC"){
    groups<-cluster_output$groups_DHC
  } else if (clusterdecision == "nbclusters"){
    groups<-cluster_output$groups
  } else if (clusterdecision == "NbClust"){
    groups<-cluster_output$groups_NbClust
  } else if (clusterdecision == "Kmeans"){
    groups<-cluster_output$groups_Kmeans
  } else if (clusterdecision == "Kmeans_alg"){
    groups<-cluster_output$groups_Kmeans_alg
  }
  
  if(form != "spectra"){
    write.table(as.data.frame(round(pca$loadings,2)), file=paste("PCA_loadings_",".csv",sep=paste(name,form,sep="")), sep =";", col.names=NA)
  }
  
  colors <- create_color_list(colorscheme)
  
  # ---------------
  # 3d pca movie
  # Ellipses cannot be drawn in 3D when there are only 2 observations in a group
  # https://stackoverflow.com/questions/53489831/3d-pca-plotting-error-the-leading-minor-of-order-3-is-not-positive-definite
  wd = getwd()
  count <- table(groups)
  ellips_labels <- paste (rep(text,dim(count)), as.factor(group_order_nb), sep="")
  info <- as.data.frame(cbind(group_order_nb,ellips_labels,color=colors[1:length(group_order_nb)]))
  info$group_order_nb <- as.numeric(as.character(info$group_order_nb))
  info_ordered <- arrange(info, group_order_nb)
  
  if (is.element(2, count) != TRUE){
    pca3d(pca,group=as.factor(groups), biplot=T,biplot.vars = 6,
          shape = 1,palette=colors,
          show.ellipses = T, ellipse.ci=CI, show.plane=F)
    legend3d("topleft",legend = info_ordered$ellips_labels, pch=16,col=as.character(info_ordered$color),box.lty=0)
    # snapshotPCA3d(file="3dPCA_ePFT_04.png")
    
  } else {  
  obs2=c()
  rmcolor=c()
  for (i in 1:max(groups)){
    if (count[i] == 2 ) {
      obs <- which(groups == names(count)[i])
      obs2 <- c(obs2, obs)
      rmcolor <- c(rmcolor,as.integer(names(count)[i]))
    }
  }
  groups2 <- groups[-obs2]
  pca2 <- pca
  pca2$x <- pca2$x[-obs2,]
  colors2 <- colors[-rmcolor]
  info2 <- info[-unique(groups[obs2]),]
  info_ordered2 <- arrange(info2, group_order_nb)
  
  pca3d(pca2,group=as.factor(groups2), biplot=T,biplot.vars = 6, 
        shape = 1,palette=colors2,
        show.ellipses = T, ellipse.ci=CI, show.plane=F)
  # plot3d(pca2$x[obs2,'PC1'],pca2$x[obs2,'PC2'],pca2$x[obs2,'PC3'], type = "p", col=colors[rmcolor],add=T)
  legend3d("topleft",legend = info_ordered2$ellips_labels, pch=16,col=as.character(info_ordered2$color),box.lty=0)
  # snapshotPCA3d(file="3dPCA_ePOT_04.png")
  
   }
  makeMoviePCA(dir=wd, clean=T, type=".gif",
               movie=paste(paste("Clusters_on_PCA_",clusterdecision,sep=name),"_movie",sep=toString(CI)))
  
  # ----------
  # 2d pca
  # PC 1 and 2 with and without variable names 
  x11(width=1500,height=720)
  par(mfrow=c(1,2))
  g <- myggbiplot(pca, obs.scale = 1, var.scale = 1,ellipse.prob=CI, #choices = 2:3,
                  ellipse = T, circle = FALSE, var.axes=T, groups=as.factor(groups), #labels=label,as.factor(membership$groups)
                  col_arrows = "gray18",var.labels=T) 
  print(g)  
  g <- g + scale_color_manual(name='Species',values=colors) +
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  g <- g + theme(legend.position = "none")
  g_reorder <- move_layers(g,idx=2L,position="bottom")
  print(g_reorder)
  ggsave(paste(paste("ClusterOnPCA_ellipses_",name,sep=toString(CI)),"_PC12.jpeg",sep=clusterdecision),width = 20,height = 15,units=("cm"),dpi=600)
  
  
  
  x11(width=1500,height=720)
  par(mfrow=c(1,2))
  g <- myggbiplot(pca, obs.scale = 1, var.scale = 1,ellipse.prob = CI, #choices = 2:3,
                  ellipse = T, circle = FALSE, var.axes=T,groups=as.factor(groups), #as.factor(membership$groups)
                  col_arrows = "gray18",var.labels=F) #grouping only works PFT_intra dataset!
  g <- g + scale_color_manual(name='Species',values=colors) +
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  g <- g + theme(legend.position = "none")
  g_reorder <- move_layers(g,idx=2L,position="bottom")
  print(g_reorder)
  ggsave(paste(paste("ClusterOnPCA_ellipses",name,sep=toString(CI)),"_PC12_withoutVarnames.jpeg",sep=clusterdecision),width = 20,height = 15,units=("cm"),dpi=600)
} 
# set clusterdecision to "nbclusters" if you want to fix the number of clusters of cluster_output to the number identified in functional clustering


# ------- Visualise clusters as dendrogram
# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
# library(dplR)
library(dplyr)
library(dendextend)
plot_dend <- function(method, cluster_output, clusterdecision, name, colorscheme) {
  if (clusterdecision == "DHC"){
    groups<-cluster_output$groups_DHC
  }  else if (clusterdecision == "NbClust"){
      groups<-cluster_output$groups_NbClust
  } else if (clusterdecision == "nbclusters"){
    groups<-cluster_output$groups
  }
  
  colors <- create_color_list(colorscheme)
  # Palette=c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", 
  #              "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
  #              "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
  #              "#8A7C64", "#599861") # https://stackoverflow.com/questions/21352683/randomising-qualitative-colours-for-large-sets-in-ggplot
  # colors=Palette[1:nbclusters]
  colors2=rep(c("white"),10)
  
  # dataset.std<-data.stand(dataset,method='standardize',margin='column',plot=F) #test<-scale(dataset)
  # dend <- dataset.std %>% 
  #   # scale %>% 
  #   dist(method = distance) %>% # euclidean or Gower --> possible to add different one? SAM
  #   hclust(method = method) %>% # Hierarchical clustering
  #   as.dendrogram
  # # dend<-as.dendrogram(fit)
  
  # Also possible to direcly use the distance matrix: more convenient because
  #   - we might be using self-defined distance metrices (SAM)
  #   - we do not always standardize data (e.g. SAM on reflectance) 
  dend_orig <- cluster_output$distancematrix %>% 
    hclust(method = method) %>% # Hierarchical clustering
    as.dendrogram
  # dend<-as.dendrogram(fit)
  
  col_GROUP <- colors
  col_GROUP <- colors[groups]
  col_GROUP <- col_GROUP[order.dendrogram(dend_orig)] 
  
  # dend <- dend %>% 
  #   set("branches_k_color", value=col_GROUP) %>%
  #   set("labels_colors", value = col_GROUP) %>%# no labels visible = make them white --> colors2 :-)
  #   set("branches_lwd", 1) %>%
  #   set("labels_cex", c(0.5))
  # 
  dend <- dend_orig %>%
    # color_branches(value=colors, k=nbclusters) %>%
    # set("branches_k_color", value=colors, k=nbclusters) %>%
    set("branches_k_color", value=col_GROUP) %>%
    set("branches_lwd", 1) %>%
    # set("labels_colors", value = colors,k=nbclusters) %>% # no labels visible = make them white --> colors2 :-)
    set("labels_colors", value = col_GROUP) %>% # no labels visible = make them white --> colors2 :-)
    set("labels_cex", c(0.5))
  # rotate(groups)
  
  dend_plot <- dend
  labels(dend_plot)<-paste(" ",substr(labels(dend_plot),1,6))
  
  x11(width=500,height=720)
  ggd1 <- as.ggdend(dend_plot)
  g<-ggplot(ggd1, horiz=T)
  # g <- g + theme(legend.position="none",
  #                                    axis.text.x = element_text(),  # show x-axis labels
  #                                    axis.ticks.x = element_line(), # show x-axis tick marks
  #                                    axis.line.x = element_line())  # show x-axis lines
  
  print(g)
  ggsave(paste(name,".jpeg",sep=clusterdecision),width = 12,height = 20,units=("cm"),dpi=600)
  
  output <- list (dend_orig = dend_orig,
                  dend_col = dend,
                  col = col_GROUP, 
                  groups = groups)
  return(output)
} # dendextend and ggplot

# ------- Visualise 2 dendrograms as tanglegram
dend_compare_tanglegram <- function(cluster1, dend1_list, dend2_list,clusterdecision, name, colorsheme1){
  # Entanglement is a measure between 1 (full entanglement) and 0 (no entanglement). A lower entanglement coefficient corresponds to a good alignment
  col_lines <- dend1_list$col
  dend1 <- dend1_list$dend_col 
  dend2 <- dend2_list$dend_col
  
  dend1_ordered <- dend1
  order.dendrogram(dend1_ordered) <- 1:length(dend1_list$groups)
  dend2_ordered <- dend2
  dend2_ordered <- match_order_by_labels(dend2_ordered, dend1_ordered) 
  
  untangle_methods <- c("labels", "ladderize","random", "step1side", "step2side", "DendSer")
  untanglement_matrix <- cbind(untangle_methods,rep(NA,6),rep(NA,6))
  colnames(untanglement_matrix) <- c("method", "untangle", "untangle_ordered")
  for (i in 1:length(untangle_methods)){
    untanglement_matrix[i,2] <- untangle(dend1, dend2, method = untanglement_matrix[i,1]) %>% entanglement
    untanglement_matrix[i,3] <- untangle(dend1_ordered, dend2_ordered, method = untanglement_matrix[i,1]) %>% entanglement # should be idem to the above
  }
  best_method <- untanglement_matrix[which(untanglement_matrix[,2] == min(untanglement_matrix[,2])),1]
  best_untangle <- untangle(dend1, dend2, method=best_method)
  best_entanglement <- as.numeric(untanglement_matrix[which(untanglement_matrix[,"method"] == best_method),2])
  
  # https://bioinformatics.stackexchange.com/questions/920/how-do-i-generate-a-color-coded-tanglegram/932#932           
  # https://www.researchgate.net/publication/299486353_Environmental_signals_in_radial_growth_stable_isotope_variations_and_nutrient_concentration_of_trees_from_different_forest_ecosystems_in_southern_Ecuador/figures?lo=1
  # labels_dend1 <- as.data.frame(labels(dend1))
  labels_dend1_untangled <- as.data.frame(labels(best_untangle[[1]]))
  Observations <- labels(cluster1$distancematrix)
  
  if (clusterdecision == "DHC"){
    Groups1<-cluster1$groups_DHC
  }  else if (clusterdecision == "NbClust"){
    Groups1<-cluster1$groups_NbClust
  } else if (clusterdecision == "nbclusters"){
    Groups1<-cluster1$groups
  }
  # Groups1 <- cluster1$groups_DHC
  colors <- create_color_list(colorsheme1)
  # col_Group1 <- colors[cluster1$groups_DHC]
  col_Group1 <- colors[Groups1]
  ObsGroups <- as.data.frame(cbind(Observations, Groups1, col_Group1))
  ObsGroups$Observations <- as.character(ObsGroups$Observations)
  ObsGroups$Groups1 <- as.integer(ObsGroups$Groups1)
  ObsGroups$col_Group1 <- as.character(ObsGroups$col_Group1)
  colors_ordered <- merge(labels_dend1_untangled,ObsGroups,by.x="labels(best_untangle[[1]])", by.y="Observations",sort=F)
  # merge function performs the same as using "order.dendrogram":
  # colors_ordered2 <- col_Group1[order.dendrogram(best_untangle[[1]])] 
  
  # labels(best_untangle[[1]]) <- substr(labels(best_untangle[[1]]),1,6)
  # labels(best_untangle[[2]]) <- substr(labels(best_untangle[[2]]),1,6)
  
  
  tiff(paste("Tanglegram_",".tiff",sep=name),res=300,width=8,height=8,units='in')
  # x11(width=720,height=720)
  best_untangle %>% plot(highlight_distinct_edges = FALSE, # Turn-off dashed lines
                         common_subtrees_color_lines = F, # Turn-off line colors
                         common_subtrees_color_branches = F, # Color common branches
                         color_lines=colors_ordered$col_Group1,
                         match_order_by_labels = T,
                         lwd=2,
                         edge.lwd = 1,
                         intersecting = F,
                         lab.cex = 0.7,
                         # k_labels=7,
                         # k_branches = 7,
                         main = paste(paste("entanglement (", ") =", sep=best_method), round(best_entanglement, 2))
  ) 
  dev.off()
  
  
  ### ----- Correlation matrix (exploration for the function "dend_compare_metrics")
  # --> Does ordering and/or untangling affect correlation meausures?
  
  # dend_list <- dendlist(dend1, dend2)
  # dend_list_ordered <- dendlist(dend1_ordered, dend2_ordered)
  # dend1_untangled <- best_untangle[[1]]
  # order.dendrogram(dend1_untangled) <- 1:length(dend1_list$groups)
  # dend2_untangled <- best_untangle[[2]]
  # dend2_untangled <- match_order_by_labels(dend2_untangled, dend1_untangled)
  # dend_list_ordered_untangled <- dendlist(dend1_untangled, dend2_untangled)
  
  # From tests it is clear that is not necessary to match_order_by_labels and/or to untangle prior to metric calculations
  # cors_coph<-cor.dendlist(dend_list, method = "cophenetic")
  # cors_baker<-cor.dendlist(dend_list, method = "baker")
  # library(corrplot)
  # x11(width=720,height=720)
  # corrplot(cors_coph, "pie", "lower",title="Cophenetic correlation",addCoef.col = "black")
  # x11(width=720,height=720)
  # corrplot(cors_baker, "pie", "lower",title="Baker correlation",addCoef.col = "black")
  
  ### Unique edges per tree (in red)
  # x11(width=720,height=720)
  # dend_diff(dend1, dend2)
}

#--------------------------------------------------------------------------------
# CLUSTERING EVALUATION

##### 1. Cluster validity and stability:
#   - agglomerative coefficient
#   - cophenetic correlation with original distance matrix
ClusterQuality <- function(cluster_output, clusterdecision, name) {
  d <- cluster_output$distancematrix
  fit <- cluster_output$fit
  
  if (clusterdecision == "DHC"){
    groups<-cluster_output$groups_DHC
  } else if (clusterdecision == "NbClust"){
    groups<-cluster_output$groups_NbClust
  } else if (clusterdecision == "nbclusters"){
    groups<-cluster_output$groups
  }
  
  ##### Amount of clustering structure
  ## Agglomerativce coefficient:
  # high value means that observations quickly agglomerate into distinct clusters that later agglomerate into a sinlge cluster at much greater dissimilarities
  # For each observation i, denote by m(i) its dissimilarity to the first cluster it is merged with, divided by the dissimilarity of the merger in the final step of the algorithm. 
  # The agglomerative coefficient is the average of all 1 - m(i). It can also be seen as the average width (or the percentage filled) of the banner plot.
  
  # http://strata.uga.edu/8370/lecturenotes/clusterAnalysis.html:
  # The agglomerative coefficient measures the dissimilarity of an object to the first cluster it joins, divided by the dissimilarity of the final merger in the cluster analysis, averaged across all samples. 
  # Low values reflect tight clustering of objects, larger values indicate less well-formed clusters. The agglomerative coefficient increases with sample sizes, making comparisons among data sets difficult.
  # --> this definition seems to be not entirely the same as the one in the help function, where the final result represent average(1-m(i)) instead of average(m(i))
  agglcoef<-coef.hclust(fit)
  
  ## Cophenetic correlation coefficient to evaluate how well the clustering represents the original dissimilarity matrix
  # cophonetic distance = intergroup dissimilarity at which two observations are first combined in a single cluster
  # Cophenetic correlation = correlation between dissimilarities in original dissimilarity matrix and cophonetic distances
  # Cophenetic correlation > .75 are considered good.
  cophcorr<-cor(d,cophenetic(fit))  
  x11(width=720,height=720)
  plot(d,cophenetic(fit))
  hclus.cophenetic(d,fit)
  
  Table <- matrix(c(agglcoef,cophcorr))
  rownames(Table)<-c("AGG","COP")
  colnames(Table)<-name
  
  output <- list ("cophcorr" = cophcorr,
                  "agglcoef" = agglcoef)
  return(output)
  
} 

#   - Jaccard similarity based on bootstrap resampling scheme
library(fpc) # clusterboot
ClusterStability <- function(cluster_output, clusterdecision) {
  
  ##### Cluster stability
  d <- cluster_output$distancematrix
  fit <- cluster_output$fit
  
  if (clusterdecision == "DHC"){
    groups<-cluster_output$groups_DHC
  } else if (clusterdecision == "NbClust"){
    groups<-cluster_output$groups_NbClust
  } else if (clusterdecision == "nbclusters"){
    groups<-cluster_output$groups
  }
  
  d.boot<-clusterboot(d, B=100, bootmethod=c('boot','subset'), bscompare=FALSE,
                      multipleboot=TRUE, clustermethod=disthclustCBI, method='ward.D2',
                      k=length(unique(groups)), cut='number', count=FALSE)
  print(d.boot)
  
  # This line will cut the dendrogram according to no. of clusters, not to what was defined by hybrid approach!
  # So ... clusterboot does not work with dynamic hybrid clustering...?
  # A quick check suggested that the dynamic cut just applies a basic cutting based on no. of clusters ... not very clear ...
  
  # d.boot<-clusterboot(d, B=100, bootmethod=c('boot','subset'), bscompare=FALSE,
  #                     multipleboot=TRUE, clustermethod=disthclustCBI, method='ward.D2',
  #                     k=length(unique(groups)), cut=cutreeDynamic(distM = as.matrix(d), minClusterSize = 2),
  #                     count=FALSE)
  # # DHC <- cutreeHybrid(dendro = fit, distM = as.matrix(d), minClusterSize = 2)
  # DHC$values 
  # print(d.boot)
  
  # The data is resampled using several schemes (bootstrap, subsetting, jittering, replacement of points by noise) and the Jaccard similarities of the
  # original clusters to the most similar clusters in the resampled data are computed. The mean over these similarities is used as an index of the stability of a cluster.
  # Guidelines for interpretation of the bootstrap(!) resampling scheme:
  #   Jaccard similarity value <= 0.5: indication of a "dissolved cluster".
  #   Jaccard similarity value <0.6: clusters should not be trusted
  #   Jaccard similarity value 0.6-0.75: indicating patterns in the data, but which points exactly should belong to these clusters is highly doubtful
  #   Jaccard similarity value => 0.75: valid, stable cluster
  #   Jaccard similarity value => 0.85: "Highly stable" clusters 
  
  ##### Silhouette
  # A value close to 1 indicates that the data point is assigned to a very
  # appropriate cluster. If Sil is close to zero, it means that that data
  # point could be assigned to another closest cluster as well because it
  # is equidistant from both the clusters. If Si is close to -1, it
  # means that data is misclassified and lies somewhere in between the clusters.
  win.graph()
  sil1 <- silhouette(groups,d)
  plot(sil1, nmax = 80, cex.names = 0.5)
  
  win.graph()
  sil2 <- silhouette(cutree(fit, k=max(groups)),d)
  plot(sil2, nmax = 80, cex.names = 0.5)
}


##### 2. Describe clusters: species membership
library(plyr) # count
membership <- function(dataset, cluster_output, clusterdecision,name){
  if (clusterdecision == "DHC"){
    groups<-cluster_output$groups_DHC
  } else if (clusterdecision == "NbClust"){
    groups<-cluster_output$groups_NbClust
  } else if (clusterdecision == "nbclusters"){
    groups<-cluster_output$groups
  }
  
  Sp_code <- as.character(dataset$meta[,c('Species.code')])
  ProfileCat<-cbind(as.data.frame(groups),Sp_code)
  ProfileCat$groups<-as.factor(ProfileCat$groups)
  
  categ_atrr <- c("Species.code","Site","Date","CSR_Hunt","Family","Order","Superorder","Group",
                   "Plant_growth_habit","Plant_lifespan_Pierce","Plant_lifespan_LEDA","Plant_growth_form_LEDA")
  ProfileCat <- cbind(ProfileCat, dataset$meta[,categ_atrr])

  # Choose which categorical attribute you would like to visualise <- ProfileCat
    attr<-ProfileCat$Sp_code
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    colors <- getPalette(length(unique(attr)))  
    
    c<-ggplot(ProfileCat,aes(groups, ..count..))+
      geom_bar(aes(fill=attr))+
      scale_fill_manual(name="Species",values = colors)+
      theme_bw()
      
    x11(width=720,height=720)
    print(c)
    ggsave(paste(name,".jpeg",sep=clusterdecision),width = 15,height = 15,units=("cm"),dpi=600)
    
    count_per_group <- plyr::count(ProfileCat,"groups")
    
    # paste CSR columns etc. for better overview
    attr_per_group <- ProfileCat[order(ProfileCat$groups),]
    
    output <- list ("ProfileCat" = ProfileCat,
                    "count_per_group" = count_per_group,
                    "attr_per_group" = attr_per_group)
    return(output)
  }


##### 3. Describe clusters: Standardized trait profile of each group
#        = Fig. 10: from Roth et al. (2016)
Fig_ProfileGroups_ordered <- function(cluster_output,clusterdecision, dataset, name, colorscheme, group_order_nb, text, trait_labels) {
  if (clusterdecision == "DHC"){
    groups<-cluster_output$groups_DHC
  } else if (clusterdecision == "nbclusters"){
    groups<-cluster_output$groups
  }
  nbclusters <- length(unique(groups))
  
  colnames(dataset) <- trait_labels
  
  
  Traits_scaled<-apply(dataset, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
  Profile_scaled<-as.data.frame(cbind(groups,Traits_scaled)) #idem to: Profile_scaled<-as.data.frame(cbind(groups,dataset.range))
  Profile_scaled.long<-melt(Profile_scaled,id.vars=c("groups"))
  Profile_scaled.long$variable = with(Profile_scaled.long, factor(variable, levels = rev(levels(variable))))
  Profile_scaled.long$groups <- as.factor(Profile_scaled.long$groups)
  
  group_order_original <- seq(1:length(group_order_nb))
  Profile_scaled$groups_ordered <- Profile_scaled$groups
  Profile_scaled.long$groups_ordered <- Profile_scaled.long$groups
  for (i in 1:length(group_order_original)){
    Profile_scaled.long[,"groups_ordered"] <- replace( Profile_scaled.long[,"groups_ordered"], Profile_scaled.long$groups==group_order_original[i],group_order_nb[i])
    Profile_scaled[,"groups_ordered"] <- replace(Profile_scaled[,"groups_ordered"], Profile_scaled$groups==group_order_original[i],group_order_nb[i])
  }
  
  colors <- create_color_list(colorscheme)
  colors<- as.data.frame(colors)
  colors$order_original <- seq(1:dim(colors)[1])
  colors <- colors[1:length(group_order_nb),]
  colors$colors_ordered <- colors$colors
  for (i in 1:length(group_order_original)){
    colors[group_order_nb[i],"colors_ordered"] <- colors[group_order_original[i],"colors"] 
  }
  colors_ordered <- as.character(colors$colors_ordered)
  
  library(plyr)
  Groups_nb_obs <- ddply(.data=Profile_scaled, 
                         .(groups_ordered), 
                         summarise, 
                         n=paste("n =", length(groups)))
  facet_labels <- paste(paste(rep(text,nrow(Groups_nb_obs)), as.factor(Groups_nb_obs$groups), sep=""),
                        paste(rep(" (",nrow(Groups_nb_obs)),Groups_nb_obs$n, sep=""),
                        paste(rep(")",nrow(Groups_nb_obs))),sep="")
  levels(Profile_scaled.long$groups_ordered) <- facet_labels
  
  s<-ggplot(Profile_scaled.long, aes(variable,value,fill=as.factor(groups_ordered)))+
    geom_boxplot(data=Profile_scaled.long[,2:3],aes(variable),fill="white",color="darkgrey",outlier.shape=NA)+
    geom_boxplot()+
    scale_fill_manual(values=colors_ordered)+
    # scale_fill_brewer(palette="Spectral")+
    coord_flip()+
    facet_wrap(~groups_ordered, ncol=5)+
    scale_y_continuous(limits=c(0,1),breaks=c(0, 0.25, 0.5, 0.75, 1),labels=c("0","0.25","0.5","0.75","1"))+
    theme_bw()+
    # geom_text(data=Groups_nb_obs, aes(x=1.8, y=5, label=n), 
    #           colour="black", inherit.aes=FALSE, parse=FALSE)+
    theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())
  # x11(width=1000,height=720)
  print(s)
  ggsave(paste(name,".jpeg",sep=clusterdecision),width = 20,height = 14,units=("cm"),dpi=600)
  
}

##### 4. Describe clusters: Spectral profile of each group
library(hsdar)
library(dplyr) # filter
Fig_ProfileGroupsRefl_ordered <- function(cluster_output,clusterdecision, datalist, name, colorscheme,group_order_nb) {
  if (clusterdecision == "DHC"){
    groups<-cluster_output$groups_DHC
  } else if (clusterdecision == "nbclusters"){
    groups<-cluster_output$groups
  }
  
  Profile <- as.data.frame(cbind(groups,labels(cluster_output$distancematrix)))
  group_order_original <- seq(1:length(group_order_nb))
  Profile$groups_ordered <- Profile$groups
  for (i in 1:length(group_order_original)){
    Profile[,"groups_ordered"] <- replace(Profile[,"groups_ordered"], Profile$groups==group_order_original[i],group_order_nb[i])
  }
  groups_ordered <- as.integer(Profile[,"groups_ordered"])
    
  dataset <- datalist$spectra
  speclib <- speclib(dataset,as.numeric(colnames(dataset)))
  SI(speclib) <- datalist$spec_meta
  idSpeclib(speclib) <- as.character(datalist$spec_meta$ID)
  
  colors <- create_color_list(colorscheme)
  colors<- as.data.frame(colors)
  colors$order_original <- seq(1:dim(colors)[1])
  colors <- colors[1:length(group_order_nb),]
  colors$colors_ordered <- colors$colors
  for (i in 1:length(group_order_original)){
    colors[group_order_nb[i],"colors_ordered"] <- colors[group_order_original[i],"colors"]
  }
  colors_ordered <- as.character(colors$colors_ordered)
  transp <- 0.3
  
  #### ------------ Plot mean spectra and sd of the different groups on separate graphs
  
  # tiff(paste(name,"_ordered.tiff",sep=clusterdecision),width = 50,height = 15,units=("cm"),res=600)
  tiff(paste(name,"_ordered.tiff",sep=clusterdecision),width = 30,height = 20,units=("cm"),res=600)
  
  # x11(width=2500,height=1000)
  # fig <- par(mfrow=c(2,5),mai = c(0.3, 0.4, 0.3, 0.3))
  fig <- par(mfrow=c(4,3),mai = c(0.3, 0.4, 0.3, 0.3))
  for (i in 1:length(unique(groups_ordered))){
    List_IDs <- dplyr::filter(as.data.frame(Profile), groups_ordered == as.character(i))$V2
    idx_int <- SI(speclib)$ID %in% List_IDs
    speclib_plot<-speclib[idx_int,]
    
    # plot(c(),c(),axes=FALSE,ylim=c(0,1.0),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
    plot(c(),c(),axes=FALSE,ylim=c(0,0.6),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
    # 1) Colour the area between the sd-borders
    wl_plot <- c(seq(400, 1340, 1), seq(1460,1780,1), seq (1970, 2400, 1))
    mean_speclib_plot <- apply(speclib_plot, FUN = mean, na.rm = TRUE)
    sd_speclib_plot <- apply(speclib_plot, FUN = sd, na.rm = TRUE) 
    spectra2plot_speclib <- rbind(spectra(mean_speclib_plot) + spectra(sd_speclib_plot),
                                  spectra(mean_speclib_plot),
                                  spectra(mean_speclib_plot) - spectra(sd_speclib_plot))
    
    xx<-c(wl_plot,(rev(wl_plot)))
    yy_speclib<-c(spectra2plot_speclib[1,],rev(spectra2plot_speclib[3,]))
    polygon(xx,yy_speclib,col=alpha(c(colors_ordered[i]),transp),border=NA,new=F)
    
    # 2) Plot mean spectrum and sd on top of the areas
    plot(speclib_plot, FUN=mean, lwd=1, ylim=c(0,1), col=c(colors_ordered[i]), new=FALSE)
    
    #3) Mask the atmospheric absorption windows
    xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
    xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
    yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
    yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
    polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
    polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
    
    # 3) Plot graphical details: axes, etc. if wanted
    abline(h=0, lwd=2)
    # axis(1, labels=F, lwd=1, lwd.ticks=0, pos=0)
    # axis(2, labels=F, lwd=1, lwd.ticks=0, pos=400)
    # axis(1,seq(350,2500,150),seq(350,2500,150),font=2)
    # axis(2,seq(0.0,1.05,0.1),seq(0.0,1.05,0.1),font=2,las=2)
    # mtext('Reflectance',side=2,font=2,line=3)
    # mtext('Wavelength (nm)',side=1,font=2,line=3)
    
    axis(1,seq(500,2500,500),seq(500,2500,500),font=2,cex.axis=1.5)
    axis(2,seq(0.0,1.05,0.2),seq(0.0,1.05,0.2),font=2,las=2,cex.axis=1.5)
    # mtext('Reflectance',side=2,font=2,line=3,cex=1.5)
    # mtext('Wavelength (nm)',side=1,font=2,line=3,cex=1.5)
  }
  dev.off()
  
  #### ------------ Plot mean spectra of the different groups on 1 graph
  
  # tiff(paste(name,"_ordered_together.tiff",sep=clusterdecision),width =25,height = 15,units=("cm"),res=600)
  tiff(paste(name,"_ordered_together.tiff",sep=clusterdecision),width =18,height = 13,units=("cm"),res=600)
  # x11(width=2500,height=1000)
  par(mar=c(5.1,5.1,2.1,2.1))
  plot(c(),c(),axes=FALSE,ylim=c(0,0.6),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
  for (i in 1:length(unique(groups_ordered))){
    List_IDs <- dplyr::filter(as.data.frame(Profile), groups_ordered == as.character(i))$V2
    idx_int <- SI(speclib)$ID %in% List_IDs
    speclib_plot<-speclib[idx_int,]
    plot(speclib_plot, FUN=mean, lwd=3, ylim=c(0,1), col=c(colors_ordered[i]), new=FALSE)
  }
  xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
  xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
  yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
  yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
  polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
  polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
  legend(x=2100, y=0.6, legend = group_order_original, col = colors_ordered, lwd=3, bty="n", cex=1.5)
  # abline(h=0, lwd=2)
  axis(1,seq(350,2500,350),seq(350,2500,350),font=2,cex.axis=1.5)
  axis(2,seq(0.0,0.6,0.2),seq(0.0,0.6,0.2),font=2,las=2,cex.axis=1.5)
  mtext('Reflectance',side=2,font=2,line=3,cex=1.5)
  mtext('Wavelength (nm)',side=1,font=2,line=3,cex=1.5)
  dev.off()
  
}


##### 5. Describe clusters: Test for significant differences between groups
# Here, you have to choose which type of test and visualisation you want to perform!
# Standardize traits prior to analysis

# Kruskal Wallis test: does not assume normality, nor homoscedasticity
KruskalWallis <- function(dataset, groups, wd1, wd2, name) {
  dataset.std<-data.stand(dataset,method='standardize',margin='column',plot=F) #test<-scale(dataset)
  groups <- groups
  
  # http://userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html
  
  kruskal.wallis.alpha=0.05
  
  kruskal.wallis.table <- data.frame()
  for (i in 1:dim(dataset.std)[2]) {
    ks.test <- kruskal.test(dataset.std[,i], g=groups)
    # Store the result in the data frame
    kruskal.wallis.table <- rbind(kruskal.wallis.table,
                                  data.frame(id=names(dataset.std)[i],
                                             p.value=ks.test$p.value
                                  ))
    # Report number of traits tested
    cat(paste("Kruskal-Wallis test for ",names(dataset.std)[i]," ", i, "/", 
              dim(dataset.std)[2], "; p-value=", ks.test$p.value,"/n", sep=""))
  }
  kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
  kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, 
                                      size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
  # kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value,
  #                                                    decreasing=FALSE), ]
  kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
  kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
  
  
  x11(width=720,height=720)
  plot(kruskal.wallis.table$p.value,
       kruskal.wallis.table$E.value,
       main='Multitesting corrections',
       xlab='Nominal p-value',
       ylab='Multitesting-corrected statistics',
       log='xy',
       col='blue',
       panel.first=grid(col='#BBBBBB',lty='solid'))
  lines(kruskal.wallis.table$p.value,
        kruskal.wallis.table$FWER,
        pch=20,col='darkgreen', type='p'
  )
  lines(kruskal.wallis.table$p.value,
        kruskal.wallis.table$q.value,
        pch='+',col='darkred', type='p'
  )
  abline(h=kruskal.wallis.alpha, col='red', lwd=2)
  legend('topleft', legend=c('E-value', 'p-value', 'q-value'), col=c('blue', 'darkgreen','darkred'), lwd=2,bg='white',bty='o')
  
  
  ## We will plot all results, also non-significant results
  # last.significant.element <- max(which(kruskal.wallis.table$q.value <= kruskal.wallis.alpha))
  # selected <- 1:last.significant.element
  # diff.cat.factor <- kruskal.wallis.table$id[selected]
  # diff.cat <- as.vector(diff.cat.factor)
  # print(kruskal.wallis.table[selected,])
  diff.cat <- colnames(dataset.std)

  #Now we plot traits significantly different between the categories
  kruskal.wallis.table$q.value_print <- paste(paste("(q =",round(kruskal.wallis.table$q.value,2)),")",sep="")
  kruskal.wallis.table[which(kruskal.wallis.table$q.value < 0.01),"q.value_print"] <- "(q < 0.01)"
  
  setwd(wd1)
  write.table(kruskal.wallis.table, file=paste("KW_table_",".csv",sep=name), sep =";", col.names=T, row.names = F) 
  # ggsave(paste(name,".jpeg",sep=clusterdecision),width = 20,height = 20,units=("cm"),dpi=600)
  
  df<-NULL
  for(i in diff.cat){
    # tmp<-data.frame(dataset.std[,i],groups,rep(paste(i," q = ",round(kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value"],5),sep=""),dim(dataset.std)[1]))
    tmp<-data.frame(dataset.std[,i],groups,rep(paste(i,kruskal.wallis.table[kruskal.wallis.table$id==i,"q.value_print"],sep=" "),dim(dataset.std)[1]))
    if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)}
  }
  colnames(df)<-c("Value","Group","Trait")
  df$Group <- as.factor(df$Group)

  # colors <- create_color_list("Dark2")
  # x11(width=1200,height=720)
  # p<-ggplot(df,aes(Type,Value,colour=Type))+ylab("Normalised")
  # p<-p+geom_boxplot()+geom_jitter()+theme_bw()+
  #   scale_color_manual(values=colors)+
  #   guides(fill = guide_legend(reverse = TRUE))+
  #   facet_wrap( ~ Trait , scales="free", ncol=3)
  # p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  # p<-p+coord_flip()
  # print(p)
  
  Profile.long <- df
  Profile.long$Group<-as.factor(Profile.long$Group)
  Profile.long$Group = with(Profile.long, factor(Group, levels = rev(levels(Group))))
  
  colors <- create_color_list("Dark2")
  colors <- colors[1:length(unique(groups))]
  
  p<-ggplot(Profile.long, aes(Group,Value, colour=Group))+ #,fill=groups
    # ylab("Normalised")+
    geom_boxplot()+
    geom_jitter()+
    # scale_fill_manual(values=rev(colors))+
    scale_colour_manual(values=rev(colors))+
    coord_flip()+
    facet_wrap(~Trait,scales="free")+
    # guides(fill = guide_legend(reverse = T))+
    theme_bw()+
    theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank()) 
  x11(width=720,height=720)
  print(p)
  setwd(wd2)
  ggsave(paste(name,"jpeg",sep="."),width = 20,height = 20,units=("cm"),dpi=600)
  
  p<-ggplot(Profile.long, aes(Group,Value, colour=Group))+ #,fill=groups
    # ylab("Normalised")+
    geom_violin()+
    geom_jitter()+
    # scale_fill_manual(values=rev(colors))+
    scale_colour_manual(values=rev(colors))+
    coord_flip()+
    facet_wrap(~Trait,scales="free")+
    # guides(fill = guide_legend(reverse = T))+
    theme_bw()+
    theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank()) 
  print(p)
  setwd(wd2)
  ggsave(paste(name,"violin.jpeg",sep="_"),width = 20,height = 20,units=("cm"),dpi=600)
  
  # x11(width=720,height=720)
  # p<-ggplot(Profile.long, aes(Value, colour=Group,fill=Group))+ #,fill=groups
  #   geom_density(adjust=2)+
  #   scale_colour_manual(values=rev(colors))+
  #   scale_fill_manual(values=alpha(rev(colors),0.3))+
  #   theme_bw()+
  #   theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank()) 
  # print(p)  
  # ggsave(paste(name,"_density.jpeg",sep=clusterdecision),width = 20,height = 20,units=("cm"),dpi=600)
  
  
  return(kruskal.wallis.table)
} 
# https://rcompanion.org/rcompanion/d_06.html
# If the Kruskal-Wallis test is significant, a post-hoc analysis can be performed 
# to determine which levels of the independent variable differ from each other level.
# options:
#   - Dunn test for multiple comparisons
#   - Nemenyi test for multiple comparisons
#   - Pairwise Mann-Whitney U-tests = Wilcoxon rank test
#  Zar (2010) states that the Dunn test is appropriate for groups with unequal numbers of observations.
library(FSA) # dunnTest
Dunn_paired_ordered <- function(dataset, groups, partitioning,group_order_nb){
  # dataset.std<-data.stand(dataset,method='standardize',margin='column',plot=F) #test<-scale(dataset)
  dataset.std <- dataset
  
  groups <- as.factor(groups)
  
  
  Profile<-as.data.frame(cbind(dataset,groups)) 
  # Profile.long<-melt(Profile,id.vars=c("groups"))
  # Profile.long$variable = with(Profile.long, factor(variable, levels = rev(levels(variable))))
  # Profile.long$groups <- as.factor(Profile.long$groups)
  
  group_order_original <- seq(1:length(group_order_nb))
  Profile$groups_ordered <- Profile$groups
  # Profile.long$groups_ordered <- Profile.long$groups
  for (i in 1:length(group_order_original)){
    # Profile.long[,"groups_ordered"] <- replace( Profile.long[,"groups_ordered"], Profile.long$groups==group_order_original[i],group_order_nb[i])
    Profile[,"groups_ordered"] <- replace(Profile[,"groups_ordered"], Profile$groups==group_order_original[i],group_order_nb[i])
  }
  groups_ordered <- Profile[,"groups_ordered"]
  
  
  dunn.alpha=0.06
  dunn.table.bonfer <- data.frame()
  dunn.table <- data.frame() 
  
  for (i in 1:dim(dataset.std)[2]) {
    d.test <- dunnTest(Profile[,i], g=groups_ordered,method="bh")
    
    sign.pairs.bonfer <- d.test$res[which(d.test$res$P.adj < dunn.alpha),]
    nb.sign.pairs.bonfer <- length(which(d.test$res$P.adj < dunn.alpha))
    if (nb.sign.pairs.bonfer != 0){
      sign.pairs.bonfer <- cbind(sign.pairs.bonfer, trait = names(dataset.std)[i])
      dunn.table.bonfer <- rbind(dunn.table.bonfer, sign.pairs.bonfer)
    } else {
      sign.pairs.bonfer <- 0
      dunn.table.bonfer <- dunn.table.bonfer
    }
    
    sign.pairs <- d.test$res[which(d.test$res$P.unadj < dunn.alpha),]
    nb.sign.pairs <- length(which(d.test$res$P.unadj < dunn.alpha))
    if (nb.sign.pairs != 0){
      sign.pairs <- cbind(sign.pairs, trait = names(dataset.std)[i])
      dunn.table <- rbind(dunn.table, sign.pairs)
    } else {
      sign.pairs <- 0
      dunn.table <- dunn.table
    }
  }
  
  dunn.table.ordered <- dunn.table[order(dunn.table[,'trait'],dunn.table[,"Comparison"]),]
  dunn.table.ordered$Comparison <- paste(rep("groups",dim(dunn.table.ordered)[1]),dunn.table.ordered$Comparison)
  dunn.table.ordered$Z <- round(dunn.table.ordered$Z, 2)
  # dunn.table.ordered$P.unadj <- round(dunn.table.ordered$P.unadj,3)
  for (i in 1:nrow(dunn.table.ordered)) {
    if (dunn.table.ordered[i,"P.unadj"] < 0.001) {
      dunn.table.ordered[i,"P.unadj"] <- c("< 0.001")
    } else {
      dunn.table.ordered[i,"P.unadj"] <- round(as.numeric(dunn.table.ordered[i,"P.unadj"]),3)
    }}
  
  dunn.table.bonfer.ordered <- dunn.table.bonfer[order(dunn.table.bonfer[,'trait'],dunn.table.bonfer[,"Comparison"]),]
  dunn.table.bonfer.ordered$Comparison <- paste(rep("groups",dim(dunn.table.bonfer.ordered)[1]),dunn.table.bonfer.ordered$Comparison)
  dunn.table.bonfer.ordered$Z <- round(dunn.table.bonfer.ordered$Z, 2)
  for (i in 1:nrow(dunn.table.bonfer.ordered)) {
    if (dunn.table.bonfer.ordered[i,"P.adj"] < 0.001) {
      dunn.table.bonfer.ordered[i,"P.adj"] <- c("< 0.001")
    } else {
      dunn.table.bonfer.ordered[i,"P.adj"] <- round(as.numeric(dunn.table.bonfer.ordered[i,"P.adj"]),3)
    }}
  for (i in 1:nrow(dunn.table.bonfer.ordered)) {
    if (dunn.table.bonfer.ordered[i,"P.unadj"] < 0.001) {
      dunn.table.bonfer.ordered[i,"P.unadj"] <- c("< 0.001")
    } else {
      dunn.table.bonfer.ordered[i,"P.unadj"] <- round(as.numeric(dunn.table.bonfer.ordered[i,"P.unadj"]),3)
    }}
  write.table(dunn.table.ordered, file=paste("Dunn_table_",".csv",sep=partitioning), sep =";", col.names=T, row.names = F) 
  write.table(dunn.table.bonfer.ordered, file=paste("Dunn_table_Bonfer_",".csv",sep=partitioning), sep =";", col.names=T, row.names = F) 
  
  f <- ggplot(dunn.table, aes(x=trait))+
    geom_bar(col="gray")+
    scale_y_continuous(name="number of significant pairs",breaks = seq(0,20,2),
                       expand = c(0,0))+
    theme_classic()+
    # theme(axis.text=element_text(size=11), plot.margin = unit(c(1,1,1,1), units = , "cm"))
    # theme_bw()+
    theme(panel.border = element_blank(), panel.background = element_rect(fill="white"),
          axis.ticks.x=element_blank(),
          axis.text=element_text(size=11), plot.margin = unit(c(1,1,1,1), units = , "cm"))
  print(f)
  ggsave(paste("Dunn_traits_","_ordered.jpeg",sep=partitioning),width = 20,height = 13,units=("cm"),dpi=600)
  
  f <- ggplot(dunn.table.bonfer, aes(x=trait))+
    geom_bar(col="gray")+
    scale_y_continuous(name="number of significant pairs",breaks = seq(0,20,2),
                       expand = c(0,0))+
    theme_classic()+
    # theme(axis.text=element_text(size=11), plot.margin = unit(c(1,1,1,1), units = , "cm"))
    # theme_bw()+
    theme(panel.border = element_blank(), panel.background = element_rect(fill="white"),
          axis.ticks.x=element_blank(),
          axis.text=element_text(size=11), plot.margin = unit(c(1,1,1,1), units = , "cm"))
  print(f)
  ggsave(paste("Dunn_traits_","_ordered_bonfer.jpeg",sep=partitioning),width = 20,height = 13,units=("cm"),dpi=600)
  
  # dunnTest.bonf.LNC <- dunnTest(dataset.std[,traitlist[1]], g=groups,method="bonferroni")
  # dunnTest.bonf.LPC <- dunnTest(dataset.std[,traitlist[2]], g=groups,method="bonferroni")
  # dunnTest.bonf.SLA <- dunnTest(dataset.std[,traitlist[3]], g=groups,method="bonferroni")
  # dunnTest.bonf.LDMC <- dunnTest(dataset.std[,traitlist[4]], g=groups,method="bonferroni")
  # dunnTest.bonf.Chltotal_area <- dunnTest(dataset.std[,traitlist[5]], g=groups,method="bonferroni")
  # dunnTest.bonf.Height <- dunnTest(dataset.std[,traitlist[6]], g=groups,method="bonferroni")
  # dunnTest.bonf.LA <- dunnTest(dataset.std[,traitlist[7]], g=groups,method="bonferroni")
  # 
  # output <- list (dunn.table = dunn.table.ordered,
  #                 dunn.table.bonfer = dunn.table.bonfer.ordered,
  #                 dunnTest.bonf.LNC_area = dunnTest.bonf.LNC_area,
  #                 dunnTest.bonf.LPC_area = dunnTest.bonf.LPC_area,
  #                 dunnTest.bonf.SLA = dunnTest.bonf.SLA,
  #                 dunnTest.bonf.LDMC = dunnTest.bonf.LDMC,
  #                 dunnTest.bonf.Chltotal_area = dunnTest.bonf.Chltotal_area,
  #                 dunnTest.bonf.Height = dunnTest.bonf.Height,
  #                 # dunnTest.bonf.Seed_weight = dunnTest.bonf.Seed_weight,
  #                 dunnTest.bonf.LA = dunnTest.bonf.LA)
  # 
  # return(output)
}
# In the above: methods was changed to BH (false discovery rate) instead of Bonferroni
# In Thomas et al. (2018) Bonferroni-adjustment of values is not explicitly coded
# In Thomas et al. (2018) traits are log transformed ("Stdvalue") and normalised (StdValueScaled), also before clustering
library(forcats) # fct_rev
Wilcoxon_paired <- function(dataset, groups){ 
  dataset.std<-data.stand(dataset,method='standardize',margin='column',plot=F) #test<-scale(dataset)
  groups <-as.factor(groups)
  
  # Wilcox.test.bonf.LNC_area <- pairwise.wilcox.test(dataset.std[,"LNC_mass"], g=groups,p.adj="bonf")
  # Wilcox.test.bonf.LPC_area <- pairwise.wilcox.test(dataset.std[,"LPC_mass"], g=groups,p.adj="bonf")
  # Wilcox.test.bonf.SLA <- pairwise.wilcox.test(dataset.std[,"SLA"], g=groups,p.adj="bonf")
  # Wilcox.test.bonf.LDMC <- pairwise.wilcox.test(dataset.std[,"LDMC"], g=groups,p.adj="bonf")
  # Wilcox.test.bonf.Chltotal_area <- pairwise.wilcox.test(dataset.std[,"Chltotal_mass"], g=groups,p.adj="bonf")
  # Wilcox.test.bonf.Height <- pairwise.wilcox.test(dataset.std[,"Height"], g=groups,p.adj="bonf")
  # # Wilcox.test.bonf.Seed_weight <- pairwise.wilcox.test(dataset.std[,"Seed_weight"], g=groups,p.adj="bonf")
  # Wilcox.test.bonf.LA <- pairwise.wilcox.test(dataset.std[,"LA"], g=groups,p.adj="bonf")
  
  
  Wilcox.test.bonf.LNC <- pairwise.wilcox.test(dataset.std[,"LNC_mass"], g=groups)
  Wilcox.test.bonf.LPC <- pairwise.wilcox.test(dataset.std[,"LPC_mass"], g=groups)
  Wilcox.test.bonf.SLA <- pairwise.wilcox.test(dataset.std[,"SLA"], g=groups)
  Wilcox.test.bonf.LDMC <- pairwise.wilcox.test(dataset.std[,"LDMC"], g=groups,p.adj="bonf")
  Wilcox.test.bonf.Chltotal <- pairwise.wilcox.test(dataset.std[,"Chltotal_mass"], g=groups)
  Wilcox.test.bonf.Height <- pairwise.wilcox.test(dataset.std[,"Height"], g=groups)
  Wilcox.test.bonf.LA <- pairwise.wilcox.test(dataset.std[,"LA"], g=groups)
  
  # for (i in 1:dim(dataset)[2]){
  #   Wilcox.test <- pairwise.wilcox.test(dataset[,i], g=groups)
  #   Wilcox.test.bonf <- pairwise.wilcox.test(dataset[,i], g=groups,p.adj="bonf")
  # }
  
  # library(rcompanion)
  # PT1 = fullPTable(Wilcox.test.bonf$p.value)
  # library(multcompView)
  # multcompLetters(PT1, 
  #                 compare="<", 
  #                 threshold=0.05,
  #                 Letters=letters,
  #                 reversed = FALSE)
  
  output <- list (Wilcox.test.bonf.LNC_mass = Wilcox.test.bonf.LNC_mass,
                  Wilcox.test.bonf.LPC_mass = Wilcox.test.bonf.LPC_mass,
                  Wilcox.test.bonf.SLA = Wilcox.test.bonf.SLA,
                  Wilcox.test.bonf.LDMC = Wilcox.test.bonf.LDMC,
                  Wilcox.test.bonf.Chltotal_mass = Wilcox.test.bonf.Chltotal_mass,
                  Wilcox.test.bonf.Height = Wilcox.test.bonf.Height,
                  # Wilcox.test.bonf.Seed_weight = Wilcox.test.bonf.Seed_weight,
                  Wilcox.test.bonf.LA = Wilcox.test.bonf.LA)
  
  return(output)
}

# Other option: assume normality but not homoscedasticity
### Welch's anova for unequal variances, but assumes normal distribution within each group (<-> KW)
# https://rcompanion.org/rcompanion/d_05.html
# http://web.pdx.edu/~newsomj/uvclass/ho_posthoc.pdf
# https://www.researchgate.net/post/How_do_I_choose_a_post_hoc_test_when_equal_variances_are_not_assumed_in_SPSS
Welch <- function(dataset, groups, trans){
  # dataset.std<-data.stand(dataset,method='standardize',margin='column',plot=F) #test<-scale(dataset)
  # dataset.std <- dataset
  if (trans == "LA"){
    dataset$LA <- log10(dataset$LA)
  } else if (trans == "LA_log"){
    dataset$LA_log <- log10(dataset$LA_log)
  }
  
  
  Welch.table <- data.frame()
  for (i in 1:dim(dataset)[2]) {
    Welch.test <- oneway.test(dataset[,i]~groups, var.equal=FALSE)
    # Store the result in the data frame
    Welch.table <- rbind(Welch.table,
                         data.frame(id=names(dataset)[i],
                                    p.value=Welch.test$p.value
                         ))
    # Report number of traits tested
    cat(paste("Welch test for ",names(dataset)[i]," ", i, "/", 
              dim(dataset)[2], "; p-value=", Welch.test$p.value,"/n", sep=""))
  }
  return(Welch.table)
} # stats package
# post-hoc test: Games-Howell test
# https://rpubs.com/aaronsc32/games-howell-test
# library(userfriendlyscience) # posthocTGH
# install.packages("remotes")
# remotes::install_github("GegznaV/BioStat")
library(biostat) #  posthoc_anova
# Biostats also returns letters:
# Differences between means that share a letter are not statistically significant (0.95 confidence level).
# https://rdrr.io/github/GegznaV/BioStat/man/make_cld.html
# https://rdrr.io/github/GegznaV/BioStat/src/R/make_cld.R
GH_paired_ordered <- function(dataset, groups, partitioning,group_order_nb, trans){
  # dataset.std<-data.stand(dataset,method='standardize',margin='column',plot=F) #test<-scale(dataset)
  # dataset.std <- dataset
  if (trans == "LA"){
    dataset$LA <- log10(dataset$LA)
  } else if (trans == "LA_log"){
    dataset$LA_log <- log10(dataset$LA_log)
  }
  
  groups <- as.factor(groups)
  
  Profile<-as.data.frame(cbind(dataset,groups)) 
  
  group_order_original <- seq(1:length(group_order_nb))
  Profile$groups_ordered <- Profile$groups
  for (i in 1:length(group_order_original)){
    Profile[,"groups_ordered"] <- replace(Profile[,"groups_ordered"], Profile$groups==group_order_original[i],group_order_nb[i])
  }
  groups_ordered <- Profile[,"groups_ordered"]
  
  GH.alpha=0.06
  GH.table <- data.frame()
  GH.table_letters <- data.frame()
  
  for (i in 1:dim(dataset)[2]) {
    # Option 1: userfriendlyscience package
    # GH.test <- posthocTGH(y=dataset[,i], x=groups_ordered, method=c("games-howell"),
    #                       conf.level = 0.95, digits=2, p.adjust="none",
    #                       formatPvalue = TRUE)$output$games.howell
    # sign.pairs.GH <- GH.test[which(GH.test$p < GH.alpha),]
    # nb.sign.pairs.GH <- length(which(GH.test$p < GH.alpha))
    # if (nb.sign.pairs.GH != 0){
    #   sign.pairs.GH <- cbind("Comparison"=rownames(sign.pairs.GH), sign.pairs.GH, trait = names(dataset)[i])
    #   GH.table <- rbind(GH.table, sign.pairs.GH)
    # } else {
    #   sign.pairs.GH <- 0
    #   GH.table <- GH.table
    # }
    
    # ----
    # Option 2: biostat package: also returns letters of significance :-)
    GH.test <- posthoc_anova(dataset[,i] ~ groups_ordered,
                             method = "Games-Howell")
    GH_matrix <- GH.test$output$result
    sign.pairs.GH <- GH_matrix[which(GH_matrix$p < GH.alpha),]
    nb.sign.pairs.GH <- length(which(GH_matrix$p < GH.alpha))
    if (nb.sign.pairs.GH != 0){
      sign.pairs.GH <- cbind(sign.pairs.GH, trait = names(dataset)[i])
      GH.table <- rbind(GH.table, sign.pairs.GH)
      
      sign.pairs.GH_letters <- cbind(make_cld(GH.test),trait = names(dataset)[i])
      GH.table_letters <- rbind(GH.table_letters, sign.pairs.GH_letters)
    } else {
      sign.pairs.GH <- 0
      GH.table <- GH.table
    }
  }
  
  for (i in 1:nrow(GH.table)) {
    if (GH.table[i,"p"] < 0.001) {
      GH.table[i,"p"] <- c("< 0.001")
    } else {
      GH.table[i,"p"] <- round(as.numeric(GH.table[i,"p"]),3)
    }}
  
  GH.table$groups <- paste(rep("groups",dim(GH.table)[1]),GH.table$groups)
  write.table(GH.table, file=paste("GH_table_",".csv",sep=partitioning), sep =";", col.names=T, row.names = F) 
  write.table(GH.table_letters, file=paste("GH_table_letters_",".csv",sep=partitioning), sep =";", col.names=T, row.names = F) 
  
  return(GH.table)
}

# Visualisation
Fig_ProfileTraits_Dunn <- function(cluster_output,clusterdecision, dataset, name, colorscheme, group_order_nb, text, facet_names, trans) {
  if (clusterdecision == "DHC"){
    groups<-cluster_output$groups_DHC
  } else if (clusterdecision == "nbclusters"){
    groups<-cluster_output$groups
  }
  groups <- as.factor(groups)
  
  
  Profile<-cbind(dataset,groups)  # dataset with real units
  # Profile<-as.data.frame(cbind(dataset,groups)) 
  
  Profile.long<-melt(Profile,id.vars=c("groups"))
  Profile.long$groups<-as.factor(Profile.long$groups)
  Profile.long$groups = with(Profile.long, factor(groups, levels = rev(levels(groups))))
  
  group_order_original <- seq(1:length(group_order_nb))
  Profile$groups_ordered <- Profile$groups
  Profile.long$groups_ordered <- Profile.long$groups
  for (i in 1:length(group_order_original)){
    Profile.long[,"groups_ordered"] <- replace( Profile.long[,"groups_ordered"], Profile.long$groups==group_order_original[i],group_order_nb[i])
    Profile[,"groups_ordered"] <- replace(Profile[,"groups_ordered"], Profile$groups==group_order_original[i],group_order_nb[i])
  }
  groups_ordered <- Profile[,"groups_ordered"]
  
  # add "ePFT" or "ePOT" to x-label
  Profile.long$groups_ordered2 = Profile.long$groups_ordered
  Profile.long$groups_ordered2 <- paste(rep(text,nrow(Profile.long)), as.factor(Profile.long$groups_ordered2), sep=" ")
  
  colors <- create_color_list(colorscheme)
  colors<- as.data.frame(colors)
  colors$order_original <- seq(1:dim(colors)[1])
  colors <- colors[1:length(group_order_nb),]
  colors$colors_ordered <- colors$colors
  for (i in 1:length(group_order_original)){
    colors[group_order_nb[i],"colors_ordered"] <- colors[group_order_original[i],"colors"] 
  }
  colors_ordered <- as.character(colors$colors_ordered)
  
  
  # ------------------
  # 1-way ANOVA with Tuckey test --> differences can be added with letters
  # https://www.researchgate.net/post/R_How_to_add_labels_for_significant_differences_on_boxplot_ggplot2
  # library(agricolae) # HSD.test
  # library(tidyverse) # as_tibble
  # library(magrittr) # %>%
  # i=1
  # Tukey_test <- aov(dataset[,i]~groups_ordered) %>%
  #   HSD.test("groups_ordered", group=T)
  #   .$groups %>%
  #   as_tibble(rownames="groups_ordered") %>%
  #   rename("Letters_Tukey"="groups_ordered") 
  #   select(-hwy) %>%
  
  # --------------------
  # Games-Howel (GH) and Dunn posthoc test, and add letters?
  # GH.alpha=0.51
  # GH.table <- data.frame()
  # GH.table_letters <- data.frame()
  
  dunn.alpha=0.051
  dunn.table <- data.frame()
  dunn.table_letters <- data.frame()
  dunn.table.bh <- data.frame()
  dunn.table.bh_letters <- data.frame()
  
  for (i in 1:ncol(dataset)) {
    # scienceuserfriendly package for GH
    # GH.test1 <- posthocTGH(y=Profile[,i], x=groups_ordered, method=c("games-howell"),
    #                       conf.level = 0.95, digits=2, p.adjust="none",
    #                       formatPvalue = TRUE) #$output$games.howell
    # GH option 2 
    #   GH.test <- posthoc_anova(Profile[,i] ~ groups_ordered,
    #                            method = "Games-Howell")
    #   GH_matrix <- GH.test$output$result
    #   sign.pairs.GH <- GH_matrix[which(GH_matrix$p < GH.alpha),]
    #   nb.sign.pairs.GH <- length(which(GH_matrix$p < GH.alpha))
    #   if (nb.sign.pairs.GH != 0){
    #     sign.pairs.GH <- cbind(sign.pairs.GH, trait = names(Profile)[i])
    #     GH.table <- rbind(GH.table, sign.pairs.GH)
    #     
    #     sign.pairs.GH_letters <- cbind(make_cld(GH.test),trait = names(Profile)[i])
    #     GH.table_letters <- rbind(GH.table_letters, sign.pairs.GH_letters)
    #   } else {
    #     sign.pairs.GH <- 0
    #     GH.table <- GH.table
    #   }
    # }
    
    # Dunn test
    d.test <- dunnTest(Profile[,i], g=groups_ordered, method="bh")
    # https://www.researchgate.net/post/How_do_I_choose_the_correct_p-adjustment_method_for_Dunns_Test
    #  Koenraad stelt BH voor ipv bonferonni = false discovery rate
    sign.pairs.bh <- d.test$res[which(d.test$res$P.adj < dunn.alpha),]
    nb.sign.pairs.bh <- length(which(d.test$res$P.adj < dunn.alpha))
    if (nb.sign.pairs.bh != 0){
      sign.pairs.bh <- cbind(sign.pairs.bh, trait = names(Profile)[i])
      dunn.table.bh <- rbind(dunn.table.bh, sign.pairs.bh)
      
      letters <- make_cld(P.adj ~ Comparison , data=d.test$res, alpha = 0.05)
      sign.pairs.bh_letters <- cbind(letters, trait = names(Profile)[i])
      dunn.table.bh_letters <- rbind(dunn.table.bh_letters, sign.pairs.bh_letters)
    } else {
      sign.pairs.bh <- 0
      dunn.table.bh <- dunn.table.bh
      dunn.table.bh_letters <- dunn.table.bh_letters
    }
    
    sign.pairs <- d.test$res[which(d.test$res$P.unadj < dunn.alpha),]
    nb.sign.pairs <- length(which(d.test$res$P.unadj < dunn.alpha))
    if (nb.sign.pairs != 0){
      sign.pairs <- cbind(sign.pairs, trait = names(Profile)[i])
      dunn.table <- rbind(dunn.table, sign.pairs)
      
      letters <- make_cld(P.unadj ~ Comparison , data=d.test$res, alpha = 0.05)
      sign.pairs_letters <- cbind(letters, trait = names(Profile)[i])
      dunn.table_letters <- rbind(dunn.table_letters, sign.pairs_letters)
    } else {
      sign.pairs <- 0
      dunn.table <- dunn.table
      dunn.table_letters <- dunn.table_letters
    }
  }
  
  # GH
  # for (i in 1:nrow(GH.table)) {
  #   if (GH.table[i,"p"] < 0.001) {
  #     GH.table[i,"p"] <- c("< 0.001")
  #   } else {
  #     GH.table[i,"p"] <- round(as.numeric(GH.table[i,"p"]),3)
  #   }}
  # 
  # GH.table$groups <- paste(rep("groups",dim(GH.table)[1]),GH.table$groups)
  # write.table(GH.table, file=paste("GH_table_",".csv",sep=name), sep =";", col.names=T, row.names = F)
  # write.table(GH.table_letters, file=paste("GH_table_letters_",".csv",sep=name), sep =";", col.names=T, row.names = F)
  # # return(GH.table)
  
  # Dunn
  dunn.table.ordered <- dunn.table[order(dunn.table[,'trait'],dunn.table[,"Comparison"]),]
  dunn.table.ordered$Comparison <- paste(rep("groups",dim(dunn.table.ordered)[1]),dunn.table.ordered$Comparison)
  dunn.table.ordered$Z <- round(dunn.table.ordered$Z, 2)
  # dunn.table.ordered$P.unadj <- round(dunn.table.ordered$P.unadj,3)
  for (i in 1:nrow(dunn.table.ordered)) {
    if (dunn.table.ordered[i,"P.unadj"] < 0.001) {
      dunn.table.ordered[i,"P.unadj"] <- c("< 0.001")
    } else {
      dunn.table.ordered[i,"P.unadj"] <- round(as.numeric(dunn.table.ordered[i,"P.unadj"]),3)
    }}
  
  dunn.table.bh.ordered <- dunn.table.bh[order(dunn.table.bh[,'trait'],dunn.table.bh[,"Comparison"]),]
  dunn.table.bh.ordered$Comparison <- paste(rep("groups",dim(dunn.table.bh.ordered)[1]),dunn.table.bh.ordered$Comparison)
  dunn.table.bh.ordered$Z <- round(dunn.table.bh.ordered$Z, 2)
  for (i in 1:nrow(dunn.table.bh.ordered)) {
    if (dunn.table.bh.ordered[i,"P.adj"] < 0.001) {
      dunn.table.bh.ordered[i,"P.adj"] <- c("< 0.001")
    } else {
      dunn.table.bh.ordered[i,"P.adj"] <- round(as.numeric(dunn.table.bh.ordered[i,"P.adj"]),3)
    }}
  for (i in 1:nrow(dunn.table.bh.ordered)) {
    if (dunn.table.bh.ordered[i,"P.unadj"] < 0.001) {
      dunn.table.bh.ordered[i,"P.unadj"] <- c("< 0.001")
    } else {
      dunn.table.bh.ordered[i,"P.unadj"] <- round(as.numeric(dunn.table.bh.ordered[i,"P.unadj"]),3)
    }}
  write.table(dunn.table.ordered, file=paste("Dunn_table_",".csv",sep=name), sep =";", col.names=T, row.names = F) 
  write.table(dunn.table.bh.ordered, file=paste("Dunn_table_BH_",".csv",sep=name), sep =";", col.names=T, row.names = F) 
  
  # ---------------
  # we now go over to visualizing the results
  
  # 1) For clarity we will transfor LA 
  if (trans == "LA"){
    Profile$LA <- log10(Profile$LA)
    Profile.long[which(Profile.long$variable == "LA"),"value"] <- log10(Profile.long[which(Profile.long$variable == "LA"),]$value)
  } else if (trans == "LA_log"){
    Profile$LA_log <- log10(Profile$LA_log)
    Profile.long[which(Profile.long$variable == "LA"),"value"] <- log10(Profile.long[which(Profile.long$variable == "LA_log"),]$value)
  }
  levels(Profile.long$variable) <- facet_names
  levels(dunn.table.bh_letters$trait) <- facet_names
  
  # 2) letters from make_cld are the labels we need (cld = compact letter display)
  # Determine positions for labels
  abs_max <- apply(Profile[,c(1:ncol(dataset))], 2, max)
  maxs <- Profile.long %>%
    group_by(groups_ordered,variable) %>%
    summarise(max=max(value)) %>%
    mutate(max2=max + 0.05 * abs_max) %>%
    mutate(value= abs_max + 0.1 * abs_max) %>%
    rename ("group" = "groups_ordered") %>%
    rename ("trait" = "variable") %>%
    arrange(factor(trait, levels=facet_names)) %>%
    left_join(dunn.table.bh_letters, by=c("group","trait"))
  maxs$groups_ordered2 = maxs$group
  maxs$groups_ordered2 <- paste(rep(text,nrow(maxs)), as.factor(maxs$groups_ordered2), sep=" ")
  
  #3) Plotting! Make sure variable names correspond
  Profile.long$trait <- Profile.long$variable
  p<-ggplot(Profile.long)+ #,fill=groups
    geom_boxplot(aes(forcats::fct_rev(factor(groups_ordered2)),value, color=groups_ordered))+
    geom_jitter(aes(forcats::fct_rev(factor(groups_ordered2)),value, color=groups_ordered))+
    scale_colour_manual(values=rev(colors_ordered))+
    geom_text(data=maxs, aes(x = groups_ordered2, y = value, label=cld), hjust="inward", color = "gray32")+
    coord_flip()+
    facet_wrap(~trait,scales="free", ncol = 2)+
    theme_bw()+
    theme(legend.position="none", axis.title.x = element_blank(), axis.title.y = element_blank())
  x11(width=720,height=720)
  print(p)
  ggsave(paste("Profile_traits_",".jpeg",sep=name),width = 16,height = 20,units=("cm"),dpi=600)
  
  # ---------------------
  # library(ggpubr) # ggboxplot
  # http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/
  
  # Each group against the mean
  # dmean <- Profile.long %>%
  #          group_by(variable) %>%
  #          summarise(MN = mean(value))
  # q <- ggboxplot(Profile.long, x="groups_ordered2",y="value",
  #                add="jitter", legend="none", color="groups_ordered2") +
  #   geom_hline(data=dmean, aes(yintercept = MN), linetype = 2)+ # Add horizontal line at base mean
  #   facet_wrap(~variable,scales="free", ncol = 2)+
  #   scale_colour_manual(values=rev(colors_ordered))+
  #   coord_flip()+
  #   stat_compare_means(method="t.test", label="p.signif", ref.group=".all.", hide.ns=T) + # Add pairwise comparisons p-value
  #   theme_bw() +
  #   theme(legend.position="none",axis.title.x = element_blank(), axis.title.y = element_blank()) 
  #   # x11(width=720,height=720)
  # print(q)
  # ggsave(paste(name,"pairwisemean.jpeg",sep=clusterdecision),width = 20,height = 20,units=("cm"),dpi=600)
  
  # Each group against each other
  # no_groups = 9
  # my_comparisons <- (combn(x = no_groups, m = 2, simplify=F))
  # q <- ggboxplot(Profile.long, x="groups_ordered",y="value",
  #                add="jitter", legend="none", color="groups_ordered") +
  #   scale_colour_manual(values=rev(colors_ordered))+
  #   facet_wrap(~variable,scales="free", ncol = 2)+
  #   stat_compare_means(comparisons = my_comparisons, method="wilcox.test", label="p.signif", hide.ns=T) # Add pairwise comparisons p-value
  # x11(width=720,height=720)
  # print(q)
  
}
Fig_DensityTraits <- function(dataset, groups, clusterdecision, name, colorscheme, group_order_nb,trans){
  # dataset.std<-data.stand(dataset,method='standardize',margin='column',plot=F) #test<-scale(dataset)
  diff.cat <- colnames(dataset)
  
  Profile<-as.data.frame(cbind(groups,dataset)) 
  if (trans == "LA"){
    Profile$LA <- log10(Profile$LA)
  } else if (trans == "LA_log"){
    Profile$LA_log <- log10(Profile$LA_log)
  }
  
  
  Profile.long<-melt(Profile,id.vars=c("groups"))
  # Profile.long$variable = with(Profile.long, factor(variable, levels = rev(levels(variable))))
  Profile.long$groups <- as.factor(Profile.long$groups)
  # Profile.long$groups = with(Profile.long, factor(groups, levels = rev(levels(groups))))
  
  group_order_original <- seq(1:length(group_order_nb))
  Profile$groups_ordered <- Profile$groups
  Profile.long$groups_ordered <- Profile.long$groups
  for (i in 1:length(group_order_original)){
    Profile.long[,"groups_ordered"] <- replace( Profile.long[,"groups_ordered"], Profile.long$groups==group_order_original[i],group_order_nb[i])
    Profile[,"groups_ordered"] <- replace(Profile[,"groups_ordered"], Profile$groups==group_order_original[i],group_order_nb[i])
  }
  
  colors <- create_color_list(colorscheme)
  colors<- as.data.frame(colors)
  colors$order_original <- seq(1:dim(colors)[1])
  colors <- colors[1:length(group_order_nb),]
  colors$colors_ordered <- colors$colors
  for (i in 1:length(group_order_original)){
    colors[group_order_nb[i],"colors_ordered"] <- colors[group_order_original[i],"colors"] 
  }
  colors_ordered <- as.character(colors$colors_ordered)
  
  p<-ggplot(Profile.long, aes(value, colour=groups_ordered,fill=groups_ordered))+ #,fill=groups
    geom_density(adjust=2)+
    scale_fill_manual(values=alpha(colors_ordered,0.3))+
    scale_colour_manual(values=colors_ordered)+
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~variable,scales="free", ncol = 2)+
    theme_classic()+
    theme(strip.text.x = element_blank())+
    theme(panel.spacing.y = unit(6, "lines"))+
    theme(panel.spacing.x = unit(2, "lines"))+
    guides(color = guide_legend(override.aes = list(size=8)))+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text=element_text(size=18),legend.key = element_blank(), plot.margin = unit(c(0.4,0.4,0.4,0.4), units = "cm"),
          legend.title = element_text(size = 20), legend.text = element_text(size = 18),legend.justification = "top")
  # x11(height=340,width=300)
  print(p)
  ggsave(paste(name,"_density.jpeg",sep=clusterdecision),width = 34,height = 30,units=("cm"),dpi=600)
  
  p<-ggplot(Profile.long, aes(value, colour=groups_ordered,fill=groups_ordered))+ #,fill=groups
    geom_density(adjust=2)+
    scale_fill_manual(values=alpha(colors_ordered,0.3))+
    scale_colour_manual(values=colors_ordered)+
    facet_wrap(~variable,scales="free", ncol = 2)+
    theme_bw()+
    # theme(strip.text.x = element_blank())+
    theme(panel.spacing = unit(2, "lines"))+
    guides(color = guide_legend(override.aes = list(size=8)))+
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.text=element_text(size=18),legend.key = element_blank(), plot.margin = unit(c(1,0.4,0.4,0.4), units = "cm"),
          legend.title = element_text(size = 20), legend.text = element_text(size = 18),legend.justification = "top")
  # x11(height=340,width=300)
  print(p)  
  ggsave(paste(name,"_density_withLabels.jpeg",sep=clusterdecision),width = 34,height = 30,units=("cm"),dpi=600)
}
Dunn_paired_refl_ordered <- function(dataset, groups, partitioning, group_order_nb){
  # dataset.std<-data.stand(dataset,method='standardize',margin='column',plot=F) #test<-scale(dataset)
  dataset.std <- dataset
  
  Profile<-as.data.frame(cbind(dataset.std,groups)) 

  group_order_original <- seq(1:length(group_order_nb))
  Profile$groups_ordered <- Profile$groups
  for (i in 1:length(group_order_original)){
    Profile[,"groups_ordered"] <- replace(Profile[,"groups_ordered"], Profile$groups==group_order_original[i],group_order_nb[i])
  }
  groups_ordered <- as.factor(Profile[,"groups_ordered"])
  

  dunn.alpha=0.06
  dunn.table.bonfer <- data.frame()
  dunn.table <- data.frame() 
  
  for (i in 1:dim(dataset.std)[2]) {
    d.test <- dunnTest(Profile[,i], g=groups_ordered,method="bonferroni")
    
    sign.pairs.bonfer <- d.test$res[which(d.test$res$P.adj < dunn.alpha),]
    nb.sign.pairs.bonfer <- length(which(d.test$res$P.adj < dunn.alpha))
    if (nb.sign.pairs.bonfer != 0){
      sign.pairs.bonfer <- cbind(sign.pairs.bonfer, wl = colnames(dataset.std)[i])
      dunn.table.bonfer <- rbind(dunn.table.bonfer, sign.pairs.bonfer)
    } else {
      sign.pairs.bonfer <- 0
      dunn.table.bonfer <- dunn.table.bonfer
    }
    
    sign.pairs <- d.test$res[which(d.test$res$P.unadj < dunn.alpha),]
    nb.sign.pairs <- length(which(d.test$res$P.unadj < dunn.alpha))
    if (nb.sign.pairs != 0){
      sign.pairs <- cbind(sign.pairs, wl = colnames(dataset.std)[i])
      dunn.table <- rbind(dunn.table, sign.pairs)
    } else {
      sign.pairs <- 0
      dunn.table <- dunn.table
    }
  }
  
  # dunn.alpha=0.05
  # dunn.table <- data.frame()
  # for (i in 1:dim(dataset.std)[2]) {
  #   d.test <- dunnTest(dataset.std[,i], g=groups,method="bonferroni")
  #   
  #   sign.pairs <- d.test$res[which(d.test$res$P.adj < dunn.alpha),]
  #   nb.sign.pairs <- length(which(d.test$res$P.adj < dunn.alpha))
  #   if (nb.sign.pairs != 0){
  #     sign.pairs <- cbind(sign.pairs, wl = names(dataset.std)[i])
  #     dunn.table <- rbind(dunn.table, sign.pairs)
  #   } else {
  #     sign.pairs <- 0
  #     dunn.table <- dunn.table
  #   }
  
  # dunn.table <- rbind(dunn.table, 
  #                     data.frame(wl = names(dataset.std)[i],nb.sign.pairs))
  
  # Report number of bands tested
  # cat(paste("Dunn test for ",names(dataset.std)[i]," ", i, "/", 
  #           dim(dataset.std)[2], "; p-value=", d.test$res$P.adj,"/n", sep=""))
  # }
  
  # for(i in 1:dim(dataset)[2]){
  #   dunnTest(dataset[,i], g=groups,method="bonferroni")
  #   dunnTest(dataset.std, g=groups,method="bonferroni")
  # }
  
  
  # Plot frequency of significantly differing pairs + add spectral signature of group 1
  wl <- colnames(dataset.std)
  group1_refl <- dataset.std[which(groups_ordered == 1),]
  group1_refl_mean <- apply(group1_refl,2,mean)
  group1_plot <- as.data.frame(cbind(wl=as.numeric(wl),refl=as.numeric(group1_refl_mean)))
  
  dunn.table$wl_cont <- as.numeric(as.character(dunn.table$wl))
  dunn.table.bonfer$wl_cont <- as.numeric(as.character(dunn.table.bonfer$wl))
  
  f <- ggplot(dunn.table,aes(x=wl_cont))+
    geom_histogram(binwidth = 1, col="gray")+
    scale_x_continuous(name="wavelength (nm)",limits = c(400,2400),breaks = seq(400,2400,200), expand=c(0,0))+
    geom_line(data=group1_plot, aes(x=wl, y=refl*15),size=0.8)+
    geom_rect(aes(xmin=1341,xmax=1459,ymin=0,ymax=12),fill="white")+
    geom_rect(aes(xmin=1781,xmax=1969,ymin=0,ymax=12),fill="white")+
    scale_y_continuous(name="number of significant pairs",breaks = seq(0,15,2),
                       sec.axis = sec_axis(~./15, name = "Reflectance"),
                       expand = c(0,0))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill="white"),
          axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))
  print(f)
  setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Dunn_test")
  ggsave(paste("Dunn_refl_","_ordered.jpeg",sep=partitioning),width = 20,height = 13,units=("cm"),dpi=600)
  
  
  f <- ggplot(dunn.table.bonfer,aes(x=wl_cont))+
    geom_histogram(binwidth = 1, col="gray")+
    scale_x_continuous(name="wavelength (nm)",limits = c(400,2400),breaks = seq(400,2400,200), expand=c(0,0))+
    geom_line(data=group1_plot, aes(x=wl, y=refl*15),size=0.8)+
    geom_rect(aes(xmin=1341,xmax=1459,ymin=0,ymax=12),fill="white")+
    geom_rect(aes(xmin=1781,xmax=1969,ymin=0,ymax=12),fill="white")+
    scale_y_continuous(name="number of significant pairs",breaks = seq(0,15,2),
                       sec.axis = sec_axis(~./15, name = "Reflectance"),
                       expand = c(0,0))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill="white"),
          axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))
  print(f)
  setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Dunn_test")
  ggsave(paste("Dunn_refl_","_ordered_bonfer.jpeg",sep=partitioning),width = 20,height = 13,units=("cm"),dpi=600)
  
  return(dunn.table.bonfer)
}



##### 6. Describe clusters: PERMANOVA: amount of functional variance explained by groups
# ------- code adapted from Thomas et al.(2019)
# Difference between MANOVA and PERMANOVA: https://rpubs.com/collnell/manova
# Thomas et al. (2018) first apply log10 and than normalisation
library(dplyr) # select
permanova <- function(dataset, groups, partitioning, nb_perm,traits_structure,traits_economics, trans){
  normalise<- function(x) (x-min(x))/(max(x)-min(x))
  data_traits<-dataset
  
  if (trans == "LA"){
    data_traits$LA <- log10(data_traits$LA)
  } else if (trans == "LA_log"){
    data_traits$LA_log <- log10(data_traits$LA_log)
  }
  
  data_traits <- normalise(data_traits) # scale between 0 and 1 --> seems necessary for the algorithm...
  data_traits_structure<-dplyr::select(data_traits, all_of(traits_structure))
  data_traits_economics<-dplyr::select(data_traits,all_of(traits_economics))

  groups<-as.factor(groups)
  
  permanova.rsq <- function(data, groups, traits, partitioning){
    ad.run <- adonis(data ~groups, permutations = nb_perm)
    rsq<-ad.run$aov.tab[1,5] #Extract R2
    pval<-ad.run$aov.tab[1,6] #Extract p-value (> F)
    rsq.table.axes<-rbind(rsq.table.axes,as.data.frame(cbind(traits,partitioning,rsq,pval))) #Combine
    return(rsq.table.axes)
  }
  rsq.table.axes <- NULL
  rsq.table.axes <- permanova.rsq(data_traits, groups, "all", partitioning)
  rsq.table.axes <- permanova.rsq(data_traits_structure, groups, "structure", partitioning)
  rsq.table.axes <- permanova.rsq(data_traits_economics, groups, "economics", partitioning)
  rsq.table.axes.output <- rsq.table.axes
  
  # --- Figure (corresponding to middle row in Fig.5 from Thomas et al. (2019))
  rsq.table.axes$rsq_perc <- 100 * as.numeric(as.character(rsq.table.axes$rsq))

  getPalette = colorRampPalette(brewer.pal(8, "Set2"))
  colors=getPalette(8)
  g <- groups.explained.fig<-ggplot(rsq.table.axes,aes(traits,rsq_perc))+
    geom_bar(stat = "identity",colour="black",fill=colors[1:3],width=0.8)+ #colors[1:4]
    coord_flip()+
    theme_bw()+
    ylim(0,100)+
    ggtitle("")+
    # theme(plot.title = element_text(hjust = 0.5))+
    # theme(plot.title = element_text(margin = margin(t = 10,r = 10,b = 10, l = 10)))+
    # theme(plot.title = element_text(size=15, vjust=1), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
    labs(y="Variance Explained (%)",x="Trait Type")
  print(g)  
  ggsave(paste("Permanova_",".jpeg",sep=partitioning),width = 13,height = 13,units=("cm"),dpi=600)
  
  rsq.table.axes.output$rsq<- round(as.numeric(as.character(rsq.table.axes.output$rsq)), 2)
  # rsq.table.axes.output$pval <- round(as.numeric(as.character(rsq.table.axes.output$pval)), 2)
  # write.table(rsq.table.axes.output[,1:4], file=paste("Permanova_",".csv",sep=partitioning), sep =";", col.names=T, row.names = F)
  
  return(rsq.table.axes.output)
}


##### 7. Describe clusters: PERMANOVA: amount of functional variance explained between varying amount and choice of traits explained by groups
# ------- code adapted from Thomas et al.(2019)
# Build permanova loops to calculate variance explained by groups, for all trait combinations
permanova_loop <- function(dataset, groups, partitioning, nb_perm){
  normalise<- function(x) (x-min(x))/(max(x)-min(x))
  data_traits<-dataset
  data_traits <- normalise(data_traits) # scale between 0 and 1 --> seems necessary for the algorithm...
  traits <- colnames(data_traits)
  
  groups<-as.factor(groups)
  data_traits$Group <- groups
  
  rsq.table=NULL #Output table
  for(i in 2:length(traits)){ #Start at two as won't work for just one trait
    combinations<-combn(traits, i) #extract unique combination of traits for i samples
    combinations<-rbind(combinations,"Group") #Add Group column
    
    for(j in 1:ncol(combinations)){ #For number of unique combinations
      trait.data<-data_traits[ , names(data_traits) %in% as.vector(combinations[,j])] #Extract trait columns
      trait.data <- subset(trait.data,complete.cases(trait.data[,c(1:i)])) #Extract only complete cases
      
        mydata<-trait.data[,-(i+1)] #extract numerical data
        
        #Conduct permanova
        ad.run<-adonis(mydata ~ groups, permutations=nb_perm) #permanova - how much variation explained by groups
        
        #Extract variables
        rsq<-ad.run$aov.tab[1,5] #Extract R2
        nvar<-i #Note number of traits
        vars<-paste(unlist(colnames(mydata)), sep=",", collapse=", ") #Note trait identities
        
        #Combine in table
        rsq.table<-rbind(rsq.table,as.data.frame(cbind(nvar,rsq,vars))) #Combine
  
      }
  }
  
  #Convert to numbers
  rsq.table$rsq<-as.numeric(as.character(rsq.table$rsq))
  se <- function(x) sd(x)/sqrt(length(x))
  rsq.table<-rsq.table %>%
    group_by(nvar,vars) %>%
    summarise (mean = mean(rsq),
               sd = sd(rsq),
               se = se(rsq))
  rsq.table$vars<-as.character(rsq.table$vars)
  
  rsq.table$w.LeafN<-apply(rsq.table,1,function(x) ifelse(grepl("Leaf N",x[2]),1,0))
  rsq.table$w.SLA<-apply(rsq.table,1,function(x) ifelse(grepl("SLA",x[2]),1,0))
  rsq.table$w.LDMC<-apply(rsq.table,1,function(x) ifelse(grepl("LDMC",x[2]),1,0)) 
  rsq.table$w.Chltotal<-apply(rsq.table,1,function(x) ifelse(grepl("Chl",x[2]),1,0))
  rsq.table$w.Height<-apply(rsq.table,1,function(x) ifelse(grepl("Plant Height",x[2]),1,0))
  rsq.table$w.LA<-apply(rsq.table,1,function(x) ifelse(grepl("area",x[2]),1,0))
  
  rsq.table$Type<-apply(rsq.table,1,function(x) ifelse(sum(c(as.numeric(x["w.LeafN"]),as.numeric(x["w.SLA"]),as.numeric(x["w.LDMC"]),as.numeric(x["w.Chltotal"])))==0,"structure",
                                                ifelse(sum(c(as.numeric(x["w.LeafN"]),as.numeric(x["w.SLA"]),as.numeric(x["w.LDMC"])))==0 & as.numeric(x["w.Chltotal"])==1,"structure with pigments",
                                                ifelse(sum(c(as.numeric(x["w.Height"]),as.numeric(x["w.LA"]),as.numeric(x["w.Chltotal"])))==0,"economics", # without pigments
                                                ifelse(sum(c(as.numeric(x["w.Height"]),as.numeric(x["w.LA"])))==0 & as.numeric(x["w.Chltotal"])==1,"economics with pigments",
                                                "combi")))))
  
  getPalette = colorRampPalette(brewer.pal(8, "Set2"))
  colors=getPalette(8)
  
  #Plot change in explanatory power of functional group with trait information
  a<-ggplot(rsq.table,aes(nvar,mean*100,colour=factor(Type), shape=factor(Type)))+
      geom_point(cex=2.5)+
      theme_bw()+
      scale_shape_manual(values=c(1,19,1,19,1),
                         name="Trait Combination",
                         breaks = c("combi","economics","economics with pigments","structure","structure with pigments"),
                         labels = c("Combination","Only LES traits","LES traits with pigments","Only Morphological traits","Morphological traits with pigments"))+
        scale_colour_manual(values=c("black",colors[1],colors[1],colors[2],colors[2]),
                          name="Trait Combination",
                          breaks = c("combi","economics","economics with pigments","structure","structure with pigments"),
                          labels = c("Combination","Only LES traits","LES traits with pigments","Only Morphological traits","Morphological traits with pigments"))+
      theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
      labs(x = "Number of traits", y = expression(paste("Variance Explained (%)")))+
      theme(axis.title.x = element_blank())
  print(a)
  ggsave(paste("Permanova_loop_LESvsMorpho_",".jpeg",sep=partitioning),width = 20,height = 13,units=("cm"),dpi=600)
  
  
  b<-ggplot(rsq.table,aes(nvar,mean*100,colour=factor(w.Chltotal),shape=factor(w.Chltotal)))+
    geom_point(cex=2.5)+
    theme_bw()+
    scale_shape_manual(values=c(1,19),
                       breaks = c(1,0),
                       name="Trait Combination",
                       labels=c("Includes pigments", "Excludes pigments"))+
    scale_colour_manual(values=c("black",colors[1]),
                        breaks = c(1,0),
                        name="Trait Combination",
                        labels=c("Includes pigments", "Excludes pigments"))+
    theme(legend.title=element_text(size=15) , legend.text=element_text(size=12), plot.title = element_text(size=18, vjust=1), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
    labs(x = "Number of traits", y = expression(paste("Variance Explained (%)")))
  print(b)
  ggsave(paste("Permanova_loop_Chl_",".jpeg",sep=partitioning),width = 20,height = 13,units=("cm"),dpi=600)
  
  
  
  c<-ggplot(rsq.table,aes(nvar,mean*100,colour=factor(w.LDMC),shape=factor(w.LDMC)))+
    geom_point(cex=2.5)+
    theme_bw()+
    scale_shape_manual(values=c(1,19),
                       breaks = c(1,0),
                       name="Trait Combination",
                       labels=c("Includes LDMC", "Excludes LDMC"))+
    scale_colour_manual(values=c("black",colors[2]),
                        breaks = c(1,0),
                        name="Trait Combination",
                        labels=c("Includes LDMC", "Excludes LDMC"))+
    theme(legend.title=element_text(size=15) , legend.text=element_text(size=12), plot.title = element_text(size=18, vjust=1), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), axis.line.x = element_line(color="black"),axis.line.y = element_line(color="black"), axis.text=element_text(size=12), axis.title=element_text(size=15,vjust = -0.5),axis.text.x=element_text(hjust = -0.05),legend.key = element_blank(), plot.margin = unit(c(1,1,1,1), units = , "cm"))+
    labs(x = "Number of traits", y = expression(paste("Variance Explained (%)")))
  print(c)
  ggsave(paste("Permanova_loop_LDMC_",".jpeg",sep=partitioning),width = 20,height = 13,units=("cm"),dpi=600)
  
}
 


#--------------------------------------------------------------------------------
# CLUSTERING CORRESPONDENCE
# Comparing 2 dendograms:
# https://cran.r-project.org/web/packages/dendextend/vignettes/introduction.html
# http://www.sthda.com/english/wiki/print.php?id=237

##### 1. Tanglegram and entanglement of two dendrograms
# See up

##### 2. Cophenetic correlationa and Baker's Gamma Index between dendrograms
library(tibble) # tibble:lst
library(corrplot) # cor.dendlist
dend_compare_metrics <- function(...){ # dend1_list, dend2_list,
  # https://cran.r-project.org/web/packages/dendextend/vignettes/introduction.html
  
  input_list <- tibble::lst(...)
  output_list <- lapply(X=input_list, function(x) {x[["dend_col"]]})
  dend_list <-do.call(dendlist, output_list)  

  
  ### 1. Global comparison: edge lengths do not make sense when different distance metrics used?
  all.equal(dend_list, use.edge.length = FALSE)

  ### 2. Correlation matrix
  # Check website for bootstrapping --> p-value and confidence interval
  
  # The cophenetic distance between two observations that have been clustered is defined to be
  # the inter-group dissimilarity at which the two observations are first combined into a single cluster. 
  # The cophenetic correlation (see sokal 1962) is the correlation between two cophenetic distance
  # matrices of two trees.
  # The value can range between -1 to 1. With near 0 values meaning that the two trees are not statistically similar. 
  cors_coph<-cor.dendlist(dend_list, method = "cophenetic")
  x11(width=720,height=720)
  corrplot(cors_coph, "pie", "lower",title="Cophenetic correlation",addCoef.col = "black")
  
  # Baker's Gamma Index (see baker's paper from 1974) is a measure of association (similarity) 
  # between two trees of Hierarchical clustering (dendrograms). It is defined as the 
  # rank correlation between the stages at which pairs of objects combine in each of the two trees.
  # Notice that this measure is not affected by the height of a branch but only of its
  # relative position compared with other branches.
  # The value can range between -1 to 1. With near 0 values meaning that the two trees are not statistically similar. 
  cors_baker<-cor.dendlist(dend_list, method = "baker") # idem to cor_bakers_gamma(dend_list)
  x11(width=720,height=720)
  corrplot(cors_baker, "pie", "lower",title="Baker' Gamma",addCoef.col = "black")
}

##### 3. Adjusted Rand Index and Fowlkes-Mallows Index between clusters
# Check the following link:
# https://scikit-learn.org/stable/modules/clustering.html#clustering-performance-evaluation
# !!! Wagner and Wagner (2007, page 16) and Meila (2003, page 7) provide strong arguments against FMI and ARI
library(dendextend) # FM_index
library(mclust) # adjustedRandIndex
# library(rcompanion) # cramerV
# library(clues) # adjustedRand

cluster_compare_metrics <- function(dataset,cluster_output, clusterdecision, cluster_output2, clusterdecision2, name){

  if (is.list(cluster_output) == T){
    if (clusterdecision == "DHC"){
      groups<-cluster_output$groups_DHC
    } else if (clusterdecision == "NbClust"){
      groups<-cluster_output$groups_NbClust
    } else if (clusterdecision == "nbclusters"){
      groups<-cluster_output$groups
    }
    Observations <- labels(cluster_output$distancematrix)
    names(groups) <- Observations
  } else {
    groups <- cluster_output
    names(groups) <- dataset$meta[,"ID"]
  }

    if (clusterdecision2 == "DHC"){
      groups2<-cluster_output2$groups_DHC
    } else if (clusterdecision2 == "NbClust"){
      groups2<-cluster_output2$groups_NbClust
    } else if (clusterdecision2 == "nbclusters"){
      groups2<-cluster_output2$groups
    }
    Observations2 <- labels(cluster_output2$distancematrix)
    names(groups2) <- Observations2    
  

  
  # (The observations of the two group vectors are sorted)

  ### 1.The Fowlkes-Mallows Index (see fowlkes 1983) (FM Index, or Bk) is a measure of similarity
  # between two clusterings. The FM index ranges from 0 to 1, a higher value indicates a
  # greater similarity between the two clusters.
  # FMI_sp <- FM_index(groups,dataset$meta[,c('Species.code')],assume_sorted_vectors = T)
  # FMI_type <- FM_index(groups, groups2, include_EV = T, assume_sorted_vectors = F) # vectors are sorted, so T or F gives the same result
  # # Output includes the attributes E_FM and V_FM for the relevant expectancy and variance 
  # under the null hypothesis of no-relation.
  # The E_FM and V_FM are the values expected under the null hypothesis that the two trees 
  # have the same topology but one is a random shuffle of the labels of the other 
  # (i.e.: "no connection" between the trees).
  
  ### 2. Adjusted Rand Index
    # https://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.htmlà
    # labelings that have pure clusters with members coming from the same classes but unnecessary splits are penalized
    # So, if different number of groupings, ARI will automatically be lower
    
  # ARI_sp <- adjustedRandIndex(groups,dataset$meta[,c('Species.code')])
  # ARI_type <- adjustedRandIndex(groups, groups2)

  # the "adjustedRand" function from "clues" package can calculate different metric simultaneously
  # the results correspond with those obtained with FM_index and adjustedRandIndex
  # FMI_type2 <-adjustedRand(groups, groups2, randMethod="FM") # this function assumes vectors are sorted.
  # FMI_type2 <-adjustedRand(groups, groups2) #  HA = Adjusted Rand Index
  
  ### 3. Cramer's V = correlation between two nominal/categorical variables: 
  # http://rcompanion.org/handbook/H_10.html
  # Cramer_sp <- cramerV(as.character(groups), as.character(dataset$meta[,c('Species.code')]))
  # Cramer_type <- cramerV(as.character(groups), as.character(groups2))
  
  # Table<-as.data.frame(cbind(FMI_sp, ARI_sp, Cramer_sp, FMI_type, ARI_type, Cramer_type))
  
  
  comparisons <- c("Species.code","Site","Date","Family","Order","Superorder","Group",
                   "Plant_growth_habit","Plant_lifespan_Pierce","Plant_lifespan_LEDA","Plant_growth_form_LEDA","CSR_Hunt")
  Table <- matrix(ncol = 2, nrow=length(comparisons)+1)
  rownames(Table) <- c(comparisons,"Emergent functional groups")
  # colnames(Table) <- c("ARI","FMI")
  colnames(Table) <- c("ARI","FMI")
  for (i in 1: length(comparisons)){
    Table[i,"ARI"] <- adjustedRandIndex(groups,dataset$meta[,comparisons[i]])
    Table[i,"FMI"] <- FM_index(groups,dataset$meta[,comparisons[i]],assume_sorted_vectors = T)
  }
  # Table["Emergent functional groups",] <- adjustedRand(groups, groups2, randMethod=c("HA","FM"))
  Table["Emergent functional groups","ARI"] <- adjustedRandIndex(groups,groups2)
  Table["Emergent functional groups","FMI"] <- FM_index(groups,groups2,assume_sorted_vectors = F)
  Table <- round(Table,2)
  
  write.table(Table, file=paste("Cluster_Correspondence_",".csv",sep=name), sep =";", col.names=NA) 
  
  return(Table)
}


##### 4. Procrustes analysis on 2 pca's
# By default PROTEST performs 999 permutations, we keep default settings
VI_pca_procrustes_ordered <- function (pca1, pca2, nb_comp, cluster1, cluster2, name, colorscheme1, group1_order_nb, group2_order_nb) {
  # Two ordinations can be very similar but this similarity may be masked as a result of the two ordinations
  # having different scalings, orientations and signs
  # the main difference between basic Procrustes and PROTEST is that protest is always symmetric whereas procrustes defaults to non-symmetric
  
  
  # pca2 will be rotated to match pca1. Arrows run thus from pca1 to pca2
  proc <- procrustes(pca1$x[,1:nb_comp], pca2$x[,1:nb_comp],symmetric=F)
  prot <- protest(pca1$x[,1:nb_comp], pca2$x[,1:nb_comp])
  
  # Simple plots
  win.graph()
  par(mfrow=c(1,2))
  plot(proc,kind="1")
  plot(proc,kind="2") # residual plot
  
  # Advanced plots: separate colour for each species (https://stackoverflow.com/questions/30325739/ggplot2-for-procrustes-rotation-in-vegan)
  # Procrustes residual plots
  
  # library(stringr) #str_sub
  # species <- str_sub(labels(pca1$x)[[1]],1,6)
  
  groups1 <- as.factor(cluster1$groups_DHC)
  groups1_order_original <- seq(1:nlevels(groups1))
  groups1_ordered <- groups1
  for (i in 1:length(group1_order_nb)){
    groups1_ordered <- replace(groups1_ordered, groups1 == groups1_order_original[i], group1_order_nb[i])
  }
  
  groups2 <- as.factor(cluster2$groups_DHC)
  groups2_order_original <- seq(1:nlevels(groups2))
  groups2_ordered <- groups2
  for (i in 1:length(group2_order_nb)){
    groups2_ordered <- replace(groups2_ordered, groups2 == groups2_order_original[i], group2_order_nb[i])
  }

  
  plotpro.sym <- data.frame(pca1x=prot$X[,1],
                            pca1y=prot$X[,2],
                            pca2x=prot$Yrot[,1],
                            pca2y=prot$Yrot[,2],
                            groups1=groups1_ordered, groups2=groups2_ordered)
  plotpro.asym <- data.frame(pca1x=proc$X[,1],
                             pca1y=proc$X[,2],
                             pca2x=proc$Yrot[,1],
                             pca2y=proc$Yrot[,2],
                             groups1=groups1_ordered, groups2=groups2_ordered)
  
  
  palette_color <- create_color_list(colorscheme1)
  palette_shape <- seq(0,14)
  
  colors<- as.data.frame(palette_color)
  colors$order_original <- seq(1:dim(colors)[1])
  colors <- colors[1:length(group1_order_nb),]
  colors$colors_ordered <- colors$palette_color
  for (i in 1:length(groups1_order_original)){
    colors[group1_order_nb[i],"colors_ordered"] <- colors[groups1_order_original[i],"palette_color"]
  }
  palette_color <- as.character(colors$colors_ordered)
 
  
  proc_plot <- function(plotpro_data,protest,sym){
    # x11(width=720,height=720)
    u <- ggplot(plotpro_data) +
      geom_point(aes(x=pca1x, y=pca1y, colour=groups1),size=3) +
      geom_point(aes(x=pca2x, y=pca2y, shape=groups2),size=3) +
      geom_segment(aes(x=pca1x,y=pca1y,xend=pca2x,yend=pca2y,colour=groups1),arrow=arrow(length=unit(0.2,"cm")),size=1)+
      scale_color_manual(name='ePFTs',values=palette_color) +
      scale_shape_manual(name='ePOTs',values=palette_shape)+
      annotate("text", x = Inf, y = -Inf,hjust=1.1,vjust=-3,size=5,
               label = paste("correlation",round(protest$t0,digits=2),sep=" = "), parse=F)+ # Add r and signifcance: find out how they are stored
      annotate("text", x = Inf, y = -Inf,hjust=1.1,vjust=-1,size=5,
               label = paste("p value",round(protest$signif,digits=3),sep=" = "), parse=F)+ # Add r and signifcance: find out how they are stored
      labs(x="PC1", y="PC2") +
      theme_bw()+ theme(text = element_text(size=16))+
      # theme(legend.key.size=unit(2,'lines'), legend.box="horizontal") +
      # theme(legend.position="none")+
      theme(legend.key.size=unit(2,'lines')) +
      theme(legend.box = "horizontal")+
      coord_equal()
    
    if (sym == "asym") {
      u <- u + scale_x_continuous(breaks=seq(floor(min(plotpro_data$pca1x,plotpro_data$pca2x))+floor(min(plotpro_data$pca1x,plotpro_data$pca2x))%%2,ceiling(max(plotpro_data$pca1x,plotpro_data$pca2x))+ceiling(max(plotpro_data$pca1x,plotpro_data$pca2x))%%2,2))+
        scale_y_continuous(breaks=seq(floor(min(plotpro_data$pca1y,plotpro_data$pca2y))+floor(min(plotpro_data$pca1y,plotpro_data$pca2y))%%2,ceiling(max(plotpro_data$pca1y,plotpro_data$pca2y))+ceiling(max(plotpro_data$pca1y,plotpro_data$pca2y))%%2,2))
    } 
    else if (sym == "sym") {
      u <- u + scale_x_continuous(breaks=seq(round(min(plotpro_data$pca1x,plotpro_data$pca2x),digits=1),max(round(plotpro_data$pca1x,plotpro_data$pca2x),digits=1),0.1))+
        scale_y_continuous(breaks=seq(round(min(plotpro_data$pca1y,plotpro_data$pca2y),digits=1),round(max(plotpro_data$pca1y,plotpro_data$pca2y),digits=1),0.1))
    }
    print(u)
  }
  
  proc_plot(plotpro.sym,prot,sym="sym")
  ggsave(paste("procrustes_sym_","_ordered.jpeg", sep=name),width = 23,height = 15,units=("cm"),dpi=600)
  proc_plot(plotpro.asym,prot,sym="asym")
  ggsave(paste("procrustes_asym_","_ordered.jpeg", sep=name),width = 23,height = 15,units=("cm"),dpi=600)
  
  # x11(width=720,height=720)
  # u <- ggplot(plotpro.sym) +
  #   geom_point(aes(x=pca1x, y=pca1y, colour=groups1),size=3) +
  #   geom_point(aes(x=pca2x, y=pca2y, shape=groups2),size=3) +
  #   geom_segment(aes(x=pca1x,y=pca1y,xend=pca2x,yend=pca2y,colour=groups1),arrow=arrow(length=unit(0.2,"cm")),size=1)+
  #   scale_color_manual(name='ePFTs',values=palette_color) +
  #   scale_shape_manual(name='ePOTs',values=palette_shape)+
  #   theme_bw()+ theme(text = element_text(size=16))+
  #   theme(legend.key.size=unit(2,'lines')) +
  #   theme(legend.box = "horizontal")+
  #   coord_equal()
  # library(cowplot)
  # legend<-get_legend(u)
  # ggdraw(plot_grid(legend,ncol=1))
  # setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Procrustes")
  # ggsave(paste("Procrustes_legend_","_ordered.jpeg",sep=name),width = 10,height = 15,units=("cm"),dpi=600)
}