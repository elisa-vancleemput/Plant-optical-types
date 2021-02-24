######################################################################
###                                                                ###
###     Multivariate trait estimation from hyperspectral data      ###
###                                                                ###
######################################################################

### The objective of this script is to estimate traits from hyperspectral data.
# We will use PLSR and assess how wel the estimates reflect conventionally measured traits.

### The input of these steps are the average spectra of the bouquets/patches of 
# - the unmixed spectra
# - with or without brightness normalization (_"br"): this is performed in "Load_PFT_spec"

# Script by Elisa Van Cleemput, 2017-20
#--------------------------------------------------------------------------------
# Clean workspace

rm(list=ls())

#--------------------------------------------------------------------------------
Load_PFT_spec <- function(proc,traitset,LUTdata, width) {
  
  #--------------------------------------------------------------------------------
  # 1) Conventionally measured functional traits
  
  PFT <- read.csv("Species_x_traits.csv",sep=";",header=TRUE)
  PFT <- PFT[order(PFT$ID),] # order the observations according to ID
  target <- which(names(PFT) == 'ID')[1]
  PFT_meta<-PFT[,1:target]
  PFT_traits<-PFT[,-c(1:target)]
  
  PFT_traits$LNC_mass <- 10 * PFT_traits$LNC # % dry mass to mg/g
  PFT_traits$Chltotal_mass <- PFT_traits$Chltotal
  
  PFT_traits$LNC_mass_log <- log10(PFT_traits$LNC_mass)
  PFT_traits$SLA_log <- log10(PFT_traits$SLA)
  PFT_traits$LDMC_log <- log10(PFT_traits$LDMC)
  PFT_traits$Chltotal_mass_log <- log10(PFT_traits$Chltotal_mass)
  PFT_traits$LA_log <- log10(PFT_traits$LA)
  PFT_traits$Height_log <- log10(PFT_traits$Height)
  
  PFT_traits$LNC_mass_sqrt <- sqrt(PFT_traits$LNC_mass)
  PFT_traits$SLA_sqrt <- sqrt(PFT_traits$SLA)
  PFT_traits$LDMC_sqrt <- sqrt(PFT_traits$LDMC)
  PFT_traits$Chltotal_mass_sqrt <- sqrt(PFT_traits$Chltotal_mass)
  PFT_traits$LA_sqrt <- sqrt(PFT_traits$LA)
  PFT_traits$Height_sqrt <- sqrt(PFT_traits$Height)

  #--------------------------------------------------------------------------------
  # 2) Reflectance
  
  spec <- read.csv(paste("Species_x_reflectance_",sep=proc,".csv"),sep=";",header=TRUE)
  target <- which(names(spec) == 'nm350')[1]
  spectra<-spec[,target:ncol(spec)]
  meta<-spec[,1:target-1]
  rownames(meta)=meta[,"ID"]
  
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
  
  # speclib_smooth<- smoothSpeclib(speclib,method="sgolay", n=51)
  
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
    colnames(spectra) <- wl
  }  else if (width == "1"){
    ### 3b) Removal of noise bands in case of no spectral binning
    mask(speclib_smooth)<-c(349,400,1340,1460,1780,1970,2400,2501) # for visualisation
    spectra <- speclib_smooth@spectra@spectra_ma # smoothed spectra in matrix
    wl<-speclib_smooth@wavelength
    colnames(spectra) <- wl
  }
  
  #### 4) Brightness normalization (optional) (only for raw and unmixed spectra, not 1deriv and 2deriv)
  if (isTRUE(grepl("br",LUTdata,fixed=T))){
    # Binned or detailed spectra
    brightness<-sqrt(apply(spectra^2,1,sum))
    spectra<-spectra/brightness
    speclib_smooth<-speclib(spectra,wl)
    SI(speclib_smooth)=meta_ordered
    mask(speclib_smooth)<-c(349,400,1340,1460,1780,1970,2400,2501) # for visualisation, this is binned in case of binning!
    colnames(spectra) <- wl
  } 
  
  output <- list ("PFT" = PFT,
                  "PFT_traits" = PFT_traits,
                  "PFT_meta" = PFT_meta,
                  "spectra" = spectra,
                  "spec_meta" = meta_ordered,
                  "wl" = wl)
  return(output)
}

setwd("C:/Users/u0091812/Box Sync/04. Optical types/Data and analyses/Raw data")
data_PFT_spec <- Load_PFT_spec(proc="unmixed", LUTdata = "_refl_sg", width = 1)
data_PFT_spec_br <- Load_PFT_spec(proc="unmixed", LUTdata = "_refl_sg_br", width = 1)

#--------------------------------------------------------------------------------
# Exploration: identify outliers in the dataset
#----------------------------
# http://r-statistics.co/Outlier-Treatment-With-R.html:
# For a given continuous variable, outliers are those observations that lie outside 1.5 * IQR,
# where IQR, the 'Inter Quartile Range' is the difference between 75th and 25th quartiles. 

traits <- c("LNC_mass","LNC_mass_log","LNC_mass_sqrt", 
            "SLA","SLA_log", "SLA_sqrt", 
            "LDMC", "LDMC_log", "LDMC_sqrt",
            "Chltotal_mass","Chltotal_mass_log","Chltotal_mass_sqrt",
            "LA","LA_log","LA_sqrt",
            "Height","Height_log","Height_sqrt")

my_outliers <- function(data, traits, name){
  PFT <- data$PFT_traits
  meta <- data$PFT_meta
  
  trait_outliers <- matrix(nrow=10, ncol=length(traits))
  colnames(trait_outliers) <- traits
  
  tiff(paste("Boxplot_traits_",".tiff",sep=name),res=300,width=15,height=8,units='in')
  # x11(width=1280,height=720)
  fig <- par(mfrow=c(3,4))
  
  for (i in 1:length(traits)){
    outlier_values <- boxplot.stats(PFT[,traits[i]])$out
    boxplot(PFT[,traits[i]], main=traits[i], boxwex=0.1)

    outlier_idx <- rownames(PFT[which(PFT[,traits[i]] %in% outlier_values),])
    outlier_ID <- meta[outlier_idx,"ID"]
    
    trait_outliers[,i]<-c(paste(as.character(outlier_ID), round(outlier_values,2), sep=": "),
                          rep(NA,10-length(outlier_ID)))
    
    # mtext(paste("Outliers: ", paste(outlier_ID, collapse=", ")), cex=0.6)
    # mtext(paste("Outliers: ", paste(outlier_values, collapse=", ")), cex=0.6)
  }
  dev.off()
  
  write.table(trait_outliers, file=paste("Outliers_traits_",".csv",sep=name), sep =";", row.names=F, col.names=T) 
  
  return(trait_outliers)
}

setwd("C:/Users/u0091812/Box Sync/04. Optical types/Data and analyses/PLSR results")
outliers_WithEQU <- my_outliers(data_PFT_spec, traits, name="WithEQU")

#--------------------------------------------------------------------------------
# PLSR : iterative stepwise procedure
#----------------------------
# PLSR step 1: determine optimal number of components (= latent variables, LV's)
  # "Cross-validation is commonly used to determine the optimal number of components to take into account"
  #   Mevik & Wehrens (2018); Salter (2018)
  # https://cran.r-project.org/web/packages/pls/vignettes/pls-manual.pdf
# PLSR step 2: Fit model and reduce number of predictor variables (= spectral bands)
# PLSR step 3: Fit model with reduced number of predictor variables, and again find optimal number of LV's
# Repeat steps 2 and 3 until RMSEP does not improve
# Execute this stepwise procedure several times (rep) to account for different data combinations in the CV procedure.
# Result: Caculate median predictive trait values across the repetitions

# By default, pls scales the data

# Input: - data = list containing "PFT_traits", "PFT_meta", "spectra", "spec_meta" and "wl"
#        - traits =  vector with trait names
#        - ncomp_traits_start = set a random initial number of no. of LVs to limit the extent of the plot (the same for each trait)
#        - rep = number of repetitions of step2 and 3
#        - wd = working directory to store the results
#        - name = string to store the results with a unique name


library(pls)
my_plsr <- function(data, trait, ncomp_traits_start, rep, wd, name)  {
  PFT <- data$PFT_traits
  rownames(PFT) <- data$spec_meta$ID
  spectra <- data$spectra
  rownames(spectra) <- data$spec_meta$ID
  wl <- data$wl
  
  # Create empty matrices that get to be filled with the plsr results
  traits_predicted_cv <- matrix(nrow = nrow(PFT), ncol = length(traits))
  colnames(traits_predicted_cv) <- traits
  rownames(traits_predicted_cv) <- rownames(PFT)
  traits_predicted_cal <- traits_predicted_cv
  
  traits_pred_acc <- matrix(nrow = 15, ncol = length(traits))
  rownames(traits_pred_acc) <- c("RMSEP_CV","nRMSEP_CV","RMSEP_CV_adj","nRMSEP_CV_adj",
                                 "RMSE_train","nRMSE_train",
                                 "R2_CV","R2_train",
                                 "TheilU_CV_1","TheilU_CV_2", "TheilU_train_1", "TheilU_train_2",
                                 "ncomp","nbands","nloops")
  colnames(traits_pred_acc) <- traits
  
  traits_regr_coeff_median <- matrix (nrow = length(traits), ncol =dim(spectra)[2])
  colnames(traits_regr_coeff_median) <- wl
  traits_regr_coeff_mean <- traits_regr_coeff_median
  traits_regr_coeff_min <- traits_regr_coeff_median
  traits_regr_coeff_max <- traits_regr_coeff_median
  traits_regr_coeff_nb <- traits_regr_coeff_median
  rownames(traits_regr_coeff_median) <- traits
  rownames(traits_regr_coeff_mean) <- traits
  rownames(traits_regr_coeff_min) <- traits
  rownames(traits_regr_coeff_max) <- traits
  rownames(traits_regr_coeff_nb) <- traits
  
  # -------------------------
  for (j in 1 : length(traits)){
# for (j in 1){
    
    # For every repetition: calculate predictions, statistic metrics and regression coefficients 
    # Create empty matrices for storing them
    trait_pred_acc <- matrix(nrow = 15, ncol = rep)
    rownames(trait_pred_acc) <- c("RMSEP_CV","nRMSEP_CV","RMSEP_CV_adj","nRMSEP_CV_adj",
                                   "RMSE_train","nRMSE_train",
                                   "R2_CV","R2_train",
                                   "TheilU_CV_1","TheilU_CV_2", "TheilU_train_1", "TheilU_train_2",
                                   "ncomp","nbands","nloops")
    trait_coeff <- matrix (nrow = rep, ncol =dim(spectra)[2])
    colnames(trait_coeff) <- wl
    trait_pred_cv <- matrix(nrow = nrow(PFT), ncol = rep)
    trait_pred_cal <- matrix(nrow = nrow(PFT), ncol = rep)
    
    for (i in 1 : rep){
  # for (i in 2){
        ### ----- create initial values of while-loop
      # 1) Start = model with all wavelenghts and all LV's
      # PLSR.trait_full<-plsr(trait ~ spectra, data = train, ncomp = ncomp_traits_start, validation="CV", jackknife = T)
      PLSR.trait_full<-plsr(PFT[,traits[j]]~spectra, ncomp = ncomp_traits_start, validation="CV", jackknife = T)
      spectra.sign.wl.update <- spectra
      
      # 2) Reduce no. of LV's before identifying teh number of significant bands
      ncomp.press.update <- ncomp.press <- which.min(PLSR.trait_full$validation$PRESS)
      # model with reduced no. of LV's = start data for looping
      PLSR.trait_LV<-plsr(PFT[,traits[j]] ~ spectra, ncomp = ncomp.press, validation="CV", jackknife = T)
      RMSEP_CV_new <- RMSEP(PLSR.trait_LV)$val[1,,ncomp.press+1]
      PLSR.trait <- PLSR.trait_LV
      
      # 3) Set other initial parameters
      RMSEP_CV_ini = Inf
      loop = 0  # record number of loops
      
      ### ----- Run while-loop
      while (RMSEP_CV_new < RMSEP_CV_ini && ncomp.press.update > 1 && dim(spectra.sign.wl.update)[2] != 0) {
        RMSEP_CV_ini <- RMSEP_CV_new
        
        PLSR.trait <- PLSR.trait_LV # this is the final model that we want to retain
        ncomp.press <- ncomp.press.update # this is the final number of components 
        spectra.sign.wl <- spectra.sign.wl.update # this are the final spectral bands
        
        #identify significant bands
        jack <- jack.test(PLSR.trait) 
        spectra.sign.wl.update <-spectra[,which(jack$pvalues < 0.05),drop=F]
        ncomp_traits <- ncomp_traits_start # let the no. of LV's vary freely again
        if (dim(spectra.sign.wl.update)[2] != 0){
          if (ncomp_traits > dim(spectra.sign.wl.update)[2]){
            ncomp_traits <- dim(spectra.sign.wl.update)[2]
          }
          # Fit model with only significant wavelengths and free no. of LV's
          PLSR.trait_wl<-plsr(PFT[,traits[j]] ~ spectra.sign.wl.update, ncomp = ncomp_traits, validation="CV", jackknife = T)
          # Again reduce no. of LV's
          ncomp.press.update <- which.min(PLSR.trait_wl$validation$PRESS)
          PLSR.trait_LV<-plsr(PFT[,traits[j]] ~ spectra.sign.wl.update, ncomp = ncomp.press.update, validation="CV", jackknife = T)
          
          RMSEP_CV_new <- RMSEP(PLSR.trait_LV)$val[1,,ncomp.press.update+1]
        }
        loop = loop+1
      }
      
      ### ----- Store statistics
      nRMSEP_CV     <- round(RMSEP(PLSR.trait)$val[1,,ncomp.press+1]/PLSR.trait$Ymean,2)
      # nRMSEP_CV     <- round(RMSEP(PLSR.trait)$val[1,,ncomp.press+1]/(max(PLSR.trait$fitted.values)-min(PLSR.trait$fitted.values)),2)
      # nRMSEP_CV     <- round(RMSEP(PLSR.trait)$val[1,,ncomp.press+1]/(max(PFT[,traits[i]])-min(PFT[,traits[i]])),2)
      nRMSEP_CV_adj     <- round(RMSEP(PLSR.trait)$val[2,,ncomp.press+1]/PLSR.trait$Ymean,2)
      # nRMSEP_CV_adj     <- round(RMSEP(PLSR.trait)$val[2,,ncomp.press+1]/(max(PFT[,traits[i]])-min(PFT[,traits[i]])),2)
      nRMSE_train  <- round(RMSEP(PLSR.trait, estimate="train")$val[1,,ncomp.press+1]/PLSR.trait$Ymean,2)
      # nRMSE_train  <- round(RMSEP(PLSR.trait, estimate="train")$val[1,,ncomp.press+1]/(max(PLSR.trait$fitted.values)-min(PLSR.trait$fitted.values)),2)
      # nRMSE_train  <- round(RMSEP(PLSR.trait, estimate="train")$val[1,,ncomp.press+1]/(max(PFT[,traits[i]])-min(PFT[,traits[i]])),2)
      trait_pred_acc["RMSEP_CV",i] <- round(RMSEP(PLSR.trait)$val[1,,ncomp.press+1],2)
      trait_pred_acc["nRMSEP_CV",i] <- nRMSEP_CV
      trait_pred_acc["RMSEP_CV_adj",i] <- round(RMSEP(PLSR.trait)$val[2,,ncomp.press+1],2)
      trait_pred_acc["nRMSEP_CV_adj",i] <- nRMSEP_CV_adj
      trait_pred_acc["RMSE_train",i] <- round(RMSEP(PLSR.trait, estimate="train")$val[1,,ncomp.press+1],2)
      trait_pred_acc["nRMSE_train",i] <- nRMSE_train
      trait_pred_acc["R2_CV",i] <-  round(R2(PLSR.trait)$val[1,,ncomp.press+1],2)
      trait_pred_acc["R2_train",i] <- round(R2(PLSR.trait, estimate="train")$val[1,,ncomp.press+1],2)
      
      # Extra's from Schweiger et al. 2017
      # ---- Theil's U
      obs_pred_cv <- as.data.frame(predplot(PLSR.trait, ncomp = ncomp.press, which="validation"))
      obs_pred_cal <- as.data.frame(predplot(PLSR.trait, ncomp = ncomp.press, which="train"))
      # (rm_int <- RMSEP (PLSR.trait, ncomp=ncomp_traits[i], intercept=T, estimate="CV"))
      # (rm_noint <- RMSEP (PLSR.trait, ncomp=ncomp_traits[i], intercept=F, estimate="CV")) # result also contained in rm_int
      # a <- sqrt(mean(obs_pred$measured^2))
      # b <- sqrt(mean(obs_pred$predicted^2))
      # (theil_u <- rm_int$val[[2]] /(a+b))
      # (theil_u <- rm_noint$val[[1]] /(a+b))
      library(DescTools) # TheilU
      # a <- obs_pred$measured
      # b <- obs_pred$predicted
      # TheilU(a, b, type = 1)
      # TheilU formula's on github:
      # n <- length(a)
      # res1 <- sqrt(sum((a-b)^2/n))/(sqrt(sum(a^2)/n) + sqrt(sum(b^2)/n)) # type 1
      # res2 <- sqrt(sum((a-b)^2))/(sqrt(sum(a^2))) # type 2
      trait_pred_acc["TheilU_CV_1",i] <- round(TheilU(obs_pred_cv$measured, obs_pred_cv$predicted, type = 1),2)
      trait_pred_acc["TheilU_CV_2",i] <- round(TheilU(obs_pred_cv$measured, obs_pred_cv$predicted, type = 2),2)
      trait_pred_acc["TheilU_train_1",i] <- round(TheilU(obs_pred_cal$measured, obs_pred_cal$predicted, type = 1),2)
      trait_pred_acc["TheilU_train_2",i] <- round(TheilU(obs_pred_cal$measured, obs_pred_cal$predicted, type = 2),2)
      
      trait_pred_acc["ncomp",i] <- ncomp.press
      trait_pred_acc["nbands",i] <- dim(spectra.sign.wl)[2]
      trait_pred_acc["nloops",i] <- loop
      
      
      ### ----- Store cross-validated fitted values
      # trait_pred_cv <- matrix(nrow = nrow(PFT), ncol = rep)
      trait_pred_cv[,i] <- as.data.frame(predplot(PLSR.trait, ncomp = ncomp.press, which="validation"))[,"predicted"]
      trait_pred_cal[,i] <- as.data.frame(predplot(PLSR.trait, ncomp = ncomp.press, which="train"))[,"predicted"]
      
      ### ----- Store regression coefficients
      # Difference between loadings and coefficients:
      #   - interpreting individual components: use loadings (or loading.weights) = regression coefficients of each latent variable
      #   - interpreting the whole model: use regression coefficients or VIP scores = combination of regression coefficients of latent variables and regression coefficients of overall model
      # Trait = a*(LV1)+b*(LV2)+c*(LV3)+d*(LV4)
      # LV1 = e*(wl1)+f*(wl2)+g*(wl3)+h*(wl4)+...
      regr_coeff <-  coef(PLSR.trait)
      trait_coeff[i, colnames(trait_coeff) %in% rownames(regr_coeff)] <- regr_coeff[rownames(regr_coeff) %in% colnames(trait_coeff), , ]
      }
    
    # Calculate average statistics and regression coefficients over all repetitions
    # Calculate median predictive values
    traits_pred_acc[,j] <- apply(trait_pred_acc, 1, mean)
    traits_predicted_cv[,j] <- apply(trait_pred_cv, 1, median)
    traits_predicted_cal[,j] <- apply(trait_pred_cal, 1, median)
    traits_regr_coeff_median[j,] <- apply(trait_coeff, 2, median, na.rm=T)
    traits_regr_coeff_mean[j,] <- apply(trait_coeff, 2, mean, na.rm=T)
    traits_regr_coeff_min[j,] <- apply(trait_coeff, 2, min, na.rm=T)
    traits_regr_coeff_max[j,] <- apply(trait_coeff, 2, max, na.rm=T)
    traits_regr_coeff_nb[j,] <- apply(trait_coeff, 2, function(x) length(which(!is.na(x))))
    }
 
  setwd(wd)
  traits_pred_acc <- t(traits_pred_acc)
  write.table(traits_pred_acc, file=paste(name,".csv",sep="Accuracy"), sep =";", row.names=T, col.names=NA) 
  write.table(traits_predicted_cv, file=paste(name,".csv",sep="predictions_cv"), sep =";", row.names=T, col.names=NA) 
  write.table(traits_predicted_cal, file=paste(name,".csv",sep="predictions_cal"), sep =";", row.names=T, col.names=NA) 
  write.table(traits_regr_coeff_median, file=paste(name,".csv",sep="regr_coeff_median"), sep =";", row.names=T, col.names=NA) 
  write.table(traits_regr_coeff_mean, file=paste(name,".csv",sep="regr_coeff_mean"), sep =";", row.names=T, col.names=NA) 
  write.table(traits_regr_coeff_min, file=paste(name,".csv",sep="regr_coeff_min"), sep =";", row.names=T, col.names=NA) 
  write.table(traits_regr_coeff_max, file=paste(name,".csv",sep="regr_coeff_max"), sep =";", row.names=T, col.names=NA) 
  write.table(traits_regr_coeff_nb, file=paste(name,".csv",sep="regr_coeff_nb"), sep =";", row.names=T, col.names=NA) 
  
     
  list <- list("traits_pred_acc" = traits_pred_acc,
               "traits_predicted_cv" = traits_predicted_cv,
               "traits_predicted_cal" = traits_predicted_cal,
               "traits_regr_coeff_median" = traits_regr_coeff_median,
               "traits_regr_coeff_mean" = traits_regr_coeff_mean,
               "traits_regr_coeff_min" = traits_regr_coeff_min,
               "traits_regr_coeff_max" = traits_regr_coeff_max,
               "traits_regr_coeff_nb" = traits_regr_coeff_nb)
  return(list)
  }
  
wd = "C:/Users/u0091812/Box Sync/04. Optical types/Data and analyses/PLSR results/a) PLSR models"
plsr_refl <- my_plsr(data_PFT_spec, traits, ncomp_traits_start=20, rep=20, wd=wd,  name = "PLSRmodel_refl_")
plsr_refl_br <- my_plsr(data_PFT_spec_br, traits, ncomp_traits_start=20, rep=20, wd=wd,  name = "PLSRmodel_refl_br_")

#--------------------------------------------------------------------------------
# Find best predictive model per trait: with or without brightness normalization?
#----------------------------
# For each trait, the model_best_results function compares the model statistics obtained by different plsr modelling approaches
# The output is an overview of the best model per trait, based on the different model statistics

# Input: - traits =  vector with trait names
#        - data = list containing different plsr approaches for each trait
#        - name = string to store the results with a unique name

model_best_results <- function(traits, data, wd, name){
    traits_pred_acc_best <- matrix(nrow = length(traits)*2, ncol = 12)
    traits_pred_acc_best <- as.data.frame(traits_pred_acc_best)
    colnames(traits_pred_acc_best) <- c("RMSEP_CV","nRMSEP_CV","RMSEP_CV_adj","nRMSEP_CV_adj","RMSE_train","nRMSE_train","R2_CV","R2_train","TheilU_CV_1","TheilU_CV_2", "TheilU_train_1", "TheilU_train_2")

  for (i in 1:length(traits)){
    for (j in 1:6){
      temp <- sapply(data, function(x) x[i,j])
      traits_pred_acc_best[i*2-1,j] <- min(temp)
      rownames(traits_pred_acc_best)[i*2-1] <- traits[i]
      rownames(traits_pred_acc_best)[i*2] <- paste("techn",traits[i],sep=".")
      if(length(labels(temp)[which.min(temp)]) == 0){
        traits_pred_acc_best[i*2,j] <- NA
      } else {
        traits_pred_acc_best[i*2,j] <- labels(temp)[which.min(temp)]
      }}
    for (j in 7:8){
      temp <- sapply(data, function(x) x[i,j])
      traits_pred_acc_best[i*2-1,j] <- max(temp)
      # colnames(traits_pred_acc_best[i,j*2-1]) <- traits[j]
      # colnames(traits_pred_acc_best[i,j*2]) <- paste(techn,traits[j],sep="_")
      if(length(labels(temp)[which.max(temp)]) == 0){
        traits_pred_acc_best[i*2,j] <- NA
      } else {
        traits_pred_acc_best[i*2,j] <- labels(temp)[which.max(temp)]
      }}
    for (j in 9:ncol(traits_pred_acc_best)){
      temp <- sapply(data, function(x) x[i,j])
      traits_pred_acc_best[i*2-1,j] <- min(temp)
      rownames(traits_pred_acc_best)[i*2-1] <- traits[i]
      rownames(traits_pred_acc_best)[i*2] <- paste("techn",traits[i],sep=".")
      if(length(labels(temp)[which.min(temp)]) == 0){
        traits_pred_acc_best[i*2,j] <- NA
      } else {
        traits_pred_acc_best[i*2,j] <- labels(temp)[which.min(temp)]
      }}
    }
  
  # Save overview
  setwd(wd)
  write.table(traits_pred_acc_best, file=paste(name,".csv",sep="best_models"), sep =";", row.names=T, col.names=NA) 
  
  return(traits_pred_acc_best)
}

plsr_list <- list("Unmixed" = plsr_refl$traits_pred_acc, 
                  "Unmixed_br" = plsr_refl_br$traits_pred_acc)
wd = "C:/Users/u0091812/Box Sync/04. Optical types/Data and analyses/PLSR results/b) PLSR best models"
PLSR_best_model_per_trait <- model_best_results(traits, plsr_list, wd, name="PLSR_")

#--------------------------------------------------------------------------------
# Store best predictions per trait
#----------------------------
# The function extract_best_results consults the previously stored model predictions and
# concatenates the results fom the best models in 1 table

# Within this function, the results generated through the my_plsr function are loaded --> manually defined)

# Input
# Inspect the table generated by the model_best_results function,
# and indicate for each trait which transformation of the trait and which plsr approach 
# you want to maintain for further ecological analyses

# ---------------------- Concatenate best predictions
extract_best_results <- function (traits, best_data, wd){
  setwd(wd)
  Pred_cal_unmixed <- read.csv("PLSRmodel_refl_predictions_cal.csv",sep=";",header=TRUE, row.names=1)
  Pred_cal_unmixed_br <- read.csv("PLSRmodel_refl_br_predictions_cal.csv",sep=";",header=TRUE, row.names=1)
  Pred_cv_unmixed <- read.csv("PLSRmodel_refl_predictions_cv.csv",sep=";",header=TRUE, row.names=1)
  Pred_cv_unmixed_br <- read.csv("PLSRmodel_refl_br_predictions_cv.csv",sep=";",header=TRUE, row.names=1)
  
  Pred_unmixed_stats <- read.csv("PLSRmodel_refl_Accuracy.csv",sep=";",header=TRUE, row.names=1)
  Pred_unmixed_br_stats <- read.csv("PLSRmodel_refl_br_Accuracy.csv",sep=";",header=TRUE, row.names=1)
  
  Pred_refl_coeff_median <- read.csv("PLSRmodel_refl_regr_coeff_median.csv",sep=";",header=TRUE, row.names=1)
  Pred_refl_coeff_mean <- read.csv("PLSRmodel_refl_regr_coeff_mean.csv",sep=";",header=TRUE, row.names=1)
  Pred_refl_coeff_min <- read.csv("PLSRmodel_refl_regr_coeff_min.csv",sep=";",header=TRUE, row.names=1)
  Pred_refl_coeff_max <- read.csv("PLSRmodel_refl_regr_coeff_max.csv",sep=";",header=TRUE, row.names=1)
  Pred_refl_coeff_nb <- read.csv("PLSRmodel_refl_regr_coeff_nb.csv",sep=";",header=TRUE, row.names=1)
  Pred_refl_br_coeff_median <- read.csv("PLSRmodel_refl_br_regr_coeff_median.csv",sep=";",header=TRUE, row.names=1)
  Pred_refl_br_coeff_mean <- read.csv("PLSRmodel_refl_br_regr_coeff_mean.csv",sep=";",header=TRUE, row.names=1)
  Pred_refl_br_coeff_min <- read.csv("PLSRmodel_refl_br_regr_coeff_min.csv",sep=";",header=TRUE, row.names=1)
  Pred_refl_br_coeff_max <- read.csv("PLSRmodel_refl_br_regr_coeff_max.csv",sep=";",header=TRUE, row.names=1)
  Pred_refl_br_coeff_nb <- read.csv("PLSRmodel_refl_br_regr_coeff_nb.csv",sep=";",header=TRUE, row.names=1)
  
  # Create empty matrices to fill
  best_predictions_cal <- matrix(nrow=nrow(Pred_cal_unmixed), ncol=length(traits))
  rownames(best_predictions_cal) <- rownames(Pred_cal_unmixed)
  colnames(best_predictions_cal) <- traits
  best_predictions_cv <- best_predictions_cal
  best_predictions_stats <- as.data.frame(matrix(nrow=length(traits), ncol = ncol(Pred_unmixed_stats)))
  rownames(best_predictions_stats) <- traits
  colnames(best_predictions_stats) <- colnames(Pred_unmixed_stats)
  
  best_predictions_regr_coeff_median <- data.frame()
  best_predictions_regr_coeff_mean <- data.frame()
  best_predictions_regr_coeff_min <- data.frame()
  best_predictions_regr_coeff_max <- data.frame()
  best_predictions_regr_coeff_nb <- data.frame()
  
  
  for (i in 1:length(traits)){
    if (best_data[i] == "unmixed"){
      best_predictions_cal[,i] <- Pred_cal_unmixed[,c(traits[i])]
      best_predictions_cv[,i] <- Pred_cv_unmixed[,c(traits[i])]
      best_predictions_stats[i,] <- Pred_unmixed_stats[c(traits[i]),]
      best_predictions_regr_coeff_median <- rbind(best_predictions_regr_coeff_median,Pred_refl_coeff_median[c(traits[i]),])
      best_predictions_regr_coeff_mean <- rbind(best_predictions_regr_coeff_mean,Pred_refl_coeff_mean[c(traits[i]),])
      best_predictions_regr_coeff_min <- rbind(best_predictions_regr_coeff_min,Pred_refl_coeff_min[c(traits[i]),])
      best_predictions_regr_coeff_max <- rbind(best_predictions_regr_coeff_max,Pred_refl_coeff_max[c(traits[i]),])
      best_predictions_regr_coeff_nb <- rbind(best_predictions_regr_coeff_nb,Pred_refl_coeff_nb[c(traits[i]),])
    } else if (best_data[i] == "unmixed_br"){
      best_predictions_cal[,i] <- Pred_cal_unmixed_br[,c(traits[i])]
      best_predictions_cv[,i] <- Pred_cv_unmixed_br[,c(traits[i])]
      best_predictions_stats[i,] <- Pred_unmixed_br_stats[c(traits[i]),]
      best_predictions_regr_coeff_median <- rbind(best_predictions_regr_coeff_median,Pred_refl_coeff_median[c(traits[i]),])
      best_predictions_regr_coeff_mean <- rbind(best_predictions_regr_coeff_mean,Pred_refl_coeff_mean[c(traits[i]),])
      best_predictions_regr_coeff_min <- rbind(best_predictions_regr_coeff_min,Pred_refl_coeff_min[c(traits[i]),])
      best_predictions_regr_coeff_max <- rbind(best_predictions_regr_coeff_max,Pred_refl_coeff_max[c(traits[i]),])
      best_predictions_regr_coeff_nb <- rbind(best_predictions_regr_coeff_nb,Pred_refl_coeff_nb[c(traits[i]),])
    }
  }
  
  best_predictions_stats <- cbind(best_predictions_stats,best_data)
  best_predictions_stats.print <- best_predictions_stats[,c("R2_train","R2_CV", "nRMSE_train", "nRMSEP_CV", "ncomp", "nbands", "nloops","best_data")]

  output <- list("best_predictions_cal" = best_predictions_cal,
                 "best_predictions_cv" = best_predictions_cv,
                 "best_predictions_stats.print" = best_predictions_stats.print,
                 "best_predictions_regr_coeff_median" = best_predictions_regr_coeff_median,
                 "best_predictions_regr_coeff_mean" = best_predictions_regr_coeff_mean,
                 "best_predictions_regr_coeff_min" = best_predictions_regr_coeff_min,
                 "best_predictions_regr_coeff_max" = best_predictions_regr_coeff_max,
                 "best_predictions_regr_coeff_nb" = best_predictions_regr_coeff_nb)
  
  return(output)
}

best_trait_transf <- c("SLA_log","LDMC_log","LNC_mass_log","Chltotal_mass_sqrt","LA_log","Height_sqrt")
best_spec_transf <- c("unmixed_br","unmixed","unmixed_br","unmixed_br","unmixed","unmixed_br")
wd = "C:/Users/u0091812/Box Sync/04. Optical types/Data and analyses/PLSR results/a) PLSR models"

best_pred <- extract_best_results(best_trait_transf, best_spec_transf, wd)

# Store the output
setwd("C:/Users/u0091812/Box Sync/04. Optical types/Data and analyses/PLSR results/b) PLSR best models")
write.table(best_pred$best_predictions_cal, file="PLSR_best_predictions_cal.csv", sep =";", row.names=T, col.names=NA)
write.table(best_pred$best_predictions_cv, file="PLSR_best_predictions_cv.csv", sep =";", row.names=T, col.names=NA)
write.table(best_pred$best_predictions_stats.print, file="PLSR_best_predictions_accuracy.csv", sep =";", row.names=T, col.names=NA)
write.table(best_pred$best_predictions_regr_coeff_median, file="PLSR_best_predictions_regr_coeff_median.csv", sep =";", row.names=T, col.names=NA)
write.table(best_pred$best_predictions_regr_coeff_mean, file="PLSR_best_predictions_regr_coeff_mean.csv", sep =";", row.names=T, col.names=NA)
write.table(best_pred$best_predictions_regr_coeff_min, file="PLSR_best_predictions_regr_coeff_min.csv", sep =";", row.names=T, col.names=NA)
write.table(best_pred$best_predictions_regr_coeff_max, file="PLSR_best_predictions_regr_coeff_max.csv", sep =";", row.names=T, col.names=NA)
write.table(best_pred$best_predictions_regr_coeff_nb, file="PLSR_best_predictions_regr_coeff_nb.csv", sep =";", row.names=T, col.names=NA)

#--------------------------------------------------------------------------------
# Visualise regression coefficients
#----------------------------

# It seemed that in the case of Leaf N, there are some regression coefficient outliers,
# and therefore the standardized figure does not show relative differences well.
# We will therefore replace them by the max/min values within the common range

library(reshape2)
library(ggplot2)
wl <- as.numeric(sub("X","",names(best_pred$best_predictions_regr_coeff_median),fixed = TRUE))
regr_coeff_median <- data.frame(cbind(wl, t(best_pred$best_predictions_regr_coeff_median)))
colnames(regr_coeff_median) <-c("wl","SLA","LDMC","Leaf N", "Leaf Chl a+b","Leaf area", "Height")

# Leaf N has 4 clear outliers! 3 positive and 1 negative
# Replace them by the 4th maximal and 2nd minimal value respectively
boxplot(regr_coeff_median$`Leaf N`)
library(Rfast)
max1_rownb <- Rfast::nth(regr_coeff_median$`Leaf N`, 1, descending = T, index.return = T)
max2_rownb <- Rfast::nth(regr_coeff_median$`Leaf N`, 2, descending = T, index.return = T)
max3_rownb <- Rfast::nth(regr_coeff_median$`Leaf N`, 3, descending = T, index.return = T)
max4 <- Rfast::nth(regr_coeff_median$`Leaf N`, 4, descending = T)
regr_coeff_median[c(max1_rownb,max2_rownb,max3_rownb),c("Leaf N")] <- max4
min1_rownb <- Rfast::nth(regr_coeff_median$`Leaf N`, 1, descending = F, index.return = T)
min2 <- Rfast::nth(regr_coeff_median$`Leaf N`, 2, descending = F)
regr_coeff_median[c(min1_rownb),c("Leaf N")] <- min2

regr_coeff_mean <- data.frame(cbind(wl, t(best_pred$best_predictions_regr_coeff_mean)))
colnames(regr_coeff_mean) <-c("wl","SLA","LDMC","Leaf N", "Leaf Chl a+b","Leaf area", "Height")
regr_coeff_nb <- data.frame(cbind(wl, t(best_pred$best_predictions_regr_coeff_nb)))
colnames(regr_coeff_nb) <-c("wl","SLA","LDMC","Leaf N", "Leaf Chl a+b","Leaf area", "Height")

stand <- function(x) {
  stand_factor <- max(abs(x))
  return(x/stand_factor)
}
regr_coeff_median_stand <- as.data.frame(apply(regr_coeff_median, 2, stand))
regr_coeff_median_stand$wl <- wl
regr_coeff_mean_stand <- as.data.frame(apply(regr_coeff_mean, 2, stand))
regr_coeff_mean_stand$wl <- wl

stats.melted <- melt(regr_coeff_median_stand, id.vars = "wl")
colnames(stats.melted)[3] <- "coeff_median"
stats.melted = cbind(stats.melted,
                     melt(regr_coeff_mean_stand, id.vars = "wl")[3],
                     melt(regr_coeff_nb, id.vars = "wl")[3]
                     )
colnames(stats.melted)[4:5] <- c("mean","nb")
stats.melted$refl_br <- apply(data_PFT_spec_br$spectra, 2, mean) # load data at line 21 if necessary
stats.melted$refl <- apply(data_PFT_spec$spectra, 2, mean)


# Without standardization:
  # median <- best_pred$best_predictions_regr_coeff_median
  # wl <- as.numeric(sub("X","",names(best_pred$best_predictions_regr_coeff_median),fixed = TRUE))
  # d <- data.frame(cbind(wl, t(best_pred$best_predictions_regr_coeff_median)))
  # stats.melted <- melt(median, id.vars = "wl")
  # stats.melted$value <- round(as.numeric(stats.melted$value), digits=2)
  # colnames(stats.melted)[3] <- "median"
  # stats.melted = cbind(stats.melted, 
  #                      round(melt(data.frame(cbind(wl,t(best_pred$best_predictions_regr_coeff_mean))), id.vars = "wl")[3], digits = 2),
  #                      round(melt(data.frame(cbind(wl,t(best_pred$best_predictions_regr_coeff_min))), id.vars = "wl")[3], digits = 2),
  #                      round(melt(data.frame(cbind(wl,t(best_pred$best_predictions_regr_coeff_max))), id.vars = "wl")[3], digits = 2),
  #                      melt(data.frame(cbind(wl,t(best_pred$best_predictions_regr_coeff_nb))), id.vars = "wl")[3])
  # colnames(stats.melted)[4:7] <- c("mean","min","max","nb")
  # stats.melted$refl_br <- apply(data_PFT_spec$spectra, 2, mean)
  # stats.melted$refl <- apply(data_PFT_spec_br$spectra, 2, mean)

setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/PLSR regression coefficients/")
x11(width=1000,height=720)
ggplot(stats.melted) + 
  geom_point(aes(x = wl, y = coeff_median, group=1, color = nb), size=0.5) +
  geom_point(aes(x = wl, y = refl*2, group=1), size=0.5, color = "gray50") +
  scale_y_continuous(sec.axis = sec_axis(~./2*100, name = "Reflectance [%]", breaks = c(0, 25, 50))) +
  geom_point(aes(x = wl, y = coeff_median, group=1, color = nb), size=0.5) +
  scale_colour_gradient(low = "thistle1",high = "orangered4")+
  ylab("standardized regression coefficients") +
  xlab("Wavelength (nm)") +
  facet_wrap(~ variable) +
  theme_bw()+
  theme(plot.margin = margin(10, 15, 10, 10))
ggsave("Regr_coefficients_median.jpeg",width = 25,height = 15,units=("cm"),dpi=600)

ggplot(stats.melted) +
  geom_point(aes(x = wl, y = mean, group=1, color = nb), size=0.5) +
  geom_point(aes(x = wl, y = refl*2, group=1), size=0.5, color = "gray50") +
  scale_y_continuous(sec.axis = sec_axis(~./2*100, name = "Reflectance [%]", breaks = c(0, 25, 50))) +
  geom_point(aes(x = wl, y = mean, group=1, color = nb), size=0.5) +
  scale_colour_gradient(low = "thistle1",high = "orangered4")+
  ylab("standardized regression coefficients") +
  xlab("Wavelength (nm)") +
  facet_wrap(~ variable) +
  theme_bw()+
  theme(plot.margin = margin(10, 15, 10, 10))
ggsave("Regr_coefficients_mean.jpeg",width = 25,height = 15,units=("cm"),dpi=600)

