# Run scripts:
# - 02_PFT_clustering.R
# - 03_PFT_clustering_RunFunctions.R: until line 54



library(hsdar)
library(dplyr) # filter

# Code extracted from the function "Fig_ProfileGroupsRefl_ordered" in 02_PFT_clustering.R

datalist <- data_PFT_spec
dataset <- datalist$spectra 
speclib <- speclib(dataset,as.numeric(colnames(dataset)))
SI(speclib) <- datalist$spec_meta
idSpeclib(speclib) <- as.character(datalist$spec_meta$ID) # contains 73 spectra
# Let's plot a set of randomly selected signatures
sign_nb <- c(2,8,16,20,35,44,51,59,64)
sign_nb <- c(2,8,20,35,64)
sign_nb <- c(2,9,21,36,65)

setwd("C:/Users/u0091812/Box Sync/04. Optical types/Figures/Conceptual figure")
tiff("Reflectance.tiff",width =18,height = 13,units=("cm"),res=600)
x11(width=2500,height=1000)
par(mar=c(5.1,5.1,2.1,2.1))
plot(c(),c(),axes=FALSE,ylim=c(0,0.6),xlim=c(350,2500),ylab='',xlab='', xaxs="i", yaxs="i")
for (i in 1:length(unique(sign_nb))){
  plot(speclib, FUN=i, lwd=3, ylim=c(0,1), col="#008B8B", new=FALSE)
  }
xx_mask1<-c(seq(1340,1460,1),rev(seq(1340,1460,1)))
xx_mask2<-c(seq(1780,1970,1),rev(seq(1780,1970,1)))
yy_mask1<-c(seq(0,0,length.out=length(xx_mask1)/2),seq(1,1,length.out=length(xx_mask1)/2))
yy_mask2<-c(seq(0,0,length.out=length(xx_mask2)/2),seq(1,1,length.out=length(xx_mask2)/2))
polygon(xx_mask1,yy_mask1,col="white",border=NA,new=F)
polygon(xx_mask2,yy_mask2,col="white",border=NA,new=F)
# abline(h=0, lwd=2)
axis(1,seq(350,2500,350),seq(350,2500,350),font=2,cex.axis=1.5)
axis(2,seq(0.0,0.6,0.2),seq(0.0,0.6,0.2),font=2,las=2,cex.axis=1.5)
mtext('Reflectance',side=2,font=2,line=3,cex=1.5)
mtext('Wavelength (nm)',side=1,font=2,line=3,cex=1.5)
dev.off()





# Fig_ProfileGroupsRefl_ordered(cluster_PFT, "DHC", data_PFT_spec, name = "ProfileGroupsRefl_PFT_", colorscheme = "Dark2",  
                              
cluster_output <- cluster_PFT
groups <- cluster_output$groups_DHC
datalist <- data_PFT_spec
group_order_nb=c(7,10,3,4,1,2,8,5,9,6) # this is not important because we will only use 1 color


Fig_ProfileGroupsRefl_ordered <- function(cluster_output,clusterdecision, datalist, name, colorscheme,group_order_nb) {

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


