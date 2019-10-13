
#######################################################################################################
#######################################################################################################
#
# Code set #1
#
# Informing trait-based ecology by assessing remotely-sensed functional diversity across a 4 broad tropical temperature gradient
# SM Duran, RE Martin, S Diaz, BS Maitner, Y Malhi, N Salinas, A Shenkin, MR Silman, DJ Wieczynski, GP Asner, LP Bentley, VM Savage, BJ Enquist.
# This code has been designed for calculating functional diversity area curves from raster layers of plant chemical traits 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#Load necessary packages 
library(raster)
library(doParallel)

#Load necessary functions, these functions are provided in the code folder
#The spatial_contiguous_nested function was designed following Type I spatial sampling approach by Scheiner (2003, Ref.#51)
source("code/spatial_contiguous_nested_FDAR.R")

#The non_spatial non_contiguous random function was designed following Type IIIB noncontiguous, nonspatial by Scheiner (2003, Ref.#51)
source("code/nonspatial_noncontiguous_random_FDAR.R")

#Load list of plots
plot_names<-strsplit(list.files("C:/Users/Owner/Dropbox/PDF/Carnegie/data/CAO_data/"),split = "_")
plot_names<-unlist(unique(lapply(X = plot_names,FUN = function(x){paste(x[1],x[2],sep ="_")})))

#Set parameters for FDARs
nreps=1000 #number of replicated FDARs to produce
parallelize<-TRUE

#Iterate through plots, calculating FDAR stuff
fdar_output <- NULL


for( i in 1:length(plot_names)){
  
  site<-plot_names[i]    
  
  #Load raster for site i
  #raster file name in this example is lidar_chem
  lidar_chem<-brick(paste("data/CAO_data3/",site,"_obs_lidar_chem",sep = ""))
  
   #Run FDAR code
  
  #First, scale chemical trait values in raster layer (lidar_chem) before calculating functional diversity indices
  lidar_chem<-scale(lidar_chem)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  #Run through spatial, contiguous, nested fdar fx (i.e. expanding window)
  if(!parallelize){
    scn<-scn_fdar(raster_brick = lidar_chem,nreps = nreps)
    scn<-cbind(site,"scn",scn)
  }
  
  
  if(parallelize){
    
    
    cl<-makePSOCKcluster(detectCores())
    doParallel::registerDoParallel(cl)
    system.time( scn<-foreach(i=nreps,.combine = rbind, .packages =  c('raster','geometry')) %dopar% (
      
      scn <- cbind(site,"scn",i,scn_fdar(raster_brick = lidar_chem,nreps = nreps))
      
      #scn <- cbind(site,"scn",scn)  
      
    )#dopar
    
    )#timing
    
    
    scn<-scn[,-which(colnames(scn)=='replicate')]  
    colnames(scn)[which(colnames(scn)=='i')]<-'replicate'
    colnames(scn)[which(colnames(scn)=='"scn"')]<-"fdar_type"
    stopCluster(cl)
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#  
  #Run through non-spatial, non-contiguous, random fdar (i.e. NULL)
  
  if(!parallelize){
    nnr<-nnr_fdar(raster_brick = lidar_chem,nreps = nreps,cells_to_add = 8)
    nnr<-cbind(site,"nnr",nnr)
  }
  
  if(parallelize){
    
    
    cl<-makePSOCKcluster(detectCores())
    doParallel::registerDoParallel(cl)
    system.time( nnr<-foreach(i=nreps,.combine = rbind, .packages =  c('raster','geometry')) %dopar% (
      
      nnr <- cbind(site,"nnr",i,nnr_fdar(raster_brick = lidar_chem,nreps = nreps,cells_to_add = 8))
      
      
    )#dopar
    
    )#timing
    
    
    nnr<-nnr[,-which(colnames(nnr)=='replicate')]  
    colnames(nnr)[which(colnames(nnr)=='i')]<-'replicate'
    colnames(nnr)[which(colnames(nnr)=='"nnr"')]<-"fdar_type"
    stopCluster(cl)
  }
  
  
  
  #Save output
  
  
  fdar_output<-rbind(fdar_output,rbind(scn,nnr))    
  #rm(scn,nnr)
  
}#end for i loop


#Cleanup
rm(i,plot_names,site,cl,nreps,parallelize)
colnames(fdar_output)[2]<-"fdar_type"


