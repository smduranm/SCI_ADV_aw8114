
#######################################################################################################
#######################################################################################################
#
# Code set #3
#
# Informing trait-based ecology by assessing remotely-sensed functional diversity across a 4 broad tropical temperature gradient
# SM Duran, RE Martin, S Díaz, BS Maitner, Y Malhi, N Salinas, A Shenkin, MR Silman, DJ Wieczynski, GP Asner, LP Bentley, VM Savage, BJ Enquist.
# This code has been designed for calculating functional diversity area curves from raster layers of plant chemical traits 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#Code to calculate a non-spatial FDAR (i.e. a null expectation)
#The non_spatial non_contiguous random function was designed following Type IIIB noncontiguous, nonspatial by Scheiner (2003, Ref.#51)


library(geometry)
library(raster)

#Calculates a spatially contiguous functional diversity area curves by randomly adding adjacent cells
#nreps = number of replicates (should be less than or equal to number of cells with values)
#raster_brick...the name of this has to be specified in the file FDAR_curves_from_raster.R before running
#cells_to_add: the number of cells to add at each timestep.  Higher number speeds process

nnr_fdar<-function(raster_brick,nreps,cells_to_add=1){

  
  one_dim_ch<-function(raster_brick,occupied_cells_sampled){
    
    vals<-extract(x=raster_brick,y = occupied_cells_sampled)  
    return(abs(max(vals)-min(vals)))
  }
  
    
  output<-NULL
  
  valid_cells<-which(!is.na(getValues(raster_brick[[1]])))
  
    if(length(names(raster_brick))>1){  
      
      for(n in 1:nreps){
        
        cells_sampled<-NULL  
        
        #Initiate by choosing a cell at random
        
        cells_sampled<-sample(x = valid_cells,size = cells_to_add)
        
        #continue through until all cells sampled
        while(length(cells_sampled)<length(valid_cells)){
          
          unsampled_cells<-valid_cells[which(!valid_cells%in%cells_sampled)]
          
          if(length(unsampled_cells) < cells_to_add){
          cells_sampled<-valid_cells  
            
          }else{
            
            cells_sampled<-c(cells_sampled,sample(x = unsampled_cells,size = cells_to_add,replace = F))   
          }
          
           
          
          
          
          print(paste("Replicate ",n," ",round(length(cells_sampled)/length(valid_cells)*100,digits = 2),"% complete"))
          
          if(length(cells_sampled)>length(names(raster_brick))){ #we need at least n+1 points to define a n dim volume
            
            #calc convex hull volume    
            
            convhull<-NA  
            try(convhull<-convhulln(p = extract(x = raster_brick,y = cells_sampled),options = "FA")$vol,silent = T)
            
            
            #Calc other stuff
            
            #FDiv
            
            #traits_df<-extract(x = raster_brick,y=replicate(expr = cells_sampled[1],n = 10)) #Sanity check for FDiv
            traits_df<-extract(x = raster_brick,y = cells_sampled)
            centroid<-colMeans(traits_df)
            centroid_distances<-apply(X = traits_df,MARGIN = 1,FUN = function(x){sum((x-centroid)^2)^.5})
            dGbar<-mean(centroid_distances)
            deltad <- sum( (1/length(centroid_distances))*abs(centroid_distances - dGbar) )
            FDiv <- dGbar / (deltad + dGbar)
            
             
            
            #calculate area (only including cells with values)
            area<-length(cells_sampled)
            
            #replicate
            replicate<-n
            
            #write iteration level output
            
            iter_output<-cbind(replicate,area,convhull,FDiv)
            
            output<-rbind(output,iter_output)
            
          }#if enough points to calc a volume  
          
          
          
          
          
          
        }#while statement  
        
        
        
        
        
        
      }#nreps loop
    
    }#if more than one raster layer  
  
  
    if(length(names(raster_brick))==1){
      
      for(n in 1:nreps){
        
        cells_sampled<-NULL  
        
        #Initiate by choosing a cell at random
        
        cells_sampled<-sample(x = valid_cells,size = cells_to_add)
        
        #continue through until all cells sampled
        while(length(cells_sampled)<length(valid_cells)){
          
          unsampled_cells<-valid_cells[which(!valid_cells%in%cells_sampled)]
          
          if(length(unsampled_cells) < cells_to_add){
            cells_sampled<-valid_cells  
            
          }else{
            
            cells_sampled<-c(cells_sampled,sample(x = unsampled_cells,size = cells_to_add,replace = F))   
          }
          
          
          
          
          print(paste("Replicate ",n," ",round(length(cells_sampled)/length(valid_cells)*100,digits = 2),"% complete"))
          
          if(length(cells_sampled)>length(names(raster_brick))){ #we need at least n+1 points to define a n dim volume
            
            #calc convex hull volume    
            
            convhull<-NA  
            #try(convhull<-convhulln(p = extract(x = raster_brick,y = cells_sampled),options = "FA")$vol,silent = T)
            try(convhull<-one_dim_ch(raster_brick = raster_brick,occupied_cells_sampled = cells_sampled),silent = T)
            
            #Calc other stuff
            
            #FDiv
            
            #traits_df<-extract(x = raster_brick,y=replicate(expr = cells_sampled[1],n = 10)) #Sanity check for FDiv
            traits_df<-extract(x = raster_brick,y = cells_sampled)
            #centroid<-colMeans(traits_df)
            centroid<-mean(traits_df)
            #centroid_distances<-apply(X = traits_df,MARGIN = 1,FUN = function(x){sum((x-centroid)^2)^.5})
            centroid_distances<-unlist(lapply(X = traits_df,FUN = function(x){abs(x-centroid)}))
            dGbar<-mean(centroid_distances)
            deltad <- sum( (1/length(centroid_distances))*abs(centroid_distances - dGbar) )
            FDiv <- dGbar / (deltad + dGbar)
            
           
            #calculate area (only including cells with values)
            area<-length(cells_sampled)
            
            #replicate
            replicate<-n
            
            #write iteration level output
            
            iter_output<-cbind(replicate,area,convhull,FDiv)
            
            output<-rbind(output,iter_output)
            
          }#if enough points to calc a volume  
          
          
          
        }#while statement  
        
        
        
      }#nreps loop
    
    }#if only one raster layer  
  
  output<-as.data.frame(output)
  return(output)
  
}
