#######################################################################################################
#######################################################################################################
#
# Code set #2
#
# Informing trait-based ecology by assessing remotely-sensed functional diversity across a 4 broad tropical temperature gradient
# SM Duran, RE Martin, S Díaz, BS Maitner, Y Malhi, N Salinas, A Shenkin, MR Silman, DJ Wieczynski, GP Asner, LP Bentley, VM Savage, BJ Enquist.
# This code has been designed for calculating functional diversity area curves from raster layers of plant chemical traits 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


library(geometry)

#Calculates a spatially contiguous fdar by randomly adding adjacent cells (observed_values)
#The spatial_contiguous_nested function was designed following Type I spatial sampling approach by Scheiner (2003, Ref.#51)


scn_fdar<-function(raster_brick,nreps){
  
  one_dim_ch<-function(raster_brick,occupied_cells_sampled){
    
    vals<-extract(x=raster_brick,y = occupied_cells_sampled)  
    return(abs(max(vals)-min(vals)))
  }
  
  
  output<-NULL
  
  valid_cells<-which(!is.na(getValues(raster_brick[[1]])))
  all_cells<-1:ncell(raster_brick)
  
  starting_cells<-NA
  
  if(length(names(raster_brick))==1){
    
    for(n in 1:nreps){
      
      cells_sampled<-NULL  
      occupied_cells_sampled<-NULL
      
      #Initiate by choosing a cell at random
      
      
      
      cells_sampled<-sample(x = all_cells[which(!all_cells%in%starting_cells)],size = 1)
      starting_cells<-c(starting_cells,cells_sampled)
      occupied_cells_sampled<-cells_sampled[which(cells_sampled%in%valid_cells)]
      
      
      #continue through until all cells sampled
      while(length(occupied_cells_sampled)!=length(valid_cells)){
        
        unsampled_cells<-all_cells[which(!all_cells%in%cells_sampled)]
        
        cells_sampled<-unique(adjacent(x = raster_brick,cells = cells_sampled,directions = 8,include = T)[,2])
        
        occupied_cells_sampled<-cells_sampled[which(cells_sampled%in%valid_cells)]
        print(paste("Replicate ",n," ",round(length(occupied_cells_sampled)/length(valid_cells)*100,digits = 2),"% complete"))
        
        
        if(length(occupied_cells_sampled)>length(names(raster_brick))){ #we need at least n+1 points to define a n dim volume
          
          #calc convex hull volume    
          
          convhull<-NA  
          #try(convhull<-convhulln(p = extract(x = raster_brick,y = occupied_cells_sampled),options = "FA")$vol,silent = T)
          try(convhull<-one_dim_ch(raster_brick = raster_brick,occupied_cells_sampled = occupied_cells_sampled),silent = T)
          
          #Calc other stuff
          
          #FDiv
          
          #traits_df<-extract(x = raster_brick,y=replicate(expr = cells_sampled[1],n = 10)) #Sanity check for FDiv
          traits_df<-extract(x = raster_brick,y = occupied_cells_sampled)
          centroid<-mean(traits_df)
          #centroid_distances<-apply(X = traits_df,MARGIN = 1,FUN = function(x){sum((x-centroid)^2)^.5})
          centroid_distances<-unlist(lapply(X = traits_df,FUN = function(x){abs(x-centroid)}))
          
          
          dGbar<-mean(centroid_distances)
          deltad <- sum( (1/length(centroid_distances))*abs(centroid_distances - dGbar) )
          FDiv <- dGbar / (deltad + dGbar)
          
          
          #calculate area (only including cells with values)
          area<-length(occupied_cells_sampled)
          
          #replicate
          replicate<-n
          
          #write iteration level output
          
          iter_output<-cbind(replicate,area,convhull,FDiv)
          output<-rbind(output,iter_output)
          
        }#if enough points to calc a volume  
        
        
        
      }#while statement  
      
      
    }#nreps loop
    
    
  }#if only one layer in the brick
  
  
  if(length(names(raster_brick))>1){
    
    for(n in 1:nreps){
      
      cells_sampled<-NULL  
      occupied_cells_sampled<-NULL
      
      #Initiate by choosing a cell at random
      
      
      
      cells_sampled<-sample(x = all_cells[which(!all_cells%in%starting_cells)],size = 1)
      starting_cells<-c(starting_cells,cells_sampled)
      occupied_cells_sampled<-cells_sampled[which(cells_sampled%in%valid_cells)]
      
      
      #continue through until all cells sampled
      while(length(occupied_cells_sampled)!=length(valid_cells)){
        
        unsampled_cells<-all_cells[which(!all_cells%in%cells_sampled)]
        
        cells_sampled<-unique(adjacent(x = raster_brick,cells = cells_sampled,directions = 8,include = T)[,2])
        
        occupied_cells_sampled<-cells_sampled[which(cells_sampled%in%valid_cells)]
        print(paste("Replicate ",n," ",round(length(occupied_cells_sampled)/length(valid_cells)*100,digits = 2),"% complete"))
        
        
        if(length(occupied_cells_sampled)>length(names(raster_brick))){ #we need at least n+1 points to define a n dim volume
          
          #calc convex hull volume    
          
          convhull<-NA  
          try(convhull<-convhulln(p = extract(x = raster_brick,y = occupied_cells_sampled),options = "FA")$vol,silent = T)
          
          #Calc other stuff
          
          #FDiv
          
          #traits_df<-extract(x = raster_brick,y=replicate(expr = cells_sampled[1],n = 10)) #Sanity check for FDiv
          traits_df<-extract(x = raster_brick,y = occupied_cells_sampled)
          centroid<-colMeans(traits_df)
          centroid_distances<-apply(X = traits_df,MARGIN = 1,FUN = function(x){sum((x-centroid)^2)^.5})
          dGbar<-mean(centroid_distances)
          deltad <- sum( (1/length(centroid_distances))*abs(centroid_distances - dGbar) )
          FDiv <- dGbar / (deltad + dGbar)
          
          #calculate area (only including cells with values)
          area<-length(occupied_cells_sampled)
          
          #replicate
          replicate<-n
          
          #write iteration level output
          
          iter_output<-cbind(replicate,area,convhull,FDiv)
          output<-rbind(output,iter_output)
          
        }#if enough points to calc a volume  
        
        
        
      }#while statement  
      
      
      
    }#nreps loop
    
    
  }#if more than one layer
  
  
  
  output<-as.data.frame(output)
  return(output)  
  
  
}
