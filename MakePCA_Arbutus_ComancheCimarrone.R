# Setup -------------------------------------
    
    library(terra)
    library(AFFT)
    library(snowfall)
    SDMap::MakeDirIfNeeded("ComancheCimarrone","C:/TEMP")
  ###functions----
    pred1<-function(x,y,keepvec = which(cumsum(x$sdev/sum(x$sdev))<.99)){z<-round(predict(x,y)[,keepvec] * 100,0)}

  ###boundaries -----
    bnd<-vect("Y:/MPSG_VegMapping/Data/Spatial/Mapping boundaries final/MPSG_mapping_boundaries_final.shp")
    bnd0<- bnd<-project(bnd,rast("Tiles/r_fp-0.25/r_fp-0.25_196.tif"))
    bnd<-bnd[bnd$FORESTNAME %in% c("Comanche and Cimarron National Grasslands"),]#
    bnd<-buffer(bnd,100000)
  ###nlcd ----- 
    nlcdmask<-rast("Y:/MPSG_VegMapping/Data/Raster/Source/nlcd/nlcd_2021_crop.tif")
    nlcdmask<-crop(nlcdmask,bnd)
    nlcdmask<-mask(nlcdmask,bnd)
    nlcdmask<-subset(nlcdmask,1)
    nlcdmask2<-nlcdmask %in% c("Unclassified","Barren Land","Deciduous Forest","Evergreen Forest","Mixed Forest","Shrub/Scrub","Herbaceous","Hay/Pasture","Woody Wetlands","Emergent Herbaceous Wetlands")
    nlcdmask2<-mask(nlcdmask2,nlcdmask2,maskvalue = F)
  ###point samples ----- 
    samp1<-spatSample(nlcdmask2,100000,"random",xy = T,na.rm = T)
    
    fl<-list.files("C:/TEMP/ComancheCimarrone/",pattern = ".vrt",full.names = T)
    nm<-list.files("C:/TEMP/ComancheCimarrone/",pattern = ".vrt",full.names = F)
    names(fl)<-gsub(".vrt","",nm)
    v2<-rast(fl)
    xy<-samp1[,1:2]
    suppl<-rast("D:/LocalNaip/Suppl/suppl.vrt")
    suppl<-crop(suppl,nlcdmask)
    nlcdmask<-crop(nlcdmask,suppl)
    suppl<-mask(suppl,nlcdmask)
  
# Artifacts ---------------
  if(!file.exists("C:/TEMP/ComancheCimarrone/Corrected/artifactscores.csv")){
    windows(10,10);plot(v2,1);plot(bnd0,add = T);e1<-draw();dev.off()
      par(mfrow =c(2,2))
    artifacts<-data.frame(do.call(rbind,pbapply::pblapply(names(v2),function(b){
      v1<-subset(v2,b)
      plot(v1,main = b); plot(bnd0,add = T)
      plotRGB(c(v1,v1,v1),stretch = "lin")
      plot(v1,ext = e1);plot(bnd0,add = T)
      plotRGB(c(v1,v1,v1),stretch = "lin",ext = e1)
      plot(bnd0,border = "red",add = T)
      cat("\n################\n ",b,"\n###############\n\n")
      return(c(Flightline = readline("Flightline score: "),Phenology = readline("Phenology Score: "), image = b))
    })))
    write.csv(artifacts,"C:/TEMP/ComancheCimarrone/Corrected/artifactscores.csv",row.names = F)
  }else{artifacts<-read.csv("C:/TEMP/ComancheCimarrone/Corrected/artifactscores.csv")}
  

  patchy<-artifacts$image[artifacts$Phenology %in% c(4,5)]
  stripey<-artifacts$image[artifacts$Flightline %in% c(4,5)]

# Generate vrt objects for all layers, tiles --------------
  # setwd("C:/TEMP/ComancheCimarrone")
  
  # fl<-list.files("merged",full.names = T,pattern = "tif")
  # cv<-gsub("merged/","",fl);cv<-reshape2::colsplit(cv,"_",c("one","two","three"))[,c(1:2)]
  # cv<-unique(paste(cv$one,cv$two,sep = "_"))
  # v2<-rast(fl)
  # names(v2)<-cv
  #tmp<-rast("N:/mpsg_naip_afft/bri/m_4810554_se_13_060_20210613.tif")
  # 
  # source("Y:/MPSG_VegMapping/elm_dev/CodeRepo/data_pipelines/raster_management/TileRasterFunctions.R")
  # g <- st_read("Y:/MPSG_VegMapping/Data/Spatial/spatial_database.gpkg", layer = "non-overlapping tiles")
  # pbapply::pblapply(names(v2),function(x){
  #   cat(x)
  #   v0<-subset(v2,x)
  #   op<-SDMap::MakeDirIfNeeded(x, "C:/TEMP/ComancheCimarrone/Tiles",goto = F)
  #   crt <- create_raster_tiles(v0, grid = g, outpath = op,variable_name = x)
  # })
  # vl<-pbapply::pbsapply(names(v2),function(x){
  #   op<-SDMap::MakeDirIfNeeded(x, "C:/TEMP/ComancheCimarrone/Tiles",goto = F)
  #   fl<-list.files(op,full.names = T)
  #   outvrt<-paste("C:/TEMP/ComancheCimarrone/",x,".vrt",sep = "")
  #   v<-vrt(fl,filename = outvrt,overwrite = T)
  #   outvrt
  #   
  # })
  vl<-list.files("C:/TEMP/ComancheCimarrone",pattern = "vrt",full.names = T)
  nm<-list.files("C:/TEMP/ComancheCimarrone",pattern = "vrt",full.names = F);nm<-gsub(".vrt","",nm)
  names(vl)<-nm
  
# Build princomp -------------
  ### Gather Info ---------------- 
    mat2<-extract(suppl,xy)
    
    sfInit(T,5)
    sfLibrary(terra)
    sfExport("xy")
    df1<-sfClusterApplyLB(vl,function(x){extract(rast(x),xy)})
    sfStop()
    
    df2<-lapply(df1,function(x){x[,2]})
    mat1<-data.frame(do.call(cbind,df2))
    colnames(mat1)<-names(sv2)
    
    mat3<-cbind(mat1,mat2[,2:ncol(mat2)])
    mat1<-mat3[apply(mat3,1,function(x){!any(is.na(x))}),]
    rm(list =c("mat2","mat3"))
  
  ### Build model -------------------
    mat3<-mat1
    drop.me<-list()
    
    drop.me[[1]]<-c("jdate","sl","imgyr","stateid",unique(c(patchy,stripey)))
 #   drop.me[[2]]<-c('ndng_f-3','ndvi_kurt','ndgr_f-0.25','r_f-0.25','g_f-0.25','ndvi_f-0.25')
 #   drop.me[[3]]<-c('r_med','ndgr_med','g_Q025','r_Q025','ndgr_Q025','ndng_Q025','ndng_med','bri_Q025','n_med','g_f-1.5')
    mat3<-mat1[,!colnames(mat1) %in% do.call(c,drop.me)]
    pca1<-princomp(mat3)
     
    save(pca1,file = "C:/TEMP/ComancheCimarrone/PCA/pca1.RData")

# Make iamges -------------------    

    tl<-substr(list.files("Tiles/bri_f-0.25"),12,14)
    fll<-list.files("Tiles",recursive = T, full.names = T,pattern = "tif")
    
    tids<-substr(reshape2::colsplit(fll,"_",letters[1:4])[,4],1,3)
    fllt<-split(fll,tids)
    
    sfInit(T,20)
    sfLibrary(AFFT)
    sfLibrary(terra)
    sfExport("pred1","pca1")
  
  
    fl<-sfClusterApplyLB(fllt,function(x){
      y<-rast(x)
      y<-subset(y,names(pca1$center))
      if(!any(is.na(y@ptr@.xData$range_max))){
        tid<-substr(x[1],nchar(x[1])-6,nchar(x[1])-4)
        fn<-paste("C:/TEMP/ComancheCimarrone/PCA/pca1_",tid,".tif",sep = "")
        predict(y,pca1,fun = pred1,filename = fn,overwrite = T)
        return(fn)
      }
      return(NULL)
    })

# Review results ------------ 
    v1<-vrt(unique(do.call(c,fl)))
    par(mfrow =c(2,2))
    i<-1;for(i in seq(0,dim(v1)[3]-3,by = 3)){plotRGB(subset(v1,1:3+i),stretch = "lin",main = paste(1:3+i,collapse = "-"))}
  