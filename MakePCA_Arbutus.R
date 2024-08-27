## Setup -------------------------------------
  setwd("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full2")
  library(terra)
  library(AFFT)
  library(snowfall)
 fp<-GetFootprints(fl<-list.files("r/r_tiles",full.names = T),return.as.list = T,quiet = T)
  fp<-vect(fp)
  bnd<-vect("Y:/MPSG_VegMapping/Data/Spatial/Mapping boundaries final/MPSG_mapping_boundaries_final.shp")
  bnd0<- bnd<-project(bnd,fp)
  bnd<-bnd[bnd$FORESTNUMB %in% c("15","06","10"),]#
  bnd<-buffer(bnd,100000)
  
  fp$TileID<- substr(fl,18,20)
  curtiles<-intersect(fp,bnd)

  nlcdmask<-rast("Y:/MPSG_VegMapping/Data/Raster/Source/nlcd/nlcd_2021_crop.tif")
  nlcdmask<-crop(nlcdmask,curtiles)
  nlcdmask<-mask(nlcdmask,curtiles)
  nlcdmask<-subset(nlcdmask,1)
  nlcdmask2<-nlcdmask %in% c("Unclassified","Barren Land","Deciduous Forest","Evergreen Forest","Mixed Forest","Shrub/Scrub","Herbaceous","Hay/Pasture","Woody Wetlands","Emergent Herbaceous Wetlands")
  nlcdmask2<-mask(nlcdmask2,nlcdmask2,maskvalue = F)
  
  samp1<-spatSample(nlcdmask2,100000,"random",xy = T,na.rm = T)
  
  fl<-list.files("C:/TEMP/RockyMountains/Corrected",pattern = ".tif",full.names = F)
  fl<-fl[!grepl("artifact",fl)]
  
  bb<-unique(apply(reshape2::colsplit(fl,"_",as.character(1:4))[,1:2],1,function(x){paste(x,collapse = "_")}))

  fl<-list.files("C:/TEMP/RockyMountains/Corrected",pattern = ".tif",full.names = T)
  vl<-pbapply::pblapply(bb,function(b){
    fl2<-fl[grepl(b,fl)]
    v1<-vrt(fl2,filename = paste(b,".vrt",sep = ""),overwrite = T)
    v1
  })  
  v2<-rast(vl)

  
  xy<-samp1[,1:2]
  suppl<-rast("D:/LocalNaip/Suppl/suppl.vrt")
  suppl<-crop(suppl,nlcdmask)
  nlcdmask<-crop(nlcdmask,suppl)
  suppl<-mask(suppl,nlcdmask)
  
  
  
## identify problematic bands ----------------------

  ### Noise ------------
  # if(!exists("C:/TEMP/RockyMountains/PCA/noisescores.csv")){
  #   vll<-sapply(vl,sources);names(vll)<-sapply(vl,names)
  #   sfInit(T,3)
  #   sfLibrary(terra);sfLibrary(AFFT)
  #   sfExport("vll")
  #   write("Scoring Noise","D:/LocalNaip/Progress.txt",append = F)
  #   write("","C:/TEMP/RockyMountains/Corrected/noises.csv",append = F)
  #   noises<-sfClusterApplyLB(vll,function(x){
  #     y<-rast(x)
  #     z<-ScoreNoise(y)
  #     n<-c(File = x,z)
  #     write(paste(n,collapse = ","),"C:/TEMP/RockyMountains/Corrected/noises.csv",append = T)
  #     write(paste(x,":",which(vll == x),"out of",length(vll)),"D:/LocalNaip/Progress.txt",append = T);n
  #   })
  #     nm<-names(noises[[1]])
  #   
  #   noises<-read.csv("C:/TEMP/RockyMountains/Corrected/noises.csv",header = F);colnames(noises)<-nm
  #   write.csv(noises,"C:/TEMP/RockyMountains/Corrected/noisescores.csv",row.names = F)
  #   sfStop()
  # 
  # }else{noises<-read.csv("C:/TEMP/RockyMountains/Corrected/noisescores.csv")}
  # ### Artifacts ---------------
  if(!file.exists("C:/TEMP/RockyMountains/Corrected/artifactscores.csv")){
    v1<-vrt(fl[grepl(bb[1],fl)]);windows(10,10);plot(v1);plot(bnd0,add = T);e1<-draw();dev.off()
    artifacts<-data.frame(do.call(rbind,pbapply::pblapply(bb,function(b){
      fl2<-fl[grepl(b,fl)]
      v1<-vrt(fl2)
      par(mfrow =c(2,2))
      plot(v1,main = b); plot(bnd0,add = T)
      plotRGB(c(v1,v1,v1),stretch = "lin")
      plot(v1,ext = e1);plot(bnd0,add = T)
      plotRGB(c(v1,v1,v1),stretch = "lin",ext = e1)
      plot(bnd0,border = "red",add = T)
      cat("\n################\n ",b,"\n###############\n\n")
      return(c(Phenology = readline("Phenology Score: "), Flightline = readline("Flightline score: "),image = b))
    })))
    write.csv(artifacts,"C:/TEMP/RockyMountains/Corrected/artifactscores.csv",row.names = F)
  }else{artifacts<-read.csv("C:/TEMP/RockyMountains/Corrected/artifactscores.csv")}
  
  
  # rn<-gsub("C:/TEMP/RockyMountains/Corrected/","",noises$File)
  # rn<-gsub(".vrt","",rn)
  # rownames(noises)<-rn
  # noises<-data.frame(noises)
  # nl<-lapply(rownames(noises),function(n){par(mfrow =c(2,2),mar =c(0,0,0,0));y<-subset(v2,n)
  #   plot(y,main = n,legend = F);plot(vect(e1),add = T)
  #   plot(y,ext = e1,legend = F);plot(vect(e2),add = T)
  #   plot(y,ext = e2,legend = F);plot(vect(e3),add = T);
  #   plot(y,ext = e3,legend = F);readline(n)})

  
  patchy<-artifacts$image[artifacts$Phenology %in% c(4,5)]
  stripey<-artifacts$image[artifacts$Flightline %in% c(4,5)]
  #noisy<-rownames(noises)[noises$X10. > 650]
  drop.me<-unique(c(patchy,stripey))#,noisy
  drop.me<-gsub(".vrt","",drop.me)
  drop.me<-gsub("C:/TEMP/RockyMountains/Corrected/","",drop.me)
  drop.me<-unique(drop.me)
  
  # nl<-lapply(drop.me,function(n){par(mfrow =c(2,2),mar =c(0,0,0,0));y<-subset(v2,n)
  #    plot(y,main = n,legend = F,ext = e1);plot(vect(e2),add = T)
  #    plot(y,ext = e2,legend = F)
  #    plot(y,ext = e3,legend = F);plot(vect(e4),add = T)
  #    plot(y,ext = e4,legend = F);readline(n)})
   
  drop.me0<-c('bri_f-0.25','bri_mean','bri_med','bri_Q95','n_f-0.25','r_f-1.5')
## Scroll through all vrt layers, build rast object vrt containing all layers --------------
  setwd("C:/TEMP/RockyMountains")
  
  cv<-rownames(noises)[!rownames(noises) %in% drop.me];names(cv)<-cv
  fl<-list.files("Corrected",full.names = T,pattern = "vrt")
  cv<-gsub("Corrected/","",fl);cv<-reshape2::colsplit(cv,"_",c("one","two","three"))[,c(1:2)]
  cv<-unique(paste(cv$one,cv$two,sep = "_"))
  v2<-rast(fl)

 
## Build princomp -------------
  ### Gather Info ---------------- 
  xy<-df[,1:2]
    system.time(mat2<-extract(suppl,xy))
    
    sv2<-sources(v2);names(sv2)<-names(v2)
    
    sfInit(T,5)
    sfLibrary(terra)
    sfExport("xy")
    df1<-sfClusterApplyLB(sv2,function(x){extract(rast(x),xy)})
    
    df2<-lapply(df1,function(x){x[,2]})
    mat1<-data.frame(do.call(cbind,df2))
    colnames(mat1)<-names(sv2)
    
    
    mat3<-cbind(mat1,mat2[,2:ncol(mat2)])
    mat1<-mat3[apply(mat3,1,function(x){!any(is.na(x))}),]
    rm(list =c("mat2","mat3"))
  
  ### Build PCA -------------------
    mat3<-mat1
    drop.me<-list()
    
    drop.me[[1]]<-c("jdate","sl","imgyr","stateid")
    drop.me[[2]]<-c("ndvi_Q025","ndng_Q025","g_Q95","n_med","r_Q95","n_Q95","g_mean","g_f-1.5","g_f-0.25","bri_f-1.5","ndvi_f-0.25","g_f-3","r_f-3","ndng_f-1.5","ndgr_f-3","ndgr_f-12","ndgr_f-6","ndng_f-3","ndvi_f-3")
    drop.me[[3]]<-c('ndvi_f-6','n_f-3','ndvi_f-1.5','bri_f-6','ndgr_f-0.25','ndgr_f-1.5','g_f-6','n_f-1.5','bri_f-3','ndng_f-0.25','r_f-6','r_f-0.25')
    drop.me[[4]]<-c('ndng_f-6','ndvi_f-12','ndgr_fp-0.25','bri_kurt','r_kurt','g_kurt','g_kurt','bri_kurt','r_kurt')
    drop.me[[5]]<-c("r_mean","g_med")
    drop.me[[6]]<-c("n_Q025")
    drop.me[[7]]<-c("bri_f-0.25","bri_Q95","bri_med","g_Q025","bri_mean")
   # drop.me[[8]]<-c("n_f-0.25","r_f-1.5")
    #drop.me[[9]]<-c("r_med","ndvi_Q95","r_Q025","bri_Q025")#"n_mean"
    
    mat3<-mat1[,!colnames(mat1) %in% do.call(c,drop.me)]
    pca1<-princomp(mat3)
     
    save(pca1,file = "C:/TEMP/RockyMountains/PCA/pca1.RData")
    load("C:/TEMP/RockyMountains/PCA/pca1.RData")
## Go through tile IDs, stack, predict PCA to tile.--------------
  fl<-list.files("C:/TEMP/RockyMountains/Corrected",recursive = T,full.names = T,pattern = '.tif')
  fl<-fl[SDMap::grepll(names(pca1$center),fl,any)]
  
  pred1<-function(x,y,keepvec = which(cumsum(x$sdev/sum(x$sdev))<.991)){z<-round(predict(x,y)[,keepvec] * 100,0)}
  ct0<-unique(curtiles$TileID)
  ct0<-ct0[!ct0 %in% c("181","213")]
  sfInit(T,19)
  sfLibrary(terra)
  
  fs<-file.size(list.files("C:/TEMP/RockyMountains/PCA",pattern = "a.tif",full.names = T))
  names(fs)<-substr(list.files("C:/TEMP/RockyMountains/PCA",pattern = "a.tif",full.names = F),6,8)
  fs<-fs[order(fs,decreasing = T)]
  ct0<-names(fs)
  
  sfExport(list = c("fl","pred1","pca1","ct0"))
  
 sfClusterApplyLB(ct0,function(x){
    cat(x,"\n")
    terraOptions(memfrac = 1/21)
    fl2<-fl[grepl(x,fl)]
    y<-rast(fl2)
    on1<-paste("C:/TEMP/RockyMountains/PCA/tile_",x,"a.tif",sep = "")
    y2<-predict(y,pca1,fun = pred1,filename = on1,overwrite = T)
    write(paste(x,":",which(ct0 == x),"out of",length(ct0)),"D:/LocalNaip/Progress.txt",append = T)
  })
  
  sfStop()
 