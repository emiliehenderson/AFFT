# Setup -------------------------------------
    
    library(terra)
    library(AFFT)
    library(snowfall)
    SDMap::MakeDirIfNeeded("RockyMountains","C:/TEMP")
    
  ###functions----
    pred1<-function(x,y,keepvec = which(cumsum(x$sdev/sum(x$sdev))<.99)){z<-round(predict(x,y)[,keepvec] * 100,0)}

  ###boundaries -----
    bnd<-vect("Y:/MPSG_VegMapping/Data/Spatial/Mapping boundaries final/MPSG_mapping_boundaries_final.shp")
    bnd0<- bnd<-project(bnd,rast(list.files("Y:/MPSG_VegMapping/Data/Raster/Predictors/elevation",full.names = T,pattern = "tif")[1]))
    bnd<-bnd[bnd$FORESTNUMB %in% c("15","10","06"),]#c("18")
    bnd<-buffer(bnd,100000)
    bnd<-aggregate(bnd)
    fl<-list.files("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/r/r_tiles",full.names = T)
    names(fl)<-substr(list.files("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/r/r_tiles",full.names = F),8,10)
  
    fp<-GetFootprints(fl,return.as.list = F)
    fp$TileID<-names(fl)
    crs(fp)<-crs(bnd)
    bnd1<-crop(fp,bnd)
  ###nlcd ----- 
    nlcdmask<-rast("Y:/MPSG_VegMapping/Data/Raster/Source/nlcd/nlcd_2021_crop.tif")
    nlcdmask<-crop(nlcdmask,bnd)
    nlcdmask<-mask(nlcdmask,bnd)
    nlcdmask<-subset(nlcdmask,1)
    nlcdmask2<-nlcdmask %in% c("Unclassified","Barren Land","Deciduous Forest","Evergreen Forest","Mixed Forest","Shrub/Scrub","Herbaceous","Hay/Pasture","Woody Wetlands","Emergent Herbaceous Wetlands")
    nlcdmask2<-mask(nlcdmask2,nlcdmask2,maskvalue = F)
  ###point samples ----- 
    samp1<-spatSample(nlcdmask2,150000,"random",xy = T,na.rm = T)
    
    if(!any(grepl("vrt",list.files("C:/TEMP/RockyMountains/Corrected")))){
      fl<-list.files("C:/TEMP/RockyMountains/Corrected",pattern = ".tif",full.names = T)
      nm<-gsub("C:/TEMP/RockyMountains/Corrected/","",fl)
      nm<-reshape2::colsplit(nm,"_",c(1:4))
      nm<-paste(nm[,1],nm[,2],sep = "_")
      fll<-split(fl,nm)
      nm<-names(fll);names(nm)<-nm
      v1<-pbapply::pblapply(nm,function(x){cat(x);vrt(fll[[x]],filename = paste("C:/TEMP/RockyMountains/Corrected/",x,".vrt",sep = ""),overwrite = T)})
    }
    fl<-list.files("C:/TEMP/RockyMountains/Corrected",pattern = ".vrt",full.names = T)
    
      v2<-rast(fl)
      xy<-samp1[,1:2]
     suppl<-rast("D:/LocalNaip/Suppl/suppl.vrt")
     suppl<-crop(suppl,nlcdmask)
     nlcdmask<-crop(nlcdmask,suppl)
     #suppl<-mask(suppl,nlcdmask)
     suppl<-c(suppl,nlcdmask)
# Artifacts ---------------
  if(!file.exists("C:/TEMP/Dakotas/Corrected/artifactscores.csv")){
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
    write.csv(artifacts,"C:/TEMP/RockyMountains/Corrected/artifactscores.csv",row.names = F)
  }else{artifacts<-read.csv("C:/TEMP/RockyMountains/Corrected/artifactscores.csv")}
  

  patchy<-artifacts$image[artifacts$Phenology %in% c(4,5)]
  stripey<-artifacts$image[artifacts$Flightline %in% c(4,5)]


# Build princomp -------------
  ### Gather Info ---------------- 
    
    sfInit(T,5)
    sfLibrary(terra)
    fl<-list.files("C:/TEMP/RockyMountains/Corrected",pattern = ".tif",full.names = T)
    nm<-list.files("C:/TEMP/RockyMountains/Corrected",pattern = ".tif",full.names = F)
    nm<-reshape2::colsplit(nm,"_",c(1:3))[,3]
    nm<-gsub(".tif","",nm)
    fll<-split(fl,nm)
    pl<-extract(bnd1,xy)
    pl$TileID[is.na(pl$TileID)]<--1
    pll<-split(pl,pl$TileID)
    xyl<-lapply(pll,function(x){xy[x[,1],]})
    ss<-sources(suppl)
    sfExport("fll");sfExport("xyl");sfExport("ss")
    write("",file = "C:/TEMP/Progress.txt",append = F)
    df1<-sfClusterApplyLB(tids,function(x){
      e0<-extract(rast(fll[[x]]),xyl[[x]])
      e0<-e0[!apply(e0,1,function(x){any(is.na(x))}),]
      e1<-extract(rast(ss[1]),xyl[[x]])
      e1<-e1[!apply(e1,1,function(x){any(is.na(x))}),]
      e2<-merge(e1,e0,by = "ID" )
      write(x,file = "C:/TEMP/Progress.txt",append = T)
      e2
    })
    sfStop()
    mat1<-data.frame(do.call(rbind,df1))
    mat1<-unique(mat1)
    colnames(mat1)<-colnames(df1[[1]])
    mat1<-mat1[,2:ncol(mat1)]
  ### Build model -------------------
    drop.me<-list()
    
    drop.me[[1]]<-c("jdate","sl","imgyr","stateid",unique(c(patchy,stripey)))
    drop.me[[2]]<-c("ndgr_f-0.25","ndng_f-0.25","ndvi_f-0.25")
    mat3<-mat1[,!colnames(mat1) %in% do.call(c,drop.me)]
    pca1<-princomp(mat3)
     
    save(pca1,file = "C:/TEMP/RockyMountains/PCA/pca1.RData")

# Make iamges -------------------    
setwd("C:/TEMP/RockyMountains")
    fl<-list.files("Corrected",recursive = T, full.names = T,pattern = ".tif")
    fl<-fl[!grepl(".csv",fl)]
    nm<-list.files("Corrected",recursive = T, full.names = F,pattern = ".tif")
    nm<-nm[!grepl(".csv",nm)]
    nm<-reshape2::colsplit(nm,"_",c(1:3));tids<-gsub(".tif","",nm[,3])
    
    vn<-paste(nm[,1],nm[,2],sep = "_")
    df<-data.frame(fl,tids,vn)
    df<-df[!grepl("artifact",df[,1]),]
    fllt<-tapply(df,df$tid,function(x){y<-x[,1];names(y)<-x[,3];y})
    sfInit(T,length(unique(tids)))
    sfLibrary(AFFT)
    sfLibrary(terra)
    sfExport("pred1","pca1","fllt")
  
  write("",file = "C:/TEMP/Progress.txt",append = F)
 # tids<-bnd1$TileID[order(expanse(bnd1),decreasing = T)]
  #tids<-crashtiles
  {
    fl<-sfClusterApplyLB(unique(tids),function(x){
      y<-rast(fllt[[x]])
      y<-subset(y,names(pca1$center))
      if(!any(is.na(y@ptr@.xData$range_max))){
        fn<-paste("C:/TEMP/RockyMountains/PCA/pca1",x,".tif",sep = "")
        predict(y,pca1,fun = pred1,filename = fn,overwrite = T)
        write(x,file = "C:/TEMP/Progress.txt",append = T)
        return(fn)
      }
      write(x,file = "C:/TEMP/Progress.txt",append = T)
      return(NULL)
    })
    sfStop()
}
# Review results ------------ 
    v1<-vrt(unique(do.call(c,fl)))
    par(mfrow =c(1,3))
    i<-1;for(i in seq(0,dim(v1)[3]-3,by = 3)){plotRGB(subset(v1,1:3+i),stretch = "lin",main = paste(1:3+i,collapse = "-"))}
  