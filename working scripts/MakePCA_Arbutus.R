## Setup -------------------------------------
  setwd("C:/TEMP/ComancheCimarrone")
  library(terra)
  library(AFFT)
  library(snowfall)

fp<-vect("D:/LocalNaip/3_projected/footprints.shp")
  bnd<-vect("Y:/MPSG_VegMapping/Data/Spatial/Mapping boundaries final/MPSG_mapping_boundaries_final.shp")
  bnd0<- bnd<-project(bnd,fp)
  bnd<-bnd[bnd$FORESTNAME %in% c("Comanche and Cimarron National Grasslands"),]#
  bnd<-buffer(bnd,100000)
  
  curtiles<-intersect(fp,bnd)

  nlcdmask<-rast("Y:/MPSG_VegMapping/Data/Raster/Source/nlcd/nlcd_2021_crop.tif")
  nlcdmask<-crop(nlcdmask,curtiles)
  nlcdmask<-mask(nlcdmask,curtiles)
  nlcdmask<-subset(nlcdmask,1)
  nlcdmask2<-nlcdmask %in% c("Unclassified","Barren Land","Deciduous Forest","Evergreen Forest","Mixed Forest","Shrub/Scrub","Herbaceous","Hay/Pasture","Woody Wetlands","Emergent Herbaceous Wetlands")
  nlcdmask2<-mask(nlcdmask2,nlcdmask2,maskvalue = F)
  
  samp1<-spatSample(nlcdmask2,100000,"random",xy = T,na.rm = T)
  
  fl<-list.files("C:/TEMP/ComancheCimarrone/Corrected",pattern = ".tif",full.names = F)
  fl<-fl[!grepl("artifact",fl)]
  
  bb<-unique(apply(reshape2::colsplit(fl,"_",as.character(1:4))[,1:2],1,function(x){paste(x,collapse = "_")}))
  
  fl<-list.files("C:/TEMP/ComancheCimarrone/Corrected",pattern = ".tif",full.names = T)
  fl<-fl[!grepl("artifact",fl)]
  vl<-pbapply::pblapply(bb,function(b){
    fl2<-fl[grepl(b,fl)]
    s1<-sprc(fl2)
    m1<-merge(s1,filename = paste(b,"_merge.tif"),overwrite = T)
    #v1<-vrt(fl2,filename = paste(b,".vrt",sep = ""),overwrite = T)
    m1
  })  

  fl<-list.files("C:/TEMP/ComancheCimarrone/Merged",pattern = ".tif",full.names = T)
  v2<-rast(fl)
  
  fl<-list.files("C:/TEMP/ComancheCimarrone/Merged",pattern = ".tif",full.names = F)
  fl<-gsub(" _merge.tif","",fl)
  names(v2)<-fl
  xy<-samp1[,1:2]
  suppl<-rast("D:/LocalNaip/Suppl/suppl.vrt")
  suppl<-crop(suppl,nlcdmask)
  nlcdmask<-crop(nlcdmask,suppl)
  suppl<-mask(suppl,nlcdmask)
  
  
  
## identify problematic bands ----------------------


 ### Artifacts ---------------
  if(!file.exists("C:/TEMP/ComancheCimarrone/Corrected/artifactscores.csv")){
    windows(10,10);plot(v2,4);plot(bnd0,add = T);e1<-draw();dev.off()
    artifacts<-data.frame(do.call(rbind,pbapply::pblapply(bb,function(b){

      fl2<-fl[grepl(paste("/",b," ",sep = ""),fl)]
      v1<-rast(fl2)
      par(mfrow =c(2,2))
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
  drop.me<-unique(c(patchy,stripey))
  drop.me<-unique(drop.me)
  
   # nl<-lapply(drop.me,function(n){par(mfrow =c(2,1),mar =c(0,0,0,0));y<-subset(v2,n)
   #    plot(y,main = n,legend = F);plot(vect(e1),add = T)
   #    plot(y,ext = e1,legend = F)
   #    #plot(y,ext = e2,legend = F);plot(vect(e4),add = T)
   #    #plot(y,ext = e4,legend = F)
   #    ;readline(n)})
   
## Scroll through all vrt layers, build rast object vrt containing all layers --------------
  setwd("C:/TEMP/ComancheCimarrone")
  
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
  
## Build princomp -------------
  ### Gather Info ---------------- 

    system.time(mat2<-extract(suppl,xy))
    
    sv2<-vl
    
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
    
    drop.me[[1]]<-c("jdate","sl","imgyr","stateid",unique(c(patchy,stripey)))

    mat3<-mat1[,!colnames(mat1) %in% do.call(c,drop.me)]
    pca1<-princomp(mat3)
     
    save(pca1,file = "C:/TEMP/ComancheCimarrone/PCA/pca1.RData")
    load("C:/TEMP/ComancheCimarrone/PCA/pca1.RData")

  pred1<-function(x,y,keepvec = which(cumsum(x$sdev/sum(x$sdev))<.991)){z<-round(predict(x,y)[,keepvec] * 100,0)}
  pca<-predict(subset(v2,names(pca1$center)),pca1,fun = pred1,filename = "C:/TEMP/ComancheCimarrone/PCA/pca1.tif",overwrite = T,cores = 30)
  