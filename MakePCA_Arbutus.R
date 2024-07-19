## Setup -------------------------------------
  setwd("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full")
  library(terra)
  library(AFFT)
  library(snowfall)
  bnd0<-bnd<-vect("Y:/MPSG_VegMapping/Data/Spatial/Mapping boundaries final/MPSG_mapping_boundaries_final.shp")
  fp<-GetFootprints(fl<-list.files("r/r_tiles",full.names = T),return.as.list = T,quiet = T)
  fp<-vect(fp)
  bnd<-project(bnd,fp)
  bnd0<-project(bnd0,fp)
  bnd<-bnd[bnd$FORESTNUMB %in% c("07","18"),]#c("15","06","10")
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
  

  t1<-subset(v2,"bri_f-20")
  tmp1<-extract(t1,samp1[,c(1:2)])
  df<-data.frame(samp1,tmp1)
  df$tf<-df$bri_f.20< -6000
  tmp1<-vect(df,geom =c("x","y"))
  
  samp1<-vect(df[!df$tf,c("x","y")],geom =c("x","y"))
  
  suppl<-rast("D:/LocalNaip/Suppl/suppl.vrt")
  suppl<-crop(suppl,nlcdmask)
  nlcdmask<-crop(nlcdmask,suppl)
  suppl<-mask(suppl,nlcdmask)
  
  fl<-list.files("D:/LocalNaip/Corrected",pattern = ".tif",full.names = F)
  
  bb<-unique(apply(reshape2::colsplit(fl,"_",as.character(1:4))[,1:2],1,function(x){paste(x,collapse = "_")}))

  fl<-list.files("D:/LocalNaip/Corrected",pattern = ".tif",full.names = T)
  vl<-pbapply::pblapply(bb,function(b){
    fl2<-fl[grepl(b,fl)]
    v1<-vrt(fl2,filename = paste("D:/LocalNaip/Corrected/",b,".vrt",sep = ""),overwrite = T)
    v1
  })  
  
  
  
## identify problematic bands ----------------------

  ### Noise ------------
  if(!exists("D:/LocalNaip/PCA/noisescores.csv")){
    vll<-sapply(vl,sources);names(vll)<-sapply(vl,names)
    sfInit(T,3)
    sfLibrary(terra);sfLibrary(AFFT)
    sfExport("vll")
    write("Scoring Noise","D:/LocalNaip/Progress.txt",append = F)
    write("","D:/LocalNaip/PCA/noises.csv",append = F)
    noises<-sfClusterApplyLB(vll,function(x){
      y<-rast(x)
      z<-ScoreNoise(y)
      n<-c(File = x,z)
      write(paste(n,collapse = ","),"C:/TEMP/PCA3/noises.csv",append = T)
      write(paste(x,":",which(vll == x),"out of",length(vll)),"D:/LocalNaip/Progress.txt",append = T);n
    })
      nm<-names(noises[[1]])
    
    noises<-read.csv("C:/TEMP/PCA3/noises.csv",header = F);colnames(noises)<-nm
    write.csv(noises,"C:/TEMP/PCA3/noisescores.csv",row.names = F)
    sfStop()
  
  }else{noises<-read.csv("C:/TEMP/PCA3/noisescores.csv")}
  ### Artifacts ---------------
  if(!file.exists("D:/LocalNaip/PCA/artifactscores.csv")){
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
    write.csv(artifacts,"D:/LocalNaip/PCA/artifactscores.csv",row.names = F)
  }else{artifacts<-read.csv("D:/LocalNaip/PCA/artifactscores.csv")}
  
  
  rn<-gsub("C:/TEMP/Corrected/","",noises$File)
  rn<-gsub(".vrt","",rn)
  rownames(noises)<-rn
  noises<-data.frame(noises)
  nl<-lapply(rownames(noises),function(n){par(mfrow =c(2,2),mar =c(0,0,0,0));y<-subset(v2,n)
    plot(y,main = n,legend = F);plot(vect(e1),add = T)
    plot(y,ext = e1,legend = F);plot(vect(e2),add = T)
    plot(y,ext = e2,legend = F);plot(vect(e3),add = T);
    plot(y,ext = e3,legend = F);readline(n)})

  
  patchy<-artifacts$image[artifacts$Phenology %in% c(4,5)]
  stripey<-artifacts$image[artifacts$Flightline %in% c(4,5)]
  noisy<-rownames(noises)[noises$X10. > 650]
  drop.me<-unique(c(patchy,stripey))#,noisy
  drop.me<-gsub(".vrt","",drop.me)
  drop.me<-gsub("D:/LocalNaip/Corrected/","",drop.me)
  drop.me<-unique(drop.me)
  
  drop.me0<-c('bri_fp-0.25','bri_mean','g_f-0.25','g_f-1.5','g_mean','n_f-0.25','n_fp-0.25','ndgr_f-0.25','ndgr_f-1.5','ndng_f-0.25','ndng_med','ndvi_f-0.25','ndvi_med','ndvi_Q025','ndvi_Q95','r_fp-0.25','r_Q025','r_Q95','ndvi_mean')
## Scroll through all vrt layers, build rast object vrt containing all layers --------------
  setwd("C:/TEMP")
  
  cv<-rownames(noises)[!rownames(noises) %in% drop.me];names(cv)<-cv
  fl<-list.files("Corrected",full.names = T,pattern = "tif")
  cv<-gsub("Corrected/","",fl);cv<-reshape2::colsplit(cv,"_",c("one","two","three"))[,c(1:2)]
  cv<-unique(paste(cv$one,cv$two,sep = "_"))
  
  v2<-rast(pbapply::pblapply(cv,function(x){
    vn<-paste("C:/TEMP/Corrected/",x,".vrt",sep = "")
    fl2<-fl[grepl(x,fl)]
    cl<-vrt(fl2,filename = vn,overwrite = T)
    cl
  }))
 
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
    
    drop.me2<-c("jdate","sl","imgyr","stateid")
    drop.me3<-c("ndgr_f-20","ndgr_f-12","ndgr_f-6","g_f-6","n_f-3","n_f-1.5","ndgr_f-3","ndvi_f-3","ndgr_f-6",
                "ndvi_f-20","ndng_f-6","ndvi_f-6","g_f-12","ndng_f-3","g_f-3","ndvi_f-1.5","ndgr_f-1.5",
                "n_fp-12","bri_skew","n_f-20","n_f-12"
                ,"bri_f-1.5","ndng_f-1.5","n_f-6",
                "ndng_f-12","r_f-1.5","r_Q025","ndgr_fp-0.25",
                "n_Q025","r_Q95","r_mean","bri_Q025","bri_med",
                "g_fp-0.25","ndng_fp-0.25","g_Q025"
                #"ndng_Q025",
                #"ndng_fp-6","g_fp-6",
                #"r_med","bri_mean","g_Q95",
                #"ndng_Q95","g_med","ndng_med","bri_mean",
                #"ndng_Q025","g_Q025","g_mean","g_med",
                #"ndng_mean","ndgr_fp-3"
                )
    #drop.me3<-drop.me3[!drop.me3 %in% add.back]
    drop.me4<-c('bri_f-1.5',
                'g_f-3','g_mean','g_med','g_Q025','g_Q95',
                'n_f-1.5','n_f-20','n_f-3','n_f-6',
                'ndgr_f-1.5','ndgr_f-20','ndgr_fp-0.25',
                'ndng_f-1.5','ndng_f-3','ndng_mean','ndng_med','ndng_Q025','ndng_Q95',
                'ndvi_f-1.5','r_mean','r_med')
    drop.me5<-c("bri_Q025","bri_fp-20","bri_fp-12","bri_fp-3","bri_fp-6",
                "ndgr_f-12","ndgr_f-6","ndgr_f-3","ndgr_fp1.5","ndgr_sd",
                "r_fp-20","r_fp-12","r_fp-3","r_fp-6",
                "n_fp-20","n_fp-12","n_fp-6")
    mat3<-mat3[,!colnames(mat3) %in% c(drop.me0,drop.me2,drop.me4)]
    pca1<-princomp(mat3)
     
    save(pca1,file = "D:/LocalNaip/PCA/pca1.RData")
   
## Go through tile IDs, stack, predict PCA to tile.--------------
  fl<-list.files("C:/TEMP/Corrected",recursive = T,full.names = T,pattern = '.tif')
  fl<-fl[SDMap::grepll(names(pca1$center),fl,any)]
  
  pred1<-function(x,y,keepvec = which(cumsum(x$sdev/sum(x$sdev))<.991)){z<-round(predict(x,y)[,keepvec] * 100,0)}
  ct0<-unique(curtiles$TileID)
  
  sfInit(T,13)
  sfLibrary(terra)
  
  fs<-file.size(list.files("D:/LocalNaip/PCA",pattern = "f.tif",full.names = T))
  names(fs)<-substr(list.files("D:/LocalNaip/PCA",pattern = "f.tif",full.names = F),6,8)
  fs<-fs[order(fs,decreasing = T)]
  ct0<-names(fs)
  
  sfExport(list = c("fl","pred1","pca1","ct0"))
  
  sfClusterApplyLB(ct0,function(x){
    cat(x,"\n")
    terraOptions(memfrac = 1/13)
    fl2<-fl[grepl(x,fl)]
    y<-rast(fl2)
    on1<-paste("C:/TEMP/PCA3/tile_",x,"f.tif",sep = "")
    y2<-predict(y,pca1,fun = pred1,filename = on1,overwrite = T)
    write(paste(x,":",which(ct0 == x),"out of",length(ct0)),"D:/LocalNaip/Progress.txt",append = T)
  })
  
  sfStop()
  ## Make vrt of outputs ---------------
  ## Scroll through all layers, build vrt --------------
  fl<-list.files("C:/TEMP/PCA",pattern = "f.tif");fl<-fl[!grepl("vrt",fl)];names(fl)<-fl
  
  v2<-rast(pbapply::pblapply(fl,function(x){
    browser()
    vn<-paste(x,".vrt",sep = "")
    fl2<-list.files(x,recursive = T, full.names = T, pattern = "tif")
    print(fl2)
    cv<-vrt(fl2,filename = vn,overwrite = T)
    cv
  }))
