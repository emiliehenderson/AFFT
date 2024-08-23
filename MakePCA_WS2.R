## Setup -------------------------------------
  setwd("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full2")
  library(terra)
  library(AFFT)
  library(snowfall)
  bnd0<-bnd<-vect("Y:/MPSG_VegMapping/Data/Spatial/Mapping boundaries final/MPSG_mapping_boundaries_final.shp")
  fp<-GetFootprints(fl<-list.files("r/r_tiles",full.names = T),return.as.list = T,quiet = T)
  fp<-vect(fp)
  bnd<-project(bnd,fp)
  bnd0<-project(bnd0,fp)
  bnd<-bnd[bnd$FORESTNUMB %in% c(NA),]# c("07","18")
  bnd<-buffer(bnd,100000)
  plot(bnd0)
  fp$TileID<- substr(fl,18,20)
  curtiles<-intersect(fp,bnd)
  
  curtiles<-curtiles[!curtiles$TileID %in% "196",]
  nlcdmask<-rast("Y:/MPSG_VegMapping/Data/Raster/Source/nlcd/nlcd_2021_crop.tif")
  nlcdmask<-crop(nlcdmask,curtiles)
  nlcdmask<-mask(nlcdmask,curtiles)
  nlcdmask<-subset(nlcdmask,1)
  nlcdmask2<-nlcdmask %in% c("Unclassified","Barren Land","Deciduous Forest","Evergreen Forest","Mixed Forest","Shrub/Scrub","Herbaceous","Hay/Pasture","Woody Wetlands","Emergent Herbaceous Wetlands")
  nlcdmask2<-mask(nlcdmask2,nlcdmask2,maskvalue = F)
#   
#   samp1<-spatSample(nlcdmask2,100000,"random",xy = T,na.rm = T)
#   
#   suppl<-rast("C:/TEMP/tmp.tif")
#   suppl<-crop(suppl,nlcdmask)
#   nlcdmask<-crop(nlcdmask,suppl)
#   suppl<-mask(suppl,nlcdmask)
#   
#   fl<-list.files("C:/TEMP/ComancheCimarrone/Corrected",pattern = ".tif",full.names = F)
#   
#   bb<-unique(apply(reshape2::colsplit(fl,"_",as.character(1:4))[,1:2],1,function(x){paste(x,collapse = "_")}))
#   bb<-bb[!grepl("artifact",bb)]
#   fl<-list.files("C:/TEMP/ComancheCimarrone/Corrected",pattern = ".tif",full.names = T)
#   vl<-pbapply::pblapply(bb,function(b){
#     fl2<-fl[grepl(b,fl)]
#     v1<-vrt(fl2,filename = paste("C:/TEMP/ComancheCimarrone/Corrected/",b,".vrt",sep = ""),overwrite = T)
#     v1
#   })  
#   
#   
#   
# ## identify problematic bands ----------------------
# 
#   ### Noise ------------
#   # if(!file.exists("C:/TEMP/ComancheCimarrone/Corrected/noises.csv")){
#   #   vll<-sapply(vl,sources);names(vll)<-sapply(vl,names)
#   #   sfInit(T,3)
#   #   sfLibrary(terra);sfLibrary(AFFT)
#   #   sfExport("vll")
#   #   write("Scoring Noise","D:/LocalNaip/Progress.txt",append = F)
#   #   write("","C:/TEMP/ComancheCimarrone/Corrected/noises.csv",append = F)
#   #   # noises<-sfClusterApplyLB(vll,function(x){
#   #   #   y<-rast(x)
#   #   #   z<-ScoreNoise(y)
#   #   #   n<-c(File = x,z)
#   #   #   write(paste(n,collapse = ","),"C:/TEMP/ComancheCimarrone/Corrected/noises.csv",append = T)
#   #   #   write(paste(x,":",which(vll == x),"out of",length(vll)),"D:/LocalNaip/Progress.txt",append = T);n
#   #   # })
#   #   #   nm<-names(noises[[1]])
#   #   # 
#   #   noises<-read.csv("C:/TEMP/ComancheCimarrone/Corrected/noises.csv",header = F);colnames(noises)<-nm
#   #   write.csv(noises,"C:/TEMP/ComancheCimarrone/Corrected/noisescores.csv",row.names = F)
#   #   sfStop()
#   # 
#   # }else{noises<-read.csv("C:/TEMP/ComancheCimarrone/Corrected/noisescores.csv",header = F)}
#   ### Artifacts ---------------
#   # if(!file.exists("C:/TEMP/ComancheCimarrone/Corrected/artifactscores.csv")){
#   #   v1<-vrt(fl[grepl(bb[1],fl)]);windows(10,10);plot(v1);plot(bnd0,add = T);e1<-draw();dev.off()
#   #   artifacts<-data.frame(do.call(rbind,pbapply::pblapply(bb,function(b){
#   #     fl2<-fl[grepl(b,fl)]
#   #     v1<-vrt(fl2)
#   #     par(mfrow =c(2,2))
#   #     plot(v1,main = b); plot(bnd0,add = T)
#   #     plotRGB(c(v1,v1,v1),stretch = "lin")
#   #     plot(v1,ext = e1);plot(bnd0,add = T)
#   #     plotRGB(c(v1,v1,v1),stretch = "lin",ext = e1)
#   #     plot(bnd0,border = "red",add = T)
#   #     cat("\n################\n ",b,"\n###############\n\n")
#   #     return(c(Phenology = readline("Phenology Score: "), Flightline = readline("Flightline score: "),image = b))
#   #   })))
#   #   write.csv(artifacts,"C:/TEMP/ComancheCimarrone/Corrected/artifactscores.csv",row.names = F)
#   # }else{artifacts<-read.csv("C:/TEMP/ComancheCimarrone/Corrected/artifactscores.csv")}
#   # 
#   
#   rn<-gsub("C:/TEMP/ComancheCimarrone/Corrected/","",noises$V1)
#   rn<-gsub(".vrt","",rn)
#   rownames(noises)<-rn
#   noises<-data.frame(noises)
#   v2<-rast(pbapply::pblapply(rn,function(x){
#     vn<-paste("C:/TEMP/ComancheCimarrone/Corrected/",x,".vrt",sep = "")
#     fl2<-fl[grepl(x,fl)]
#     cl<-vrt(fl2,filename = vn,overwrite = T)
#     cl
#   }))
#   
#   patchy<-artifacts$image[artifacts$Phenology %in% c(4,5)]
#   stripey<-artifacts$image[artifacts$Flightline %in% c(4,5)]
#   noisy<-rownames(noises)[noises$X10. > 650]
#   drop.me<-unique(c(patchy,stripey))#,noisy))
#   drop.me<-gsub(".vrt","",drop.me)
#   drop.me<-gsub("C:/TEMP/ComancheCimarrone/Corrected/","",drop.me)
#   drop.me<-unique(drop.me)
# 
#   windows(10,10);plot(subset(v2,drop.me[1]));plot(bnd0,add = T, border = "red");e1<-draw();e2<-draw();e3<-draw();e4<-draw();dev.off()
#   
#   nl<-lapply(drop.me,function(n){
#     layout(matrix(c(1,1,2,3,4,5),byrow = T,ncol = 3))
#     y<-subset(v2,n)
#     plot(y,main = n,legend = F);plot(bnd0,border = "orange",add = T);plot(vect(e1),add = T)
#     plot(y,ext = e1,legend = F);plot(bnd0,border = "orange",add = T);plot(vect(e2),add = T)
#     plot(y,ext = e2,legend = F);plot(bnd0,border = "orange",add = T);plot(vect(e3),add = T);
#     plot(y,ext = e3,legend = F);plot(bnd0,border = "orange",add = T);plot(vect(e4),add = T);
#     plot(y,ext = e4,legend = F);plot(bnd,border = "orange",add = T);readline(n)})
#     
# 
#     
#   
#  # drop.me0<-c('g_f-0.25','bri_fp-20','n_fp-20')
#   ## The AFFT layers identified above appear fairly useable. Keep in for the next steps.
#   
## Scroll through all vrt layers, build rast object vrt containing all layers --------------
  setwd("C:/TEMP/ComancheCimarrone")

  fl<-list.files("Corrected",full.names = T,pattern = ".tif");fl<-fl[!grepl("art",fl)]
  cv<-gsub("Corrected/","",fl);cv<-reshape2::colsplit(cv,"_",c("one","two","three"))[,c(1:2)]
  cv<-unique(paste(cv$one,cv$two,sep = "_"))
  
  v2<-rast(pbapply::pblapply(cv,function(x){
    vn<-paste("Corrected/",x,".vrt",sep = "")
    fl2<-fl[grepl(x,fl)]
    cl<-vrt(fl2,filename = vn,overwrite = T)
    cl
  }))
#  
# ## Build princomp -------------
#   ### Gather Info ---------------- 
#   xy<-samp1[,1:2]
#     system.time(mat2<-extract(suppl,xy))
#     
#     sv2<-sources(v2);names(sv2)<-names(v2)
#     
#     sfInit(T,5)
#     sfLibrary(terra)
#     sfExport("xy")
#     df1<-sfClusterApplyLB(sv2,function(x){extract(rast(x),xy)})
#     
#     df2<-lapply(df1,function(x){x[,2]})
#     mat1<-data.frame(do.call(cbind,df2))
#     colnames(mat1)<-names(sv2)
#     
#     
#     mat3<-cbind(mat1,mat2[,2:ncol(mat2)])
#     mat1<-mat3[apply(mat3,1,function(x){!any(is.na(x))}),]
#     rm(list =c("mat2","mat3"))
#   
#   ### Build PCA -------------------
#     mat3<-mat1
#     mat1<-mat1[mat1$`r_f-20`>4000,]
#     drop.me<-list()
#     
#     drop.me[[1]]<-c("jdate","sl","imgyr","stateid")
#     
#     mat3<-mat1[,!colnames(mat1) %in% do.call(c,drop.me)]
#     pca1<-princomp(mat3)
#     
   # save(pca1,file = "C:/TEMP/ComancheCimarrone/PCA/pca1.RData")
    load("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full2/ComancheCimarrone/PCA/pca1.RData")
    ## Go through tile IDs, stack, predict PCA to tile.--------------
    fl<-list.files("C:/TEMP/ComancheCimarrone/Corrected",recursive = T,full.names = T,pattern = '.tif')
    fl<-fl[SDMap::grepll(names(pca1$center),fl,any)]
    
    pred1<-function(x,y,keepvec = which(cumsum(x$sdev/sum(x$sdev))<.991)){z<-round(predict(x,y)[,keepvec] * 100,0)}
    ct0<-unique(curtiles$TileID)
    ct0<-ct0[!ct0 %in% c("181","213")]
    sfInit(T,7)
    sfLibrary(terra)
    
    fs<-file.size(list.files("C:/TEMP/ComancheCimarrone/PCA",pattern = "f.tif",full.names = T))
    names(fs)<-substr(list.files("C:/TEMP/ComancheCimarrone/PCA",pattern = "f.tif",full.names = F),6,8)
    fs<-fs[order(fs,decreasing = T)]
    ct0<-names(fs)
    
    sfExport(list = c("fl","pred1","pca1","ct0"))
    
    sfClusterApplyLB(ct0,function(x){
      cat(x,"\n")
      terraOptions(memfrac = 1/8)
      fl2<-fl[grepl(x,fl)]
      y<-rast(fl2)
      on1<-paste("C:/TEMP/ComancheCimarrone/PCA/tile_",x,"f.tif",sep = "")
      y2<-predict(y,pca1,fun = pred1,filename = on1,overwrite = T)
      write(paste(x,":",which(ct0 == x),"out of",length(ct0)),"D:/LocalNaip/Progress.txt",append = T)
    })
    
    sfStop()
    