
## Setup -------------------------------------
  setwd("D:/LocalNaip/3_projected/")
  library(readxl)
  library(terra)
  library(AFFT)
  library(snowfall)
  seamlines <- data.frame(read_excel("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full/naip_collection_year_V2.xlsx"))
  rownames(seamlines)<-nm<-seamlines$State;names(nm)<-nm
  #fp<-vect("D:/LocalNaip/3_projected/footprints.shp")#
  fl<-list.files("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/r/r_tiles",full.names = T)
  fp<-GetFootprints(fl,return.as.list = F,quiet = T)
  tids<-gsub("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/r/r_tiles/r_tile_","",fl)
  tids<-gsub(".tif","",tids)
  fp$TileID<-tids
## Buffer Model Region and Select Tiles ------------------------------

  bnd0<-bnd<-vect("Y:/MPSG_VegMapping/Data/Spatial/Mapping boundaries final/MPSG_mapping_boundaries_final.shp")
  bnd0<-bnd<-project(bnd,fp)
  bnd<-bnd0[bnd0$FORESTNUMB %in% c("15","10","06"),]#c("18")
  plot(bnd0);plot(bnd,add = T, col = "skyblue2",border = "skyblue2")
  bnd<-buffer(bnd,100000)
  bnd<-aggregate(bnd)
  
  curtiles<-intersect(fp,bnd)  ## Need to drop tile 181
  
  cat("Getting Supplemental layers             \r")
  
  suppl<-rast("D:/LocalNaip/Suppl/suppl.vrt")
  
  tids<-curtiles$TileID;names(tids)<-tids
  drop.me<-sapply(tids,function(x){
    plot(curtiles[curtiles$TileID%in%x,],col = sample(colors()[grepl("blue",colors())],1),add = T)
    text(curtiles[curtiles$TileID%in%x,],x,col = "orange")
    readline(x)})
  drop.me<-names(drop.me[drop.me %in% "x"])
  tids<-tids[!tids %in% drop.me]
  curtiles<-curtiles[curtiles$TileID %in% tids]
  suppl<-crop(suppl,curtiles)
  suppl<-mask(suppl,curtiles)
  nlcdmask<-rast("Y:/MPSG_VegMapping/Data/Raster/Source/nlcd/nlcd_2021_crop.tif")
  nlcdmask<-crop(nlcdmask,suppl)
  nlcdmask<-mask(nlcdmask,suppl)
  nlcdmask<-subset(nlcdmask,1)
  nlcdmask2<-nlcdmask %in% c("Unclassified","Barren Land","Deciduous Forest","Evergreen Forest","Mixed Forest","Shrub/Scrub","Herbaceous","Hay/Pasture","Woody Wetlands","Emergent Herbaceous Wetlands")
  suppl2<-mask(suppl,nlcdmask2,maskvalue = F)
  suppl$imgyr<-as.factor(suppl$imgyr)
  suppl$stateid<-as.factor(suppl$stateid)
  suppl<-writeRaster(suppl,"C:/TEMP/tmp.tif",overwrite = T)
  
## Sample Area of Interest ----------
  cat("Sampling AOI               \r")
  samp1<-spatSample(suppl2$jdate,500,method = "stratified",xy = T,na.rm = T)
  samp2<-spatSample(suppl2$jdate,80000,method = "random",xy = T,na.rm = T)
  samp3<-spatSample(suppl2$imgyr,500,method = "stratified",xy = T,na.rm = T)
  samp4<-spatSample(suppl2$stateid,500,method = "stratified",xy = T,na.rm = T)

  
  
  xy<-do.call(rbind,list(samp1[,1:2],samp2[,1:2],samp3[,1:2],samp4[,1:2]))
  print(system.time(mat1<-extract(suppl,xy)))
  afftnames<-names(rast(list.files("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/r/r_tiles",full.names = T)[1]))
## Band by band -----------------
  bands<-c("r","g","n","ndgr","ndng","bri","ndvi" )
  for(b in bands){
    cat("\n\n#################\n#####  ",b,"  #####\n#################\n\n")
    ### Build vrt for only selected tiles. ---------------
    tl<-unique(curtiles$TileID)
 
    fll<-paste("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/",b,"/",b,"_tiles/",b,"_tile_",tl,".tif",sep = "")
    file.remove(list.files("C:/TEMP/temptiles",full.names = T))
    fllt<-paste("C:/TEMP/temptiles/",b,"__",tl,sep = "")
    cat(" Copying files          \n")
    
    pbapply::pblapply(1:length(fll),function(i){file.copy(fll[i],fllt[i])})
    
    #cat("building vrt                \n")
    #v1<-vrt(fllt);nm<-names(v1)
    nm<-afftnames;names(nm)<-nm
    ### Extract values ---------------------
    cat("extracting                \n\n")
   # cat(system.time(mat2<-extract(v1,xy[,1:2]))[3]/60,"min.        \r")
    
  
    sfInit(T,30)
    sfExport("xy")
    sfExport("fllt")
    sfExport("afftnames")
    sfLibrary(terra)
    tmp<-sfClusterApplyLB(fllt,function(x){y<-rast(x);names(y)<-afftnames;z<-extract(y,xy[,1:2]);z<-z[!is.na(z[,2]),];z})
    mat3<-unique(do.call(rbind,tmp))
    mat3<-unique(merge(mat1,mat3,by = "ID"))
    
    mat4<-split(mat3,mat3$ID)
    mat4<-pbapply::pblapply(mat4,function(x){if(nrow(x)>1){rbind(sapply(x,function(x){if(is.factor(x)){return(x[1])}else{return(mean(x))}}))}else{return(x)}})
    mat4<-do.call(rbind,mat4)
    ## Correct Layers -----
    cat("\r Correcting layers                 \n\n")
    terraOptions(memfrac = .8)
    sfInit(T,19)
    sfLibrary(terra)
    s1<-sources(suppl)
    sfExport(list = c("mat4","b","s1","afftnames","fllt","tl","fll"))
    
    write("Correcting","D:/LocalNaip/progress.txt",append = T)
   # if(b == "ndvi"){tl<-tl[19]}
     # nb<-length(afftnames)
    tmp<-pbapply::pblapply(tl,function(tl1){
      sfExport("tl1")
      terraOptions(memfrac = .9/30)
      
      write(paste("Tile",tl1,":",which(tl==tl1),"out of",length(tl)),"D:/LocalNaip/Progress.txt",append = T)
      q<-sfLapply(afftnames,function(b2){
        gc()
        mat5<-mat4[,c(colnames(mat4)[2:5],b2)]
        colnames(mat5)[5]<-"y"
        
        mat5<-mat5[apply(mat5,1,function(x){!any(is.na(x))}),]
        suppl<-rast(s1)
       
        
        lm1<-lm(y~jdate+sl+imgyr+ jdate*sl ,data = mat5)#jdate*imgyr+imgyr+stateid+
        
        fll2<-list.files("C:/TEMP/temptiles",full.names = T)
        
        fll2<-rast(fll2[grepl(tl1,fll2)])
        names(fll2)<-afftnames
        fll2<-subset(fll2,b2)
        suppl3<-crop(suppl,fll2)
       # suppl3$imgyr<-as.factor(suppl3$imgyr)
        
        #suppl3<-trim(suppl3)
        fll2<-crop(fll2,suppl3)
        fll2<-rast(list(y = fll2,suppl3))
        
        outname<-paste(b,"_",b2,sep = "")
        outfile<-paste("C:/TEMP/RockyMountains/Corrected/",outname,"_",tl1,".tif",sep = "")
        
        names(fll2)[1]<-"y"
       predfun<-function(mod = lm1,newdata){
         df<-data.frame(newdata)
         newdata$imgyr<-factor(newdata$imgyr,levels =c("19","20","21"))
         predict(lm1,newdata)
       }
        pred1<-predict(fll2,lm1,fun = predfun,na.rm = T)
        pred2<-fll2$y - pred1
        names(pred2)<-outname
        
        pred2<-writeRaster(pred2,filename = outfile,overwrite = T)#
        gc()
        sources(pred2)
      })
    })
    #rm("v1")
    file.remove(fllt)
    sfStop()
  }  
    
#### Unify tiles to a vrt -----------------------

